#!/usr/bin/env python3
"""Benchmark BBMap vs minimap2 on a paired-end FASTQ + reference.

Measures wall time, peak RSS, BAM record count, and output BAM size for both
aligners using the same flag set the dumpling snakemake pipeline uses. Built
to compare aligners on production-scale fixtures (e.g. lmna.bb.r{1,2}.fastq.gz
~8GB pair) without having to spin up a full snakemake run.

Index building happens once per aligner per (reference, k-mer) combo and is
excluded from the timing — only the per-sample mapping work is measured.
Re-runs reuse cached indices under --workdir.

Usage:
    python tools/bench_aligners.py \\
        --r1 lmna.bb.r1.fastq.gz \\
        --r2 lmna.bb.r2.fastq.gz \\
        --ref references/lmna_ref.fasta \\
        --threads 16 --mem 16

    # Skip an aligner if needed:
    python tools/bench_aligners.py ... --skip-bbmap
    python tools/bench_aligners.py ... --skip-minimap2

Run from a shell where `dumpling_env` is active so bbmap.sh, minimap2 and
samtools are all on PATH.

What it measures:
- Wall: python's time.perf_counter() around the full shell pipeline
  (aligner | samtools view -b -o file.bam).
- Peak RSS: polls /proc/<pid>/status VmRSS every 100ms for the duration of
  the run. Catches the alignment process's peak; brief spikes shorter than
  100ms may be missed, but for multi-second alignment jobs this is fine.
- BAM record count + BAM size: read after each run.

Caveats:
- On tiny DMS references (~1-10 kb), BBMap's perfect-hash index allocates
  GB-scale memory regardless of reference size while minimap2's minimizer
  index scales with the reference. So the wall/RSS ratios are dramatic on
  small refs and shrink on whole-genome refs. The ratios reported here
  are valid for the input you provide — do not extrapolate across
  reference-size regimes.
- Both aligners' BAM record counts will differ slightly because they make
  different soft-clip/multi-mapping decisions. This script reports the
  delta but doesn't gate on it. For biological concordance, run both
  through GATK ASM and compare scores (Pearson r > 0.99 is the gate the
  audit set; see tools/ for related comparison scripts if added).
"""

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path


def _read_vmrss_kb(pid: int) -> int:
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmRSS:"):
                    return int(line.split()[1])
    except (FileNotFoundError, ProcessLookupError):
        pass
    return 0


def _descendants(root_pid: int) -> list[int]:
    """All descendants of `root_pid` via /proc/<pid>/task/<pid>/children,
    including the root. Misses processes that exit between iterations."""
    out = []
    queue = [root_pid]
    while queue:
        pid = queue.pop()
        out.append(pid)
        try:
            with open(f"/proc/{pid}/task/{pid}/children") as f:
                queue.extend(int(c) for c in f.read().split())
        except (FileNotFoundError, ProcessLookupError):
            pass
    return out


def poll_peak_rss_kb(root_pid: int, stop: threading.Event, interval: float = 0.1) -> int:
    """Poll the SUM of VmRSS across `root_pid` and all of its descendants.
    Tracks the maximum sum seen. Bash pipelines spawn the aligner and
    samtools as siblings under the bash subprocess, so this catches both."""
    peak = 0
    while not stop.is_set():
        current = sum(_read_vmrss_kb(p) for p in _descendants(root_pid))
        if current > peak:
            peak = current
        if current == 0 and not stop.is_set():
            # Root and all descendants gone; nothing more to measure.
            break
        time.sleep(interval)
    return peak


def run_with_measurement(cmd: str, label: str, log_path: Path) -> tuple[int, float, int]:
    """Run a shell command, measure wall and peak RSS of the immediate child.
    Returns (exit_code, wall_seconds, peak_rss_kb)."""
    print(f"\n=== {label}: starting ===", flush=True)
    print(f"  cmd: {cmd}", flush=True)
    print(f"  log: {log_path}", flush=True)
    start = time.perf_counter()
    with open(log_path, "wb") as log_file:
        proc = subprocess.Popen(
            cmd, shell=True, executable="/bin/bash",
            stdout=log_file, stderr=subprocess.STDOUT,
        )
        stop = threading.Event()
        peak_holder = [0]

        def poll():
            peak_holder[0] = poll_peak_rss_kb(proc.pid, stop)

        t = threading.Thread(target=poll, daemon=True)
        t.start()
        rc = proc.wait()
        stop.set()
        t.join()
    wall = time.perf_counter() - start
    print(f"  exit: {rc}, wall: {wall:.2f}s, peak_rss: {peak_holder[0] / 1024:.0f} MB", flush=True)
    return rc, wall, peak_holder[0]


def build_bbmap_index(ref: Path, index_dir: Path, mem_gb: int, threads: int, log: Path) -> None:
    """Build BBMap index (untimed, one-time cost)."""
    if (index_dir / "ref" / "genome" / "1" / "summary.txt").exists():
        print(f"=== BBMap index already at {index_dir}, skipping build ===", flush=True)
        return
    print(f"=== Building BBMap index → {index_dir} ===", flush=True)
    cmd = (
        f"bbmap.sh -Xmx{mem_gb}g ref={ref} path={index_dir} build=1 k=15 "
        f"32bit=t rebuild=t 2>&1 | tee {log}"
    )
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


def build_minimap2_index(ref: Path, index_file: Path, threads: int, log: Path) -> None:
    """Build minimap2 short-read index (untimed)."""
    if index_file.exists():
        print(f"=== minimap2 index already at {index_file}, skipping build ===", flush=True)
        return
    print(f"=== Building minimap2 index → {index_file} ===", flush=True)
    index_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = f"minimap2 -x sr -t {threads} -d {index_file} {ref} 2>&1 | tee {log}"
    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


def count_bam_records(bam: Path) -> int:
    out = subprocess.check_output(["samtools", "view", "-c", str(bam)])
    return int(out.strip())


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--r1", required=True, help="paired R1 FASTQ (.fastq.gz)")
    p.add_argument("--r2", required=True, help="paired R2 FASTQ (.fastq.gz)")
    p.add_argument("--ref", required=True, help="reference FASTA")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--mem", type=int, default=16, help="BBMap -Xmx in GB")
    p.add_argument("--workdir", default=None,
                   help="Where to keep indices + output BAMs. Default: a tmpdir under /tmp. "
                        "Reused across runs for cached indices.")
    p.add_argument("--skip-bbmap", action="store_true")
    p.add_argument("--skip-minimap2", action="store_true")
    p.add_argument("--cleanup", action="store_true",
                   help="Delete workdir after run (default: keep so you can inspect BAMs/logs)")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.skip_bbmap and args.skip_minimap2:
        print("Both aligners skipped. Nothing to do.", file=sys.stderr)
        return 2

    r1 = Path(args.r1).resolve()
    r2 = Path(args.r2).resolve()
    ref = Path(args.ref).resolve()
    for path in (r1, r2, ref):
        if not path.exists():
            print(f"missing: {path}", file=sys.stderr)
            return 1

    work = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="bench_aligners_"))
    work.mkdir(parents=True, exist_ok=True)
    print(f"Workdir: {work}")

    # Key indices on reference content hash so a changed reference produces
    # a fresh index by construction (mirrors how the snakemake pipeline does it).
    ref_digest = hashlib.sha256(ref.read_bytes()).hexdigest()[:12]
    bbmap_index = work / f"bbmap_{ref_digest}"
    mm2_index = work / f"mm2_{ref_digest}.mmi"

    results = {}

    if not args.skip_bbmap:
        build_bbmap_index(ref, bbmap_index, args.mem, args.threads, work / "bbmap_index.log")
        bbmap_bam = work / "bbmap.bam"
        bbmap_cmd = (
            f"bbmap.sh -Xmx{args.mem}g "
            f"in1={r1} in2={r2} "
            f"sam=1.3 32bit=t path={bbmap_index} build=1 "
            f"outm=stdout.sam k=15 t={args.threads} nzo=true "
            f"| samtools view -b -o {bbmap_bam} -"
        )
        rc, wall, rss_kb = run_with_measurement(bbmap_cmd, "BBMap", work / "bbmap.log")
        if rc != 0:
            print("BBMap failed, check log", file=sys.stderr)
            return rc
        results["BBMap"] = {
            "wall_s": wall,
            "peak_rss_mb": rss_kb / 1024,
            "bam_records": count_bam_records(bbmap_bam),
            "bam_mb": bbmap_bam.stat().st_size / (2 ** 20),
        }

    if not args.skip_minimap2:
        build_minimap2_index(ref, mm2_index, args.threads, work / "mm2_index.log")
        mm2_bam = work / "mm2.bam"
        mm2_cmd = (
            f"minimap2 -ax sr --MD -t {args.threads} {mm2_index} {r1} {r2} "
            f"| samtools view -b -o {mm2_bam} -"
        )
        rc, wall, rss_kb = run_with_measurement(mm2_cmd, "minimap2", work / "mm2.log")
        if rc != 0:
            print("minimap2 failed, check log", file=sys.stderr)
            return rc
        results["minimap2"] = {
            "wall_s": wall,
            "peak_rss_mb": rss_kb / 1024,
            "bam_records": count_bam_records(mm2_bam),
            "bam_mb": mm2_bam.stat().st_size / (2 ** 20),
        }

    print("\n=== Results ===")
    cols = list(results)
    print(f"{'metric':<22}{'  '.join(f'{c:>14}' for c in cols)}")
    for metric, fmt in [
        ("wall_s",       "{:>14.2f}"),
        ("peak_rss_mb",  "{:>14.0f}"),
        ("bam_records",  "{:>14d}"),
        ("bam_mb",       "{:>14.1f}"),
    ]:
        cells = "  ".join(fmt.format(results[c][metric]) for c in cols)
        print(f"{metric:<22}{cells}")

    if len(cols) == 2 and "BBMap" in results and "minimap2" in results:
        bb = results["BBMap"]; mm = results["minimap2"]
        print(f"\nRatios (BBMap / minimap2):")
        print(f"  wall:    {bb['wall_s'] / mm['wall_s']:.1f}x faster with minimap2")
        if mm["peak_rss_mb"] > 0:
            print(f"  peak RSS: {bb['peak_rss_mb'] / mm['peak_rss_mb']:.1f}x less with minimap2")
        rec_delta_pct = 100 * (mm["bam_records"] - bb["bam_records"]) / bb["bam_records"]
        print(f"  records: {rec_delta_pct:+.2f}% (minimap2 vs BBMap)")
        print(f"\nFor concordance, run both BAMs through GATK ASM and compare Pearson r")
        print(f"on the resulting variantCounts (audit gate: > 0.99 on high-count variants).")

    if args.cleanup:
        shutil.rmtree(work)
        print(f"\nCleaned up workdir.")
    else:
        print(f"\nWorkdir preserved at {work} for inspection. Pass --cleanup to remove.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
