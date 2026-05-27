"""
Contract test: every scoring-backend install script must verify the backend
package actually loads before writing the snakemake marker file. The marker
file is what `run_*` rules consume as their dependency edge — if it gets
written when the install silently failed (e.g. renv::restore() reported
errors but didn't raise), the scoring run downstream fails with a cryptic
"could not find function" error instead of a clear "install failed" one.

The BLAS/LAPACK devcontainer gap (2026-05-27) was hidden by exactly this
class of bug: install_rosace.R lacked the check, so a missing system lib
that broke `urca`/`robustbase`/etc. left a "Rosace installed." marker
behind. install_rosace_aa.R's check was what surfaced the real failure.

This test is intentionally coarse (regex against the file contents): each
script is short enough that a stronger AST/exec-level test would be more
fragile than the contract itself. If you refactor these scripts, update
this test in lockstep.
"""

import re
from pathlib import Path

import pytest

SCRIPTS_DIR = (
    Path(__file__).resolve().parents[2] / "workflow" / "rules" / "scripts"
)

# Each entry: install script -> package name that must appear in the
# require(...) check + stop(...) sentence.
INSTALL_SCRIPTS = {
    "install_rosace.R": "rosace",
    "install_lilace.R": "lilace",
    "install_rosace_aa.R": "rosaceAA",
}


@pytest.mark.parametrize("script_name,pkg", list(INSTALL_SCRIPTS.items()))
def test_install_script_verifies_package_loads(script_name, pkg):
    src = (SCRIPTS_DIR / script_name).read_text()
    # The two-line shape we lock in:
    #   if (!require(<pkg>)) {
    #     stop("<message>")
    #   }
    # Tolerates whitespace + reasonably-formatted message strings.
    pattern = re.compile(
        rf"if\s*\(\s*!\s*require\s*\(\s*{re.escape(pkg)}\s*\)\s*\)\s*\{{"
        rf"[^}}]*?stop\s*\(",
        re.DOTALL,
    )
    assert pattern.search(src), (
        f"{script_name}: missing `if (!require({pkg})) {{ stop(...) }}` "
        f"check. The marker file gets written unconditionally at the bottom "
        f"of main(); without this gate, a silent renv::restore failure "
        f"would leave a 'installed' marker behind."
    )
