# Contributing to dumpling

## Development environment

dumpling separates its **Python test/lint layer** from the **pipeline runtime**,
and they use different tooling:

- **Python test/lint layer** (pytest, the libraries the workflow scripts import,
  snakemake for DAG construction, ruff/snakefmt) → an isolated
  [uv](https://docs.astral.sh/uv/) virtual environment. This keeps dev
  dependencies out of your conda base env.
- **Pipeline runtime** (BBTools, GATK, minimap2, samtools, R + Rosace/Lilace,
  CmdStan, Enrich2) → conda envs / the prebuilt container. These are not
  pip-installable; see the README "Installation" and "prebuilt container"
  sections.

### Running the tests

The dev/test layer is declared in `pyproject.toml` (the `dev` dependency group)
and pinned in `uv.lock`. `uv sync` creates `.venv` from it:

```bash
uv sync                # create .venv from the dev group (gitignored)
uv run pytest          # full suite, incl. the snakemake dry-run tests
uvx ruff@0.8.0 check   # lint (pinned version; matches CI)
uvx ruff@0.8.0 format  # format
```

`uv run pytest` covers **unit logic + DAG construction (dry-run)** — not a full
pipeline execution, since the bioinformatics tools live in the runtime layer.
The `dev` group is the pure-Python layer: `pytest`, `pytest-mock`, `pandas`,
`biopython`, `regex`, `jsonschema`, `pyyaml`, `mavehgvs`, `snakemake`, plus
`ruff`/`snakefmt`/`pre-commit`.

> **Why uv, not conda base?** Without a project venv, `uv run` silently falls
> back to your conda base interpreter, so dev/test deps (and their transitive
> pins) pile up in `base` and collide across projects. A project-local `.venv`
> (with `uv.lock`) makes the test environment reproducible and isolated.

### Running the cosmos-export validation test

The cosmos test loads a generated CSV into the real
[cosmos](https://github.com/pimentellab/cosmos) package. It is skip-guarded, so
it skips unless cosmos is installed. Add the optional `cosmos` extra:

```bash
uv sync --extra cosmos
uv run pytest tests/integration/test_cosmos_integration.py
```

The extra is kept out of the default `dev` group only because it pulls a heavy
stack (numpy 2.x, arviz, xarray, scikit-learn) most runs don't need — not for
any version conflict (the `.venv` is isolated, so its numpy 2.x can't collide
with a numpy-1.x tool in another env). The `matplotlib<3.11` pin in the extra
works around an `arviz` import of `matplotlib.style.core`, removed in matplotlib
3.11.
