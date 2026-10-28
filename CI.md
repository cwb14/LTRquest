# Notes for developers and coding agents

## Before every push: run what CI runs

Each push to `main` runs the jobs in `.github/workflows/ci.yml` and
`docker.yml`, and one failure marks the commit "Some checks were not
successful". A job that is already red hides the next failure in it: from
24 Sep 2026 every push failed ruff, so shellcheck, the step after it, never ran;
the newest-Nextflow stub run failed from 12 Oct, the day Nextflow 26.10 shipped;
and the image build failed from 26 Oct, so no image was published after that.

    bash dev/check.sh

runs ruff, shellcheck, the non-slow tests and the Nextflow stub runs, and exits
non-zero if one fails. A check whose tool is missing is listed as skipped, not
passed. Point it at both Nextflow versions CI uses:

    mamba create -n ltrquest_nf -c conda-forge -c bioconda nextflow shellcheck  # newest
    mamba create -n ltrquest_nf2404 -c conda-forge -c bioconda nextflow=24.04.4  # oldest allowed
    NXF_OLD=<envs>/ltrquest_nf2404/bin/nextflow NXF_NEW=<envs>/ltrquest_nf/bin/nextflow \
        PATH=$PATH:<envs>/ltrquest_nf/bin bash dev/check.sh

Run it from an environment with the package's dev extras (`pip install -e
".[dev]"`) and the tools in `environment.yml`, or the tests that need blastn or
Kmer2LTR skip.

## The jobs, and what has turned each red

- **lint: ruff** (`ruff check .`; `ruff check --fix .` sorts imports). Covers
  `tests/` too, at a Python 3.9 floor; CI installs the newest ruff.
- **lint: shellcheck** (`--severity=warning` on `src/ltrquest/scripts/*.sh`,
  `bin/entrypoint.sh`, `bin/ltrquest-container`). Runs only once ruff passes.
- **pytest** on Python 3.9, 3.11 and 3.13 (`pytest -m "not slow"`). Locally,
  put `src/` first on `PYTHONPATH`, as `dev/check.sh` does: in an environment
  with an older LTRquest installed, pytest otherwise tests that copy.
- **nextflow stub** under 24.04.4 and `latest-stable`, which picks up each
  Nextflow release the day it ships, and releases keep tightening the syntax:
  26.10 refused `return` in a workflow body (`main.nf` now uses `exit 0`).
  `nextflow lint .` under the newest Nextflow lists the errors; its warnings
  (`Channel` for `channel`, `when:` blocks, implicit `it`) are the next ones.
- **pip install** builds the wheel and runs each entry point's `--help`. Not in
  `dev/check.sh`; packaging changes need it.
- **container: build and smoke-test**. Not in `dev/check.sh`; needs Docker
  (`docker build .`). Dockerfile `RUN` lines execute under `bash -o pipefail`
  (micromamba's `SHELL`), so `tool | grep -q x` fails when the tool exits
  non-zero even though grep matched: write `{ tool 2>&1 || true; } | grep -q x`.
  Every `ARG *_REF` must be a commit that exists on GitHub, not only locally.
- **Slow tests** (`pytest -m slow`) never run in CI. Run them after touching the
  driver or a stage's outputs; `LTRQUEST_TOOLS_DIR` set to a built tools
  directory skips the clones.

## After pushing

Check the Actions tab (or `gh run list --limit 4`) and fix a red check before
the next push.
