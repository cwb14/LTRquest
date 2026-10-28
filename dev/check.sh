#!/usr/bin/env bash
# Run, before pushing, the checks GitHub runs on every push (see AGENTS.md).
#
#   bash dev/check.sh
#   NXF_OLD=/path/to/nextflow-24.04.4 NXF_NEW=/path/to/nextflow-latest bash dev/check.sh
#
# Mirrors .github/workflows/ci.yml: ruff, shellcheck, the non-slow tests, and the
# Nextflow stub runs under the oldest supported and the newest Nextflow. Tools come
# from PATH; NXF_OLD / NXF_NEW name two nextflow launchers (else `nextflow` on
# PATH runs once), each run with the Java beside it, as in its conda env. The tests
# import this checkout, whatever ltrquest is installed. A check whose tool is
# missing is reported as skipped, never as passed. The container is not built here.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
repo=$PWD
export PYTHONPATH="${repo}/src${PYTHONPATH:+:$PYTHONPATH}"
failed=()
skipped=()

check() {               # check NAME CMD...: run CMD, record it if it fails
  local name=$1
  shift
  echo "== ${name}"
  if "$@"; then echo "   ok"; else failed+=("$name"); fi
}

have() {                # have TOOL NAME: TOOL is on PATH, or NAME is recorded as skipped
  command -v "$1" > /dev/null && return 0
  skipped+=("$2 (no $1 on PATH)")
  return 1
}

stub() {                # stub NEXTFLOW: CI's stub job, in a scratch directory
  local nf wd rc
  nf=$(command -v "$1") || { echo "   not found: $1"; return 1; }
  wd=$(mktemp -d)
  (
    cd "$wd" || exit 1
    PATH="$(dirname "$nf"):$PATH"
    export PATH NXF_ANSI_LOG=false
    {
      "$nf" run "$repo" -profile test -stub-run --outdir results \
        && test -s results/test_LTRs/test_LTRs_all_depth_LTR_cleaned.gff3 \
        && test -s results/test_LTRs/test_LTRs_depth0_clean_ltr.tsv \
        && test -s results/families/test_all_ltr.consensus_id0.75_cluster.tsv \
        && "$nf" run "$repo" -profile test -stub-run --max_rounds 1 --outdir results_r1 \
        && test -s results_r1/test_LTRs/test_LTRs_all_depth_LTR_cleaned.gff3
    } > stub.log 2>&1 || { tail -n 20 stub.log; exit 1; }
  )
  rc=$?
  rm -rf "$wd"
  return "$rc"
}

have ruff ruff && check "ruff" ruff check .
have shellcheck shellcheck && check "shellcheck" shellcheck --severity=warning \
  src/ltrquest/scripts/*.sh bin/entrypoint.sh bin/ltrquest-container dev/check.sh
check "pytest -m 'not slow'" python -m pytest -m "not slow" -q -p no:cacheprovider

nxf=()
for v in NXF_OLD NXF_NEW; do [[ -n "${!v:-}" ]] && nxf+=("${!v}"); done
if (( ${#nxf[@]} == 0 )) && have nextflow "nextflow stub"; then
  nxf=(nextflow)
  skipped+=("nextflow stub under a second version (set NXF_OLD and NXF_NEW)")
fi
for nf in ${nxf[@]+"${nxf[@]}"}; do check "nextflow stub (${nf})" stub "$nf"; done

echo
if (( ${#skipped[@]} )); then printf 'skipped: %s\n' "${skipped[@]}"; fi
if (( ${#failed[@]} )); then printf 'FAILED: %s\n' "${failed[@]}"; exit 1; fi
echo "every check that ran passed"
