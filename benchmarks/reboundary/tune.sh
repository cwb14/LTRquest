#!/usr/bin/env bash
# Tuning sweeps for ltrquest-reboundary (template re-boundarying).
# Run from the LTRquest run directory (made without --reboundary). Each stage writes
# reboundary_bench/tune<STAGE>/; every cell runs B1 (n=500) and B2 with the same settings.
#   tune.sh 1                         n_templates x min_support
#   KNOBS="--set ..." tune.sh 2       min_identity                   (KNOBS: the stage-1 winner)
#   KNOBS="--set ..." tune.sh 3       max_trim                       (KNOBS: stages 1-2 winners)
# --tools-dir must hold a Kmer2LTR that ltrquest-reboundary accepts (see CHANGELOG.md).
# Finished cells are skipped, so a stage can be re-run after an interruption.
set -euo pipefail
STAGE="${1:-}"
REPO=/data2/chris/poa_LTR/LTRquest
export PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH
export PYTHONPATH=$REPO/src
P=(annua.nuclear_LTRs chaixii.nuclear_LTRs infirma.nuclear_LTRs supina.nuclear_LTRs)
G=(annua.nuclear.fa chaixii.nuclear.fa infirma.nuclear.fa supina.nuclear.fa)
OUT=reboundary_bench/tune${STAGE}
mkdir -p "$OUT"
bench() {
  python "$REPO/benchmarks/reboundary/bench.py" "$@" --run . --prefix "${P[@]}" \
    --genome "${G[@]}" --tools-dir /data2/chris/poa_LTR/ltrquest_tools --threads "${THREADS:-32}"
}
done_or() { [[ -s "$1/summary.json" ]]; }
case "$STAGE" in
  1)
    for k in 3 5 10; do for s in 1 2 3; do
      n="k${k}_s${s}"
      done_or "$OUT/b1_$n" || bench b1 --out "$OUT/b1_$n" --n 500 \
        --set n_templates=$k --set min_support=$s
      done_or "$OUT/b2_$n" || bench b2 --out "$OUT/b2_$n" \
        --set n_templates=$k --set min_support=$s
    done; done ;;
  2)
    : "${KNOBS:?set KNOBS to the stage-1 winner as --set pairs}"
    for id in 0.7 0.8 0.9; do
      n="id${id}"
      done_or "$OUT/b1_$n" || bench b1 --out "$OUT/b1_$n" --n 500 $KNOBS --set min_identity=$id
      done_or "$OUT/b2_$n" || bench b2 --out "$OUT/b2_$n" $KNOBS --set min_identity=$id
    done ;;
  3)
    : "${KNOBS:?set KNOBS to the stage-1 and stage-2 winners as --set pairs}"
    for mt in 0 5 10 20; do
      n="trim${mt}"
      done_or "$OUT/b1_$n" || bench b1 --out "$OUT/b1_$n" --n 500 $KNOBS --set max_trim=$mt
      done_or "$OUT/b2_$n" || bench b2 --out "$OUT/b2_$n" $KNOBS --set max_trim=$mt
    done ;;
  *) echo "usage: tune.sh 1|2|3" >&2; exit 2 ;;
esac
python "$REPO/benchmarks/reboundary/bench.py" collect --out "$OUT"
echo "STAGE_${STAGE}_DONE"
