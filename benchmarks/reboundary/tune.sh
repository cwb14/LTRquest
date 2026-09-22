#!/usr/bin/env bash
# Tuning sweeps for ltrquest-reboundary (design spec section 9, "Tuning order").
# Run from the LTRquest run directory. Each stage writes reboundary_bench/tune<STAGE>/.
#   tune.sh 1                         method x references            (B1 n=500 + B2)
#   METHOD=.. REFS=.. tune.sh 2       min_identity x anchor_len x max_indel  (B1 n=500)
#   KNOBS="--set ..." tune.sh 3       credit                         (B2 cached + B1)
#   KNOBS="--set ..." tune.sh 4a      max_ratio x min_copies          (B2 cached)
#   KNOBS="--set ..." tune.sh 4b      qc_min_n x untested             (B2 cached)
# Finished cells are skipped, so a stage can be re-run after an interruption.
set -euo pipefail
STAGE="$1"
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
    for m in ${METHODS:-consensus subfamily nearest}; do
      for r in modal tsd; do
        n="${m}_${r}"
        done_or "$OUT/b1_$n" || bench b1 --out "$OUT/b1_$n" --n 500 --set method=$m --set references=$r
        done_or "$OUT/b2_$n" || bench b2 --out "$OUT/b2_$n" --cache "$OUT/cache_$n.pkl" \
          --set method=$m --set references=$r
      done
    done ;;
  2)
    : "${METHOD:?set METHOD to the stage-1 winner}" "${REFS:?set REFS to the stage-1 winner}"
    for id in 0.7 0.8; do for al in 30 50 80; do for mi in 1000 5000 20000; do
      n="id${id}_al${al}_mi${mi}"
      done_or "$OUT/b1_$n" || bench b1 --out "$OUT/b1_$n" --n 500 --set method=$METHOD \
        --set references=$REFS --set min_identity=$id --set anchor_len=$al --set max_indel=$mi
    done; done; done ;;
  3)
    : "${KNOBS:?set KNOBS to the stage-2 winner as --set pairs}"
    for c in 0 50 200 1000 model; do
      done_or "$OUT/b2_credit_$c" || bench b2 --out "$OUT/b2_credit_$c" --cache "$OUT/cache.pkl" \
        $KNOBS --set credit=$c
      done_or "$OUT/b1_credit_$c" || bench b1 --out "$OUT/b1_credit_$c" --n 500 $KNOBS --set credit=$c
    done ;;
  4a)
    : "${KNOBS:?set KNOBS to the stage-3 winner as --set pairs}"
    # min_copies 5 first: its cached family phase covers the 10 and 20 runs as well.
    for mr in 1.10 1.15 1.25; do for mc in 5 10 20; do
      n="ratio${mr}_min${mc}"
      done_or "$OUT/b2_$n" || bench b2 --out "$OUT/b2_$n" --cache "$OUT/cache_ratio${mr}.pkl" \
        $KNOBS --set max_ratio=$mr --set min_copies=$mc
    done; done ;;
  4b)
    : "${KNOBS:?set KNOBS to the stage-4a winner as --set pairs}"
    for qn in 3 5 10; do for u in accept skip; do
      n="n${qn}_${u}"
      done_or "$OUT/b2_$n" || bench b2 --out "$OUT/b2_$n" --cache "$OUT/cache.pkl" \
        $KNOBS --set qc_min_n=$qn --set untested=$u
    done; done ;;
  *) echo "usage: tune.sh 1|2|3|4a|4b" >&2; exit 2 ;;
esac
python "$REPO/benchmarks/reboundary/bench.py" collect --out "$OUT"
echo "STAGE_${STAGE}_DONE"
