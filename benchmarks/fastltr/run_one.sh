#!/usr/bin/env bash
# run_one.sh GENOME CONFIG THREADS : one LTRquest benchmark run in runs/GENOME_CONFIG/
#   GENOME  athal | rice | chm13 | sim (PrinTE) | spike (spike_in.py)
#           (inputs/GENOME.fa; proteins inputs/GENOME.pep[.gz] if present)
#   CONFIG  U  (no cap, FastLTR off)   T  (cap 30 s, off)   TR (cap 30 s, rescue)
#           TA (cap 30 s, all)        UA (no cap, all)
# Plants: default rounds; chm13: --max-rounds 1 (as in the earlier CHM13 runs).
set -euo pipefail
B=${FASTLTR_BENCH:?set FASTLTR_BENCH to the benchmark directory (inputs/, regions/, runs/ ...)}
g=${1:?GENOME}; cfg=${2:?CONFIG}; t=${3:?THREADS}
# run inside the LTRquest environment (ltrquest, blastn, makeblastdb on PATH)
export LTRQUEST_TOOLS_DIR=${LTRQUEST_TOOLS_DIR:-$B/tools}
# FastLTR-on configs run the code snapshot with the tandem-array guard (v2); the
# FastLTR-off configs run the installed snapshot, which is the same code with FastLTR off.
SRC=${LTRQUEST_SRC:-}
[[ -n "$SRC" ]] && export PYTHONPATH="$SRC"
case $cfg in
  U)  opt=(--fastltr off) ;;
  T)  opt=(--fastltr off    --detect-args "ltr-timeout=30") ;;
  TR) opt=(--fastltr rescue --detect-args "ltr-timeout=30") ;;
  TA) opt=(--fastltr all    --detect-args "ltr-timeout=30") ;;
  UA) opt=(--fastltr all) ;;
  *) echo "unknown config $cfg" >&2; exit 2 ;;
esac
case $g in
  athal|sim|spike) mu=7e-9 ;; rice) mu=1.3e-8 ;;
  # CHM13 is scored on each run's first FP attempt, which LTRquest deletes on completion unless kept
  chm13) mu=2.2e-9; opt+=(--max-rounds 1 --dev-keep-fp-rounds) ;;
  *) echo "unknown genome $g" >&2; exit 2 ;;
esac
prot=()
for p in $B/inputs/$g.pep $B/inputs/$g.pep.gz; do [[ -e $p ]] && prot=(--proteins "$p"); done
d=$B/runs/${g}_${cfg}; mkdir -p "$d"; cd "$d"
ln -sfn "$B/inputs/$g.fa" "$g.fa"
echo "[run_one] $g $cfg start $(date -Is) threads=$t opts: ${opt[*]} ${prot[*]} code: $(python -c 'import ltrquest.fastltr as f; print(f.__file__)')"
rc=0
/usr/bin/time -v -o time.txt ltrquest --genome "$g.fa" "${prot[@]}" --threads "$t" \
    --mutation-rate "$mu" --no-plots "${opt[@]}" > run.log 2>&1 || rc=$?
echo "[run_one] $g $cfg done $(date -Is) rc=$rc"
