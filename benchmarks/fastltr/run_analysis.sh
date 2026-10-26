#!/usr/bin/env bash
# run_analysis.sh : score every finished benchmark run; tables land in analysis/
set -euo pipefail
B=${FASTLTR_BENCH:?set FASTLTR_BENCH to the benchmark directory (inputs/, regions/, runs/ ...)}
S=$(cd "$(dirname "$0")" && pwd)
# run inside the LTRquest environment (ltrquest, blastn, makeblastdb on PATH), after make_regions.sh
DB=${EDTA_DB:?set EDTA_DB to EDTA share/EDTA/database (curated rice7.0.0 and athrep libraries)}
mkdir -p $B/analysis && cd $B/analysis

# cost and yield of every run
python3 $S/run_stats.py $(ls -d $B/runs/*/ | sed 's|/$||') > run_stats.tsv

# plants: reference = uncapped (U); strata from the capped run T; every configuration that was run
for g in athal rice; do
  lib=$DB/athrep.updated.nonredun.fasta; [[ $g == rice ]] && lib=$DB/rice7.0.0.liban.LTR
  cfg=(); for c in U T TR TA UA; do [[ -d $B/runs/${g}_$c ]] && cfg+=($c); done
  python3 $S/capped_chunks.py $B/runs/${g}_T/run.log > ${g}_T.capped.bed
  python3 $S/curated_support.py --lib "$lib" --threads 16 --out $g.support.tsv \
      $(for c in "${cfg[@]}"; do echo $B/runs/${g}_$c/${g}_LTRs_depth*_clean_ltr.fa; done) 2> /dev/null
  python3 $S/analyze.py --genome $g --runs $B/runs --configs "${cfg[@]}" --ref U --capped-from T \
      --capped-bed ${g}_T.capped.bed --fai $B/regions/$g.fa.fai --arrays $B/regions/$g.arrays10k.bed \
      --support $g.support.tsv --out ${g}_vs_U > /dev/null
done
# centromeres and NORs
for r in athal.cen180 athal.nor rice.cenh3; do
  g=${r%%.*}
  python3 $S/region_counts.py --bed $B/regions/$r.bed --support $g.support.tsv $B/runs/${g}_{U,T,TR} > $r.txt
done

# CHM13: no uncapped run is feasible; reference = T, support = RepeatMasker LTR pairs. Both runs
# exceeded the FP threshold and re-ran on an FP-masked genome (TR's rerun was stopped), so each
# is scored on its complete first attempt, <prefix>_FPstaging/round1 (kept by --dev-keep-fp-rounds).
att=chm13_LTRs_FPstaging/round1
if [[ -d $B/runs/chm13_T/$att && -d $B/runs/chm13_TR/$att ]]; then
  mkdir -p chm13_att1
  for c in T TR; do ln -sfn $B/runs/chm13_$c/$att chm13_att1/chm13_$c; done
  awk '/re-detecting next attempt/{exit} 1' $B/runs/chm13_T/run.log > chm13_T.att1.log
  python3 $S/capped_chunks.py chm13_T.att1.log > chm13_T.capped.bed
  python3 $S/rm_support.py --rm $B/external/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed \
      --report $B/inputs/GCA_009914755.4_assembly_report.txt \
      --out chm13.support.tsv chm13_att1/chm13_{T,TR}/chm13_LTRs_depth*_clean_ltr.tsv
  python3 $S/analyze.py --genome chm13 --runs chm13_att1 --configs T TR --ref T --capped-from T \
      --capped-bed chm13_T.capped.bed --fai $B/regions/chm13.fa.fai --arrays $B/regions/chm13.arrays10k.bed \
      --support chm13.support.tsv --out chm13_vs_T > /dev/null
  python3 $S/additions.py --genome chm13 --ref chm13_att1/chm13_T --test chm13_att1/chm13_TR \
      --fai $B/regions/chm13.fa.fai --arrays $B/regions/chm13.arrays10k.bed \
      --support chm13.support.tsv > chm13_TR_additions.txt
  grep -m1 "FastLTR searched" $B/runs/chm13_TR/run.log >> chm13_TR_additions.txt || true
else
  echo "CHM13 skipped: a run lacks its first attempt ($att); rerun with run_one.sh" >&2
fi

# exact truth: PrinTE genome (no arrays) and planted elements in real arrays
python3 $S/sim_eval.py $B/inputs/sim.truth.bed $B/runs/sim_U $B/runs/sim_UA > sim_eval.txt
python3 $S/spike_eval.py $B/spike/athal_spike.truth.bed $B/runs/spike_{U,T,TR,TA} \
    --failed-log $B/runs/spike_T/run.log > spike_eval.txt
echo "done: $(ls $B/analysis)"
