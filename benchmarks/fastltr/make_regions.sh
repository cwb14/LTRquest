#!/usr/bin/env bash
# make_regions.sh : the region BEDs the analysis stratifies and counts by, into regions/
#   tools/FasTAN_stock  stock FasTAN (FASTAN 843d0ed, `make FasTAN`); tools/ONEview (ONEcode)
#   inputs/  athal.fa rice.fa chm13.fa; chm13v2.0_censat_v2.1.bed (T2T-CHM13 censat annotation);
#            GCA_009914755.4_assembly_report.txt (UCSC names -> the assembly's GenBank accessions);
#            rice.cenh3.bed (CENH3-defined centromeres of T2T-NIP); chm13.active_hor.bed (the
#            active alpha-satellite HOR array of each chromosome, censat v2.1)
# needs samtools when inputs/<genome>.fa.fai is missing
set -euo pipefail
B=${FASTLTR_BENCH:?set FASTLTR_BENCH to the benchmark directory (inputs/, regions/, runs/ ...)}
S=$(cd "$(dirname "$0")" && pwd)
mkdir -p $B/regions && cd $B/regions

for g in athal rice chm13; do
  if [[ -s $B/inputs/$g.fa.fai ]]; then cp $B/inputs/$g.fa.fai $g.fa.fai
  else samtools faidx $B/inputs/$g.fa -o $g.fa.fai; fi
done

# plants: tandem arrays >= 10 kb from a whole-genome FasTAN run; class by unit (micro < 10 bp,
# satellite < 1 kb, long >= 1 kb), median unit in column 5
for g in athal rice; do
  [[ -s $g.1aln ]] || $B/tools/FasTAN_stock -T8 -a -o$g $B/inputs/$g.fa
  python3 $S/fastan_arrays.py $B/tools/ONEview $g.1aln --min-len 10000 > $g.arrays10k.bed
done
# A. thaliana: CEN180 centromeres (170-185 bp unit, arrays >= 50 kb) and the two NORs (8-12 kb
# units in the first 6 Mb of a chromosome), each spanned first to last array per chromosome
span='{if(!($1 in s) || $2<s[$1]) s[$1]=$2; if($3>e[$1]) e[$1]=$3} END{for(c in s) print c"\t"s[c]"\t"e[c]}'
awk '$4=="satellite" && $5>=170 && $5<=185 && $3-$2>=50000' athal.arrays10k.bed | awk "$span" | sort > athal.cen180.bed
awk '$4=="long" && $5>=8000 && $5<=12000 && $3<6000000' athal.arrays10k.bed | awk "$span" | sort > athal.nor.bed
cp $B/inputs/rice.cenh3.bed rice.cenh3.bed

# CHM13: censat arrays >= 10 kb (rDNA long, the rest satellite; ct transition spans dropped)
python3 - $B/inputs/GCA_009914755.4_assembly_report.txt $B/inputs/chm13v2.0_censat_v2.1.bed <<'PY' \
    | sort -k1,1 -k2,2n > chm13.arrays10k.bed
import sys
acc = {}
for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) > 9:
        acc[f[9]] = f[4]
        acc['chr' + f[0]] = f[4]
for line in open(sys.argv[2]):
    f = line.split('\t')
    lab = f[3].split('_')[0].split('(')[0]
    s, e = int(f[1]), int(f[2])
    if lab == 'ct' or e - s < 10000 or f[0] not in acc:
        continue
    print(f"{acc[f[0]]}\t{s}\t{e}\t{'long' if lab == 'rDNA' else 'satellite'}\t0\t1")
PY
awk '$4=="long"' chm13.arrays10k.bed > chm13.rdna.bed
cp $B/inputs/chm13.active_hor.bed chm13.active_hor.bed
echo "done: $(ls $B/regions)"
