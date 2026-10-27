# FastLTR benchmarks

Does a FastLTR third caller recover what `--ltr-timeout` throws away, and what
does it cost elsewhere? How the caller works:
[docs/outputs.md §3](../../docs/outputs.md#3-usage) ("Capped chunks").

## What a cap costs

LTRharvest and LTR_FINDER write nothing until they finish, so a call killed by
`--ltr-timeout` loses every candidate of its 500 kb chunk, not only the array's.
Two CHM13 chunks capped at 30 s (chrX 57.48-57.98 Mb, an alpha-satellite edge;
chr1 1.90-2.40 Mb, one 182 kb array of a 40-bp unit) gave 0 candidates from
either tool.

## Tests

| test | genome | truth / check |
|---|---|---|
| uncapped reference | *A. thaliana* T2T (142.5 Mb), rice T2T-NIP AGIS1.0 (390 Mb) | the uncapped run's own clean calls (`U`) |
| curated support | same | an element is supported when >= 80% of its bases align at >= 80% identity (>= 80 bp HSPs) to the EDTA curated LTR library (athrep, rice7.0.0) |
| RepeatMasker support | CHM13 v2.0 (3.1 Gb) | both LTRs >= 50% covered by RepeatMasker LTR-class records sharing a family name |
| planted elements | 22 Mb of *A. thaliana* CEN180 and 45S rDNA arrays with 417 library LTR-RTs inserted (5-bp TSDs, per-base divergence 0-4%; `spike_in.py`) | exact |
| simulated genome | PrinTE burn-in (35.7 Mb, 815 intact LTR-RTs, no arrays) | exact |

Configurations (`run_one.sh`): `U` no cap, FastLTR off; `T` 30 s cap, off;
`TR` 30 s cap, `--fastltr rescue`; `TA` 30 s cap, `--fastltr all`; `UA` no cap,
`--fastltr all`. Plants run the default rounds; CHM13 one round. Elements are
matched when each LTR overlaps its counterpart by >= 50% of the shorter. The
reference's elements are stratified by what the `T` run's round-1 caps cost
their location: `lost2` (every chunk holding the element was capped for both
callers), `lost1` (for one), `open`.

The caps are wall-clock limits, so the set of capped chunks shifts a little
with machine load; capped runs of one genome were launched together.

## Results

### What the cap costs, and what rescue restores (uncapped run as reference)

A 30 s cap fired on 60 of 302 *A. thaliana* round-1 chunks (28.7 Mb, 20% of the
genome: the five CEN180 centromeres and the NOR2/NOR4 45S rDNA arrays; 14 chunks
for both callers) and on 26 of 816 rice chunks (12.7 Mb, 3.3%, mostly CentO
centromeres; 9 for both). Recall of the uncapped run's clean elements, and of the
ones a curated LTR library supports:

| *A. thaliana* T2T | wall / CPU | clean | no cap fired | one caller capped | both capped | supported, both capped |
|---|---|---|---|---|---|---|
| U uncapped (reference) | 1.48 h / 14.3 h | 1,214 | | | | |
| T cap 30 s | 0.33 h / 2.0 h | 938 | 97.5% | 46.4% | 0.0% | 0 / 97 |
| **TR cap + rescue** | 0.37 h / 2.5 h | 1,110 | 97.7% | 51.2% | 50.3% | **66 / 97** |
| TA cap + all | 0.42 h / 2.7 h | 1,252 | 94.2% | 51.2% | 50.8% | 67 / 97 |
| UA uncapped + all | 1.14 h / 12.7 h | 1,403 | 96.6% | 95.2% | 95.8% | 96 / 97 |

| rice T2T-NIP | wall / CPU | clean | no cap fired | one caller capped | both capped | supported, both capped |
|---|---|---|---|---|---|---|
| U uncapped (reference) | 1.26 h / 7.3 h | 7,444 | | | | |
| T cap 30 s | 0.81 h / 5.1 h | 7,316 | 99.6% | 79.9% | 0.9% | 1 / 76 |
| **TR cap + rescue** | 0.70 h / 4.9 h | 7,477 | 99.6% | 86.6% | 72.0% | **65 / 76** |
| UA uncapped + all | 1.00 h / 7.6 h | 8,606 | 96.9% | 96.6% | 100% | 76 / 76 |

Wall times are on a shared machine (load 200-270); CPU hours compare better.
Elements a run adds that the reference lacks, and their curated support:
*A. thaliana* T 25 (8%), TR 92 (39%), TA 263 (21%), UA 235 (16%); rice T 41
(24%), TR 117 (64%, the reference's own rate), UA 1,357 in uncapped chunks alone
(36%, against 63%). Whole-run totals, curated-supported / unsupported:
*A. thaliana* U 541 / 673, T 422 / 516, **TR 535 / 575**, UA 571 / 832; rice
U 4,742 / 2,702, T 4,648 / 2,668, **TR 4,781 / 2,696**, UA 5,134 / 3,472.

### Inside centromeres

| | *A. thaliana* CEN180 spans: elements (curated-supported; ATHILA) | rice CENH3 domains: elements (supported; CRM) |
|---|---|---|
| U uncapped | 94 (56; 57) | 155 (116; 66) |
| T cap 30 s | 47 (32; 31) | 97 (76; 40) |
| **TR cap + rescue** | 76 (55; 55) | 150 (120; 67) |

In the *A. thaliana* NORs every run's calls are `LTR/unknown` with no curated
support (U 48, T 10, TR 15): the uncapped run's extra rDNA calls are array units
the tandem filter misses, not elements rescue fails to find.

### Planted elements in real arrays (exact truth)

417 library LTR-RTs planted into 22 Mb of *A. thaliana* CEN180 and 45S rDNA
arrays; 16 sit in two chunks every run loses (four `X` in the library elements
make GenomeTools reject them), leaving 401 reachable:

| | recovered (>= 90% reciprocal overlap) | in CEN180 | in rDNA / long-unit arrays | substitution rate >= 3% | clean calls |
|---|---|---|---|---|---|
| T cap 30 s | 103 (26%) | 87 / 190 | 16 / 227 | 12 / 93 | 182 |
| **TR cap + rescue** | **352 (88%)** | 180 / 190 | 172 / 227 | 69 / 93 | 520 |
| TA cap + all | 353 (88%) | 180 / 190 | 173 / 227 | 70 / 93 | 563 |
| U uncapped | 368 (92%) | 185 / 190 | 183 / 227 | 75 / 93 | 603 |

Wall time: T 0.15 h, TR 0.22 h, U 1.59 h.

### A genome without arrays (exact truth)

PrinTE burn-in genome, 35.7 Mb, 815 intact LTR-RTs, no cap, >= 90% reciprocal
overlap: U precision 0.973, recall 0.921, F1 0.946; UA (FastLTR on every chunk)
0.968, 0.935, 0.951 (13 elements only UA finds, 2 only U).

### The tandem-array guard

A first version without it rescued, in *A. thaliana*, 80 calls in long-unit and
rDNA arrays with 3% curated support: ~1 kb pairs of sub-unit repeats recurring
every unit. Uncapped LTRharvest/LTR_FINDER make few such calls (one such chunk:
2 + 107 candidates against FastLTR's 145 pairs), since LTRharvest keeps only the
best of overlapping predictions. Dropping a FastLTR pair whose copy recurs one
spacing further out took the rescued additions over T from 349 (126 supported /
223 not) to 201 (119 / 82), and planted-element recall from 351 to 352. It costs
real elements next to a same-family copy: 66 instead of 72 of the 97 supported
fully-capped elements. Set on *A. thaliana*; on rice, not used to set it, it
flags 1.4% of the FastLTR pairs behind high-confidence calls (1.0% on *A.
thaliana*).

### Why `rescue`, not `all`, is the default

Fixed before the runs: `all` would be the default only if, in chunks no cap fired
on, the elements `UA` adds over `U` had >= 0.8x the reference's own curated
support there, in both plants. Rice 0.36 vs 0.63, *A. thaliana* 0.13 vs 0.44:
both fail. `all` raises recall (rice +392 supported elements) at a larger cost in
unsupported ones (+770).

### CHM13 (human T2T, one round)

An uncapped CHM13 run does not finish in reasonable time (single chunks take
8-40+ min), so `T` is the reference here and RepeatMasker the independent
check. On a shared machine the 30 s cap cut LTRharvest short on 678 of 6,571
round-1 chunks and LTR_FINDER on 904: far beyond the satellite arrays, since a
loaded machine pushes ordinary chunks past 30 s too. With rescue, FastLTR
searched 1,067 chunks and dropped 23,680 array-unit pairs. Both runs then
exceeded LTRquest's false-positive threshold (26% of elements in FP families)
and re-ran once on an FP-masked genome (FP fraction then 3.6% for `T`, 1.5% for
`TR`). `run_one.sh` keeps both attempts for CHM13 (`--dev-keep-fp-rounds`), so
each run is scored on its first attempt and on its final output.

| | wall / CPU | first attempt: clean | RepeatMasker-supported | final: clean | RepeatMasker-supported |
|---|---|---|---|---|---|
| T cap 30 s | 5.87 h / 108 h | 8,981 | 61% | 15,488 | 57% |
| **TR cap + rescue** | 5.92 h / 110 h | 10,181 | 58% | 16,797 | 54% |

On the first attempt TR keeps 98.7% of T's elements (99.2% of the
RepeatMasker-supported ones) and adds 1,321: 911 outside tandem arrays (41%
supported, 60 high-confidence), 25 in satellite arrays (68%), 10 near telomeres
(60%), and 375 in the 45S rDNA arrays of the five acrocentric chromosomes (all
`LTR/unknown`, low confidence, 20%). Most of those last are copies of a
structure inside the 43 kb rDNA unit: the unit is longer than FastLTR's 20 kb
reach and than the tandem filter's flank test, so neither recognizes the array.
LTRharvest and LTR_FINDER make the same calls when they get through an rDNA
chunk (LTRquest 1.0.1 at a 300 s cap: 86 such calls), so this is a gap in the
tandem filter for long-unit arrays rather than in FastLTR.

LTRquest's FP correction removes most of them: TR's final output keeps 146 rDNA
calls, and the 101 that RepeatMasker supports are copies of a provirus the rDNA
unit carries, a HERVIP10FH internal region between MER39B LTRs (5.4 kb; the
MER39B pairs recur every 38-41 kb). `T`, whose rDNA chunks were capped, has
none. In the final outputs TR keeps 98.1% of T's elements (98.6% of the
supported ones) and adds 1,608. Of these, 369 lie in chunks both runs' last
attempts capped (46% supported), and 1,038 in chunks neither capped (15%), which
FastLTR did not search on that attempt: those follow from the two runs masking
different FP families before their second attempt. A replicate `T` run, launched earlier the same day with
the same settings, differs from this `T` the same way: 1,193 elements added and
190 lost, 1,025 of the additions in chunks neither run capped (3% supported).
About 1,000 elements of a capped CHM13 final output thus change from run to
run; what rescue itself adds is the capped-chunk and array calls.

## Rerun

```bash
export FASTLTR_BENCH=/path/to/bench   # inputs/ (see make_regions.sh), tools/, external/
export EDTA_DB=/path/to/EDTA/share/EDTA/database
bash make_regions.sh                  # array and centromere BEDs into $FASTLTR_BENCH/regions
bash run_one.sh athal TR 16           # one configuration; see the header for all
bash run_analysis.sh                  # the tables above, into $FASTLTR_BENCH/analysis
```

`inputs/` holds `<genome>.fa` (and proteins, `<genome>.pep[.gz]`) for `athal`,
`rice`, `chm13`, `sim` and `spike`, the PrinTE truth `sim.truth.bed`, and the
region sources `make_regions.sh` names; `external/` the CHM13 v2.0 RepeatMasker
BED (`chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed`); `tools/` a stock FasTAN,
`ONEview`, and LTRquest's tools directory (`LTRQUEST_TOOLS_DIR`). The planted
genome (`inputs/spike.fa`, truth `spike/athal_spike.truth.bed`) comes from
`spike_in.py` with PrinTE's *A. thaliana* library
(`data/athrep.updated.fasta_062024.ltr.full.gz` in
[PrinTE](https://github.com/cwb14/PrinTE)):

```bash
B=$FASTLTR_BENCH
python spike_in.py --genome $B/inputs/athal.fa --arrays $B/regions/athal.arrays10k.bed \
    --lib athrep.updated.fasta_062024.ltr.full.gz \
    --regions Athaliana_T01:14250000-18075000 Athaliana_T02:0-5250000 \
      Athaliana_T02:9025000-11900000 Athaliana_T03:13300000-17125000 Athaliana_T04:0-3825000 \
    --out $B/spike/athal_spike --seed 11
ln -s $B/spike/athal_spike.fa $B/inputs/spike.fa
```

The tandem-array guard numbers come from runs of the build before the guard
existed, and the CHM13 replicate from an earlier `T` run; `run_analysis.sh`
reproduces neither.
