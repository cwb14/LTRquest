# Re-boundarying benchmarks

What `ltrquest-reboundary` is measured on, how to reproduce it, and the numbers
behind its defaults. Design: `docs/superpowers/specs/2026-09-22-reboundary-design.md` §9.

## The three benchmarks

- **B1 — planted obstacles (ground truth).** Elements whose called ends carry a
  Kmer2LTR TSD (independent evidence the call is right) are copied onto synthetic
  contigs with 3 kb of their real flanks. One obstacle is planted in one LTR,
  `delta` bp from an outer end: a 30 bp patch at 30% divergence, a deletion
  (10/50/200 bp), an insertion (50/300/1,000/5,000 bp of random sequence, or 300/1,000 bp
  cut from another family's element — a TE-like insertion), or nothing. The call is then
  cut at the obstacle exactly as LTRharvest / LTR_FINDER stop. Each truth element is
  left out of its own family's model (leave-one-out). Scored: recovered outer end
  exact / within 1 / within 5 bp, over-extended, wrong, left unchanged.
- **B2 — a real run.** The run's clean tables are copied and re-bounded. Measured:
  elements extended; among those with no TSD when called, the TSD rate at the new
  ends against a null with the right flank displaced by 1 kb; the same with family
  QC chosen on one half of each family's candidates and scored on the other (`cv`),
  so TSD-based QC cannot inflate its own validation; K2P shifts; time and memory.
- **B3 — false changes.** B1: untouched truth elements (no obstacle, no truncation)
  that get altered. B2: calls that already had a TSD and were altered, and how many
  of those keep a TSD.

## Decision rule (spec §9)

Maximise B1 exact-or-±1, subject to B3 false change ≤ 0.5% and B2 cross-validated
TSD gain ≥ 74.2% (the spike's figure). Ties go to the simpler, then the faster setting.

Neither numeric gate was reachable on real data at any setting tried in the sweep
below; a pre-registered fallback rule (same priorities, "not significantly worse
than the best cell" in place of the unreachable thresholds) was used instead. See
"What the benchmark could not measure" under Results.

## Reproduce

```bash
cd /data2/chris/poa_LTR/ltrquest_run          # any finished LTRquest run directory
export PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH
export PYTHONPATH=/data2/chris/poa_LTR/LTRquest/src
B=/data2/chris/poa_LTR/LTRquest/benchmarks/reboundary/bench.py
P="annua.nuclear_LTRs chaixii.nuclear_LTRs infirma.nuclear_LTRs supina.nuclear_LTRs"
G="annua.nuclear.fa chaixii.nuclear.fa infirma.nuclear.fa supina.nuclear.fa"
T=/data2/chris/poa_LTR/ltrquest_tools          # Kmer2LTR checkout at aa25f46
python $B b1 --run . --prefix $P --genome $G --tools-dir $T --out reboundary_bench/b1_NAME --n 1000
python $B b2 --run . --prefix $P --genome $G --tools-dir $T --out reboundary_bench/b2_NAME \
    --cache reboundary_bench/cache_NAME.pkl
python $B collect --out reboundary_bench
```

## Results

Data: 4 *Poa* genomes (annua, chaixii, infirma, supina), LTRquest run of 2026-09-19
(208,018 elements). Kmer2LTR aa25f46. Tables are pasted from `bench.py collect` and
from the four-stage tuning sweep (`tune.sh 1|2|3|4a|4b`, full ledger:
`.superpowers/sdd/2026-09-22-reboundary/progress.md`).

Neither of the plan's original gates (B2 cross-validated TSD gain >= 74.2%, B3 false
change <= 0.5%) was reachable at any setting tried — see "What the benchmark could
not measure" below. Every stage below therefore used a pre-registered fallback rule,
declared in the ledger *before* that stage's cells were run, that keeps the plan's
stated priorities and substitutes "statistically indistinguishable from the best
cell" for the unreachable numeric thresholds.

### Stage 1 — method x references

Six cells, B1 n=500 truths (4,486 contigs), B2 real, `min_identity=0.8
anchor_len=30 max_indel=5000 credit=200` (unchanged from the plan pending later
stages):

| cell | exact/±1 | over | wrong | B1 false | cv TSD | cv null | extended | TG..CA new |
|---|---|---|---|---|---|---|---|---|
| consensus_modal | 0.3788 | 0.0171 | 0.0163 | 0.016 | 0.5102 | 0.0092 | 18555 | 0.579 |
| consensus_tsd   | 0.3773 | 0.0176 | 0.0211 | 0.024 | 0.5263 | 0.0094 | 18051 | 0.5968 |
| subfamily_modal | 0.3976 | 0.0216 | 0.0176 | 0.024 | 0.5255 | 0.0092 | 19222 | 0.6013 |
| subfamily_tsd   | 0.3941 | 0.0186 | 0.0228 | 0.026 | 0.5274 | 0.0095 | 19402 | 0.5939 |
| nearest_modal   | 0.3966 | 0.0176 | 0.0223 | 0.018 | 0.5259 | 0.0096 | 19745 | 0.5873 |
| nearest_tsd     | 0.4029 | 0.0188 | 0.0248 | 0.022 | 0.5314 | 0.0093 | 20008 | 0.5901 |

Every cell cleared the cv-TSD >= 0.95 x max filter (>= 0.5048), so the fallback's
second step — lowest B1 false_change — was left to decide among three cells inside
one statistical sd of each other (consensus_modal 0.016, nearest_modal 0.018,
nearest_tsd 0.022; sd ~0.006 at n=500, p~0.02). Rather than call that noise a
winner, those three were re-run at n=3000 (declared *before* the re-run, so the
rule itself never moved):

| cell | exact/±1 | over | wrong | B1 false | controls |
|---|---|---|---|---|---|
| consensus_modal | 0.3247 | 0.0199 | 0.0205 | 0.0203 | 61/3000 |
| nearest_modal   | 0.3350 | 0.0190 | 0.0217 | 0.0193 | 58/3000 |
| nearest_tsd     | 0.3413 | 0.0192 | 0.0165 | 0.0187 | 56/3000 |

The B1-false ordering *reversed* between n=500 and n=3000 and all three pairs are
statistically indistinguishable there (two-proportion z, p = 0.65, 0.65, 0.87) —
direct evidence the n=500 ranking was noise. Falling through to exact/±1:
nearest_tsd and nearest_modal both beat consensus_modal (p = 0.0001, p = 0.017) but
tie with each other (p = 0.14); the plan's stated preference order (modal before
tsd) then picks nearest_modal.

**Winner: `method=nearest`, `references=modal`** — nearest posts the highest B1
recovery of any method and is confirmed significantly ahead of consensus at n=3000
(p = 0.017), while consensus — the plan's own original default — is now measurably
*worse* on recovery (p = 0.0001), a signal the noisier n=500 pass had hidden;
`modal` keeps TSD-based reference selection out of the TSD-based cv-TSD metric, so
that headline number stays independent of the knob it is validating.

### Stage 2 — min_identity x anchor_len x max_indel

18 cells, B1 n=500, method/references fixed at the stage-1 winner:

| min_identity | anchor_len | max_indel | exact/±1 | over | wrong | B1 false |
|---|---|---|---|---|---|---|
| 0.7 | 30 | 1000  | 0.3996 | 0.0238 | 0.0482 | 0.026 |
| 0.7 | 30 | 5000  | 0.4017 | 0.0236 | 0.0492 | 0.026 |
| 0.7 | 30 | 20000 | 0.4017 | 0.0236 | 0.0497 | 0.026 |
| 0.7 | 50 | 1000  | 0.4079 | 0.0236 | 0.0482 | 0.026 |
| 0.7 | 50 | 5000  | 0.4097 | 0.0238 | 0.0484 | 0.026 |
| 0.7 | 50 | 20000 | 0.4097 | 0.0238 | 0.0489 | 0.026 |
| 0.7 | 80 | 1000  | 0.4049 | 0.0241 | 0.0472 | 0.026 |
| 0.7 | 80 | 5000  | 0.4067 | 0.0243 | 0.0477 | 0.026 |
| 0.7 | 80 | 20000 | 0.4067 | 0.0243 | 0.0482 | 0.026 |
| 0.8 | 30 | 1000  | 0.3944 | 0.0176 | 0.0213 | 0.018 |
| 0.8 | 30 | 5000  | 0.3966 | 0.0176 | 0.0223 | 0.018 |
| 0.8 | 30 | 20000 | 0.3966 | 0.0176 | 0.0228 | 0.018 |
| 0.8 | 50 | 1000  | 0.4004 | 0.0173 | 0.0216 | 0.018 |
| 0.8 | 50 | 5000  | 0.4022 | 0.0178 | 0.0218 | 0.018 |
| 0.8 | 50 | 20000 | 0.4022 | 0.0178 | 0.0223 | 0.018 |
| 0.8 | 80 | 1000  | 0.4019 | 0.0173 | 0.0206 | 0.018 |
| 0.8 | 80 | 5000  | 0.4037 | 0.0178 | 0.0211 | 0.018 |
| 0.8 | 80 | 20000 | 0.4037 | 0.0178 | 0.0216 | 0.018 |

**Winner: `min_identity=0.8`, `anchor_len=30`, `max_indel=5000` (the plan's defaults
survive)** — every cell was statistically indistinguishable from the best on both
exact/±1 and B1 false at n=500, so the fallback fell through to the plan's stated
defaults; the table nonetheless shows a real, structural trade the pairwise tests
were individually underpowered to catch: relaxing `min_identity` to 0.7 buys about
+1 point of recovery but more than doubles the wrong rate (0.047-0.050 vs
0.021-0.023 for every matched 0.8 cell — a paired sign test over the 9 pairs is
9/9, p ~ 0.004) — not worth it. `anchor_len` and `max_indel` barely move recovery
among the 0.8 cells (spread 0.3944-0.4037), so 30/5000 is a free, and the more
conservative, choice.

### Stage 3 — credit

B2 real (cached family phase) + B1 n=500, method/references/min_identity/
anchor_len/max_indel fixed at the stage-1/2 winners. **cv TSD is identical
(0.5259) at every credit value in this stage** — see "What the benchmark could not
measure" for why that is expected, not a null result.

Initial grid (the plan's tested values):

| cell | cv TSD | tsd_gain | b3_real | b3_kept | reverted | extended | TG..CA new | exact/±1 | over+wrong | B1 false |
|---|---|---|---|---|---|---|---|---|---|---|
| credit_0     | 0.5259 | 0.4629 | 0.0126 | 0.2462 | 8698 | 15905 | 0.5281 | 0.2559 | 0.0402 | 0.0200 |
| credit_50    | 0.5259 | 0.4976 | 0.0141 | 0.2776 | 5931 | 18586 | 0.5613 | 0.3023 | 0.0394 | 0.0220 |
| credit_model | 0.5259 | 0.4997 | 0.0142 | 0.2797 | 5805 | 18702 | 0.5632 | 0.3164 | 0.0394 | 0.0220 |
| credit_200   | 0.5259 | 0.5224 | 0.0147 | 0.3014 | 4251 | 19745 | 0.5873 | 0.3966 | 0.0399 | 0.0180 |
| credit_1000  | 0.5259 | 0.5413 | 0.0155 | 0.3223 | 1893 | 21043 | 0.6114 | 0.4714 | 0.0349 | 0.0180 |

`credit_1000` won at the edge of the tested grid with `tsd_gain` still rising, so
rather than accept an edge value the grid was extended to 5000 and 20000 before
deciding:

| cell | tsd_gain | b3_real | reverted | extended | TG..CA new | exact/±1 | over+wrong |
|---|---|---|---|---|---|---|---|
| credit_1000  | 0.5413 | 0.0155 | 1893 | 21043 | 0.6114 | 0.4714 | 0.0349 |
| credit_5000  | 0.5534 | 0.0154 |  105 | 21505 | 0.6242 | 0.4834 | 0.0326 |
| credit_20000 | 0.5537 | 0.0154 |   22 | 21497 | 0.6246 | (tied) | (tied) |

**Winner: `credit=5000`** — the curve saturates between 1000 and 5000
(`tsd_gain` 0.5413 -> 0.5534, `extended` 21043 -> 21505, `reverted` 1893 -> 105);
20000 is indistinguishable from 5000 (the plateau), and 5000 strictly dominates
1000 on every axis at once (including `b3_real`, 0.0154 vs 0.0155 — better, not
worse), so the smaller of the two plateau values is preferred on the principle
that it leaves Kmer2LTR more veto authority for the same outcome. At credit=5000
Kmer2LTR reverts only ~0.5% of proposals (105 of ~21,600); B1 confirms withdrawing
that veto is correct, not risky — `over+wrong` *falls* as credit rises (0.0402 at
credit 0 -> 0.0326 at credit 5000), meaning the extensions Kmer2LTR was vetoing at
low credit were disproportionately the *correct* ones.

### Stage 4a — max_ratio x min_copies

B2 real (cached family phase), method/references/min_identity/anchor_len/
max_indel/credit fixed at the stage-1/2/3 winners:

| cell | cv TSD | extended | TG..CA new |
|---|---|---|---|
| max_ratio=1.10, min_copies=5 or 10   | 0.5287 | 21204 | 0.6306 |
| max_ratio=1.10, min_copies=20        | 0.5279 | 21043 | 0.6318 |
| max_ratio=1.15, min_copies=5 or 10   | 0.5259 | 21505 | 0.6242 |
| max_ratio=1.25, min_copies=5 or 10   | 0.5259 | 21505 | 0.6242 |
| max_ratio=1.15 or 1.25, min_copies=20| 0.5251 | 21345 | 0.6253 |

**Winner: `max_ratio=1.15`, `min_copies=10` (the plan's defaults survive)** — cv
TSD is statistically flat across all five distinct cells, so the fallback keeps
the plan's stated defaults; the table also shows two knobs are inert over part of
their range (see "What the benchmark could not measure"), so tuning them further
would buy nothing.

### Stage 4b — qc_min_n x untested

B2 real (cached family phase), every other knob fixed at its stage-1/2/3/4a
winner. *(This stage's table and ruling were completed directly from
`tune4b/*/summary.json` using the stage-4 decision rule the ledger declared
before stage 3 ran — the ledger recorded the stage-4b launch but the run finished
after the last write-up entry.)*

| cell | cv TSD | extended | TG..CA new | B2 already-TSD changed |
|---|---|---|---|---|
| qc_min_n=3, untested=accept  | 0.5256 | 21464 | 0.6249 | 0.0153 |
| qc_min_n=3, untested=skip    | 0.5247 | 21330 | 0.6259 | 0.0151 |
| qc_min_n=5, untested=accept  | 0.5259 | 21505 | 0.6242 | 0.0154 |
| qc_min_n=5, untested=skip    | 0.5246 | 21256 | 0.6260 | 0.0151 |
| qc_min_n=10, untested=accept | 0.5232 | 21531 | 0.6241 | 0.0154 |
| qc_min_n=10, untested=skip   | 0.5257 | 21018 | 0.6265 | 0.0150 |

**Winner: `qc_min_n=5`, `untested=accept` (the plan's defaults survive)** — cv TSD
is flat across all six cells (0.5232-0.5259, well inside noise at cv_n ~13,000-
13,600), so the fallback's second step decides: four cells fall within 1% of the
most-extended cell (`qc_min_n=10, accept` at 21531) — `qc_min_n=10/accept`,
`qc_min_n=5/accept` (21505), `qc_min_n=3/accept` (21464), `qc_min_n=3/skip`
(21330) — and the plan's stated default (`qc_min_n=5, untested=accept`) is one of
them, so it is kept.

### Chosen defaults

| knob | value | stage | deciding numbers |
|---|---|---|---|
| method | nearest | 1 | n=3000: exact/±1 0.3350 vs consensus 0.3247 (p=0.017); consensus is now significantly worse |
| references | modal | 1 | ties with `tsd` on every axis (p=0.14 on exact/±1); plan's stated preference order (modal before tsd) |
| min_identity | 0.8 | 2 | 0.7 doubles the wrong rate (0.047-0.050 vs 0.021-0.023), paired sign test 9/9, p~0.004 |
| anchor_len | 30 | 2 | spread 0.3944-0.4037 exact/±1 across all 0.8 cells — inside noise; plan default kept |
| max_indel | 5000 | 2 | spread 0.3944-0.4037 exact/±1 across all 0.8 cells — inside noise; plan default kept |
| credit | 5000 | 3 | tsd_gain 0.5413->0.5534, extended 21043->21505, reverted 1893->105 vs credit=1000; credit=20000 indistinguishable (plateau) |
| max_ratio | 1.15 | 4a | identical to 1.25 (21505 extended both) — the length-ratio gate only binds below 1.15 |
| min_copies | 10 | 4a | identical to 5 (families with 5-9 copies already fail the MIN_REFS reference floor) |
| qc_min_n | 5 | 4b | within 1% of the most-extended cell (21505 vs 21531 max); plan default kept |
| untested | accept | 4b | same tie-break as qc_min_n; within-1% survivor and the plan default |
| subfamily_jaccard | 0.5 | 1 | only used by `method=subfamily`, which lost stage 1; not independently tuned, plan default kept |

### At the chosen defaults

Real run, 4 Poa genomes, 208,018 elements
(`reboundary_bench/tune4b/b2_n5_accept/summary.json`):

candidates 65,712 -> extended 21,505; TSD at the new boundary 0.5534 vs the
displaced-flank null 0.0100; cross-validated 0.5259 vs null 0.0096 (n=13,507);
TG..CA 0.3191 -> 0.6242; extension median 316 bp, q90 1,273 bp; dK2P median
0.0018; 79 s, 4.07 GB.

| element age (K2P) | candidates | extended | TSD recovery |
|---|---|---|---|
| < 0.005 | 6,937 | 1,935 | 0.6475 |
| 0.005-0.02 | 19,215 | 6,950 | 0.6627 |
| 0.02-0.05 | 23,029 | 8,215 | 0.5585 |
| >= 0.05 | 16,531 | 4,405 | 0.3341 |

Rejections (of 65,712 candidates):

| reason | n |
|---|---|
| gate_identity | 31,618 |
| family_qc_failed | 4,634 |
| overlaps_element | 3,923 |
| kmer2ltr_not_pass | 3,475 |
| host_exceeded | 269 |
| engulfs_element | 154 |
| kmer2ltr_reverted | 105 |
| length_filter | 29 |

For comparison, the plan's original defaults (`method=consensus, credit=200`, all
else the same) gave 18,555 extended with cv TSD 0.5102 (Stage 1's
`consensus_modal` row) — the tuned defaults extend 15.9% more elements at
slightly better specificity.

### What the benchmark could not measure

These four points matter more than which cell won, because they bound how much
confidence the numbers above deserve:

- **The plan's gates were unreachable on real data.** The plan called for B2
  cross-validated TSD gain >= 74.2% and B3 false change <= 0.5%; the best measured
  across the whole sweep was cv TSD 0.531 and B3 false 0.016 — neither gate was
  ever in reach. The 0.742 figure came from the spike's QC-pass-only subset, a
  much smaller denominator than the sweep's cross-validated figure, so it was
  never a like-for-like target. A pre-registered fallback rule (highest
  recovery / lowest false-change / plan-default tie-break, each unreachable
  numeric gate replaced by "not significantly worse than the best cell") was
  declared and used instead of the plan's literal rule, and Stage 1 was re-run at
  n=3000 to resolve a genuine tie *by collecting more data*, rather than by
  amending the rule after seeing the n=500 numbers.
- **`cv TSD` cannot see past Kmer2LTR.** It is computed from `proposals.tsv` — the
  family phase's output, *before* Kmer2LTR arbitration — so it is identical
  across every value of `credit` (and any other post-arbitration knob) by
  construction. Stage 3's flat 0.5259 column is not evidence that `credit` does
  nothing; it is evidence that this particular metric cannot discriminate a
  post-arbitration knob at all. Read `tsd_gain`, `reverted`, and B1 `exact/±1`
  for that stage's real signal instead.
- **B1's `false_change` is structurally conservative.** B1's synthetic contigs
  carry no neighbouring annotations, so the conflict rules that block a
  kilobase-scale anchor jump on real data (`overlaps_element`, mutual-conflict
  resolution) cannot fire in B1 at all — a jump B2's SpanIndex would catch is
  scored as a false change in B1 regardless. Separately, a B1 control whose
  *original* finder call was already short of the true end is counted as a false
  change the moment re-boundarying extends it, even when the extension is
  biologically correct. Both effects push B1 `false_change` upward relative to
  what a real run would show.
- **Two knobs are inert over part of their tested range.** `min_copies=5` and
  `min_copies=10` give identical results at every `max_ratio` (Stage 4a):
  families with 5-9 copies already fail the `MIN_REFS` reference floor inside
  `select_references`, so lowering `min_copies` below 10 cannot admit anything.
  `max_ratio=1.15` and `max_ratio=1.25` are likewise identical: the length-ratio
  QC gate only binds below 1.15. Neither knob needs tuning again inside these
  ranges.
