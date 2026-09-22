# Family-guided LTR-RT re-boundarying (`--reboundary`) — design

Date: 2026-09-22 · Branch: `reboundary` · Status: approved design, pre-implementation

## 1. Problem

LTRharvest / LTR_FINDER extend an LTR pair outward from a seed and stop at the first
obstacle between the element's two LTRs. Kmer2LTR re-calls every pair but can only trim
inward. So when one LTR carries an indel or a mutation-dense patch near its end, the call
stops short. In a family MSA those copies "start late and end early".

Evidence from a feasibility spike on the 4-genome *Poa* run
(`/data2/chris/poa_LTR/ltrquest_run/reboundary_spike/memo.md`, 121 families, 34,625 copies):

| finding | number |
|---|---|
| copies whose true ends a family LTR model recovers (QC-pass families) | 10.9% (2,932 / 26,925) |
| TSD at called ends → at family-guided ends (those copies) | 6.2% → 72.8% (null 1.2%) |
| same, copies with no TSD when called | 74.2% gain one (null 1.2%) |
| TG..CA at called → new ends | 35.8% → 91.3% |
| median / q90 sequence recovered | 304 / 1,171 bp |
| obstacle within ±15 bp of the old end (indel or ≥20% mismatch) vs control site | 80% vs 24% |
| large gaps (≥50 bp) at the old end: deletion in one LTR / insertion | 72% / 28% |
| truncations fixed by falling back to the raw LTRharvest/LTR_FINDER call | 59 of ~3,800 ends |
| Kmer2LTR given the widened record (no credit) reverts it to the old end | 25%; family ends carry a TSD in 78% of those |
| Kmer2LTR `tsd_credit` 0 / 50 / 200 / 1000 bits → family ends kept | 53 / 71 / 77 / 86%; controls unchanged 100% |
| trimming (model shorter than the call) | wrong: TSD 45% → 12% |

## 2. Goals / non-goals

Goals
- An opt-in stage that extends truncated LTR-RT calls to the ends their family supports.
- Kmer2LTR stays the single source of truth for every Kmer2LTR column (inner LTR ends,
  divergence, K2P, time, TSD, motif, CIGAR). The family supplies evidence, not numbers.
- Runs inside the pipeline (`ltrquest --reboundary`) and post hoc on a finished run
  (`ltrquest-reboundary --posthoc`), which updates the run in place with a backup.
- Defaults chosen by benchmark (§9), not by taste.

Non-goals (v1)
- Trimming over-extended calls. Rule: extend only.
- Re-clustering families. An element is re-bounded with its own family's model, so its
  membership does not change. The pooled `*_all_ltr.*` clustering files stay as the
  family basis and are documented as pre-re-boundary.
- Re-computing nesting depth. Extensions that would change nesting are rejected with an
  explicit reason (§5.5). The spike found 0 of 692 large-gap cases needing it.
- Touching the raw `depth{N}_ltr.*` tables. They are detection outputs (`ltrquest.record`),
  so only the `_clean_` tables, `_clean_` FASTAs, GFF3s and plots change.

## 3. Terms

- **Forward frame.** Genome coordinates, 1-based inclusive. The table's `ltr5_*` / `ltr3_*`
  are the genomic-left / right LTR even when `orientation` is `-` (recover_strand.py:196).
  The module works in the forward frame throughout; biological 5′/3′ is used only when
  building models.
- **Outer ends.** The element's first and last base, `L0` and `R1`.
- **Inner ends.** `L1` (end of the left LTR) and `R0` (start of the right LTR).
- **Family model.** One or more LTR sequences with trusted termini, derived from a
  family's members (§5.1).
- **Proposal.** New outer ends for one element, plus the evidence behind them.
- **Credit.** Bits of external evidence, handed to Kmer2LTR as `tsd_credit`, that the
  proposed record termini are the element's boundaries (Kmer2LTR `align._extend`).
- **Obstacle.** The indel / mismatch patch between the element's two LTRs at the old outer
  end. Reported per element, never used to find ends.

## 4. Architecture

```
pipeline, --reboundary OFF (unchanged):
  per genome: [recover-strand align] -> annotate -> gff3 -> [recover-strand apply]

pipeline, --reboundary ON (final FP attempt only, inside the staging dir):
  per genome: [recover-strand align] -> annotate -> [recover-strand apply]
  pooled:     ltrquest-reboundary --prefix P1..Pn --genome G1..Gn      (all genomes at once:
                                                                        families are pooled)
  per genome: gff3 --reboundary-map <P>_reboundary.tsv
  (plots.sh after promotion, as today)

post hoc on a finished run directory:
  ltrquest-reboundary --posthoc --indir RUN --prefix P1..Pn --genome G1..Gn
    backup -> re-bound -> gff3 (with map) -> plots.sh
```

Moving `recover-strand apply` ahead of gff3 is safe. apply rewrites only the FASTAs and the
`orientation` column, and gff3 reads neither (workflows/ltrquest.nf says the same). It
must run before re-boundarying, because apply keys on the old element names.

### Files

| file | change |
|---|---|
| `src/ltrquest/ltr_model.py` (new) | model building (3 methods, 2 reference strategies), end polishing, family QC, placement (core + terminal anchor), proposals. Pure functions over sequences plus a genome accessor, unit-testable |
| `src/ltrquest/reboundary.py` (new) | CLI `ltrquest-reboundary`: load, pool, propose, check conflicts, Kmer2LTR arbitration, rewrite, sidecar/map, post-hoc driver |
| `src/ltrquest/kmer2ltr.py` | `api(tools_dir)`: import Kmer2LTR's Python API (`align.classify`, `genome.Window/orient/annotate/find_tsd`, `scoring.bits`), with a version/signature guard |
| `src/ltrquest/gff3.py` | `--reboundary-map`: translate lookups (families, strand provenance) for renamed elements; attributes `boundary_source`, `boundary_shift` |
| `src/ltrquest/scripts/ltrquest.sh` | `--reboundary` flag; annotation stage restructured when on; file hygiene (carry-forward, promotion) |
| `pyproject.toml` | console script `ltrquest-reboundary` |
| Nextflow | `modules/local/ltrquest/reboundary/main.nf`, `params.reboundary`, workflow wiring |
| `benchmarks/reboundary/` (new) | harness + README with results (§9) |
| docs | `docs/outputs.md` (sidecar, GFF3 attributes, clean-vs-raw note), README section, CHANGELOG |
| tests | synthetic fixtures + unit/integration tests (§10) |

No new dependencies: mafft, parasail-python, pywfa, pyfaidx and numpy are already in
`environment.yml`, `recipe/meta.yaml` and the Dockerfile.

## 5. Algorithms

### 5.1 Family models (`ltr_model.py`)

Eligibility: families with at least `min_copies` members pooled over all genomes (knob,
default 10). Sequence comes from the **original** genomes, the same rule strand recovery
follows. Genome access is random-access pyfaidx: never a whole genome in memory, never a
full streaming pass per family.

Reference copies, knob `references`:
- `modal`: stranded copies with both called LTR lengths within ±15% of the family's modal
  called LTR length. Modal length = densest 5%-wide log-length bin, smoothed by its
  neighbours; the median of that bin.
- `tsd`: stranded copies whose called ends carry a Kmer2LTR TSD (`tsd` column not `.`/`NA`).
  Falls back to `modal` below 10 copies.
- Take up to 40 youngest (lowest `k2p`) plus 40 random (seeded by `crc32(family)`).
- Exclusion list honoured (benchmark leave-one-out).

Model method, knob `method`:
- `consensus` (spike method). Each reference contributes its 5′ LTR and 3′ LTR padded by
  200 bp, oriented 5′→3′, all aligned together with `mafft --auto`. A column is LTR when
  the 5′-LTR rows and 3′-LTR rows agree on the majority base, each ≥40%. Beyond the outer
  end the 5′ rows are flank; beyond the inner end the 3′ rows are flank; either way the
  rows disagree. The span runs from the first to the last sustained agreement block
  (≥15 of 21 columns), so internal dips are ignored. Consensus = majority base over
  columns ≥50% occupied.
- `subfamily`. References (up to 300) are clustered by LTR k-mer Jaccard, greedily around
  the youngest-first centroids (threshold knob). Each cluster with ≥10 members gets a
  `consensus` model; the leftovers pool into a family model.
- `nearest`. Templates are the references' own LTRs, used only when the template's called
  outer ends agree (±2 bp) with where the family's `consensus` model places them, so a
  template's termini are confirmed by its family rather than taken on trust. Each target
  uses its top-3 templates by k-mer similarity; each outer end is the median of their
  placements.

**End polishing** (all methods). This applies the definition of an LTR end directly: the
last position where an element's two LTRs still agree. The model is placed in both LTRs of
every reference, and the base pairs just beyond each model end are compared. Extend while
the pair-agreement fraction is ≥0.7; trim while it is below. At most 8 bp per round,
3 rounds. No TSD or motif is used anywhere in building a model.

**Family QC** (per family, never per element):
1. `model_len / modal_called_len ≤ max_ratio` (knob, default 1.15). In the spike every
   failing family (8 of 8) exceeded it and gained 0% TSDs.
2. TSD enrichment. Among the family's proposals whose called ends lack a TSD, the TSD rate
   at the proposed ends must beat the displaced-flank null. That means a one-sided
   binomial test, p < 0.01, run only with ≥ `qc_min_n` proposals (default 5). TSDs are
   found with Kmer2LTR's own `genome.find_tsd`, using its k and shift settings, on windows
   cut by pyfaidx, so the semantics match the `tsd` column exactly.
   Families below `qc_min_n` are **untested**, an explicit branch with its own policy knob
   (`accept-untested` / `skip-untested`). They are never silently counted as pass or fail
   (`lessons.md`: every classifier needs an untested branch).
   TSDs are a permitted credibility signal (`lessons.md`). TG..CA is not used, because the
   finders search for it.

### 5.2 Placement

For each element, for each model (the best by k-mer prefilter; `subfamily` / `nearest`
try their top candidates):

- **Core.** Glocal alignment (parasail `sg_dx`: model end to end, genome ends free; match
  2 / mismatch −3 / open 5 / extend 2) in a window around each called LTR:
  `called LTR ± (model_len − called_len + 100)`. Orientation comes from the `strand`
  column; for `.`, the better of both.
- **Terminal anchor.** New in v1, for large indels near the end. The model's outer
  `anchor_len` bp (default 80) are aligned locally (parasail `sw`) in a window reaching
  `max_indel` bp beyond the called outer end and 200 bp inside it.
  - `max_indel` is a knob; candidate defaults are 1 kb / 5 kb / 20 kb.
  - The anchor is accepted when its bit score ≥ `anchor_bits` (knob).
  - It must be collinear with the core placement: on the outer side of it, the same
    orientation, and the gap to the core ≤ `max_indel`.
  - No annotated element may lie between them (§5.5).
  - When accepted, the anchor sets that outer end; otherwise the core does.
- **Gates** at each end that would move outward:
  - outer-`TERM` identity ≥ `min_identity` (knob 0.7 / 0.8; TERM = 30);
  - whole-LTR identity ≥ 0.6.
- **Never trim.** An outer end proposed inside the called one is set back to the called
  end. An element is a candidate only if at least one outer end moves outward by
  ≥ `min_ext` (default 5 bp).
- **Obstacle report.** The element's two LTRs are aligned pairwise (open 10 / extend 1)
  and the old end is mapped onto that alignment. Recorded: longest gap within ±15
  columns, mismatch rate in the 30 bp just outside and just inside. Reported only.

### 5.3 Kmer2LTR arbitration

Candidate record = the existing clean record, with the recovered segments spliced on and
turned to the genome-forward frame (nested inners stay masked, exactly as when the element
was detected), named `chrom:L0'-R1'#<same class suffix>`.

```
ctx    = genome.orient(seq, Window(up, head, tail, down))   # Window cut by pyfaidx, Kmer2LTR's PAD/PROBE
result = align.classify(name, seq, period_rule="outermost",
                        mutation_rate=mu, tsd_credit=credit)
result = genome.annotate(result, seq, ctx, genome.Options())
```

- Credit, knob `credit`: a constant (0 / 50 / 200 / 1000 bits), or `model`: the smaller of
  the extended ends' anchor/core terminal bit scores, via Kmer2LTR's `scoring.bits`.
- Apply detection's Step 8a filter to the result: min LTR 100, aln 90, length ratio 0.65,
  element 300.
- Accept iff status is `pass`, Kmer2LTR's settled termini are outward of the called ends
  by ≥ `min_ext` on at least one side, and not inside them on either side. Kmer2LTR may
  settle anywhere between the old and proposed ends; it has the final say.
- Rows are rebased exactly as `detect.rebase_to_trimmed`, so `assert_bounded` holds.
- `mu`: the run's `--mutation-rate` (pipeline), else `mutation_rate` from
  `<prefix>.detect.json` settings, else 3e-8. A warning is printed when defaulted.

### 5.4 Rewriting (clean side only)

- **Table row.** Columns 1–29 come from the Kmer2LTR result.
  - `orientation` is restated to the stored orientation of the existing record: a storage
    fact, kept.
  - `strand`, `family` and `domains` are preserved. Domains are absolute genomic
    coordinates and don't move.
  - `nest_status` is re-keyed.
- **nest_status.** Every other row's references to the old key are rewritten to the new
  key, in all depth tables of that genome.
- **FASTA record.** Recovered segments are spliced onto the existing record in its stored
  orientation: `+` = left-segment + record + right-segment; `-` = rc(right) + record +
  rc(left). The header becomes the new `seq_id`. The segments contain no annotated element
  (§5.5), so they need no masking.
- **Host records.** Every element whose span contains the extended element (its
  `nest-inner` hosts) gets the recovered segments painted with the extended element's
  depth letter (`reconcile.IUPAC_DEPTH_SEQ`), mirrored when the host is stored `-`. This
  matches `reconcile.apply_depth_masking`.
- **Writes** are atomic (tmp + rename), table and FASTA per depth, FASTA first, then the
  table. That is the `recover_strand.apply` order.
- **Idempotent.** A second run finds nothing to extend.

### 5.5 Conflicts (per genome, over all elements of all depths)

Let the old span be [s, e] and the new span [s′, e′] ⊇ [s, e]. Reject the extension (reason
recorded) if any other element X:
- **host_exceeded:** X contains [s, e] but not [s′, e′];
- **engulfs_element:** X lies within an added segment;
- **overlaps_element:** X partially overlaps an added segment.

Elements are processed in a fixed order (genome, chrom, start). An accepted extension
updates the interval index before the next element is checked, so two neighbours cannot
both grow into the same gap.

### 5.6 Outputs

`<prefix>_reboundary.tsv`, one row per **candidate** element (accepted or rejected):

```
#old_seq_id  new_seq_id  family  method  model_id  decision  reason
 ext5  ext3  (bp added at biological 5'/3')   end_source5  end_source3 (core|anchor)
 id_outer5  id_outer3  credit_bits  k2l_status
 tsd_called  tsd_new  tsd_null   k2p_called  k2p_new
 obstacle5  obstacle3   (e.g. gap:96 | mm:0.23 | none)
```

- `decision`: `extended` | `rejected`. `reason` is one of: `family_qc_failed`,
  `family_untested`, `gate_identity`, `host_exceeded`, `engulfs_element`,
  `overlaps_element`, `kmer2ltr_reverted`, `kmer2ltr_not_pass`, `length_filter`.
- `new_seq_id` = `.` when rejected. This file doubles as the key map for gff3.
- A per-family summary goes to the log: model length, modal length, QC outcome,
  candidate and accepted counts.
- **GFF3.** Re-bounded elements get `boundary_source=family_model` and
  `boundary_shift=<ext5>,<ext3>`. Everything else is unchanged.

## 6. Pipeline integration (`ltrquest.sh`)

- `--reboundary` (flag, default off).
- Help text sits beside `--strand-recovery`. Not part of `detection_settings`: it is
  post-detection, so genome reuse is unaffected.
- **Annotation stage when on:**
  - Loop 1, per genome: strand-recovery align, annotate, strand-recovery apply.
  - One pooled `ltrquest-reboundary` call: all prefixes, their original genomes, threads,
    mutation rate, tools dir.
  - Loop 2, per genome: gff3 with `--reboundary-map`.
  - When off, the existing loop is byte-for-byte unchanged.
- **Failure.** The re-boundary step warns and continues. Its writes are atomic per file,
  and a failure before any write leaves the clean tables untouched; a failure mid-rewrite
  restores from a staging copy it takes first.
- **Hygiene.**
  - `carry_forward_genome` skips `*_reboundary.*` and `*_pre_reboundary`.
  - `promote_fp_outputs` removes a stale `${dst}/${p}_reboundary.tsv` and
    `${dst}/${p}_pre_reboundary/` before promoting.
- The Kmer2LTR tools dir is resolved exactly as the merged stage does
  (`resolve_merged_tools_dir` / `ensure_kmer2ltr_dir`).

## 7. Post-hoc mode

```
ltrquest-reboundary --posthoc --indir RUN --prefix P1 ... Pn --genome G1 ... Gn \
    [--method ... --threads N] [--no-plots]
```

- **Preconditions (fail fast):**
  - annotated `_clean_` tables (`family`, `strand` columns) exist for every prefix;
  - all prefixes share one family namespace;
  - genomes are given one per prefix.
- **Backup.** `<prefix>_pre_reboundary/` receives the clean tables and FASTAs, both GFF3s
  and the plots directory, moved not copied, once. If the backup already exists, **inputs
  are read from it**, so re-runs with other settings always start from the originals.
- **Then:** re-bound in place → gff3 per prefix (`--consensus-cluster` / `--family-prefix`
  found as in the annotation stage, `--recovered-strands` if the sidecar exists,
  `--reboundary-map`) → `plots.sh` per prefix unless `--no-plots`.
- **Restoring:** `ltrquest-reboundary --restore --indir RUN --prefix ...` moves the backup
  back.

## 8. Nextflow

`params.reboundary` (default false). Process `LTRQUEST_REBOUNDARY` takes, collected over
all samples:
- every sample's annotated clean tables and FASTAs (after `RECOVERSTRAND_APPLY` when strand
  recovery is on);
- the genomes;
- the consensus cluster table.

It emits per-sample tables, FASTAs and sidecars, which feed `LTRQUEST_GFF3` (plus the
map) and `LTRQUEST_PLOTS`. With the flag off, channels are exactly as today. A stub block
is provided.

## 9. Benchmarks (`benchmarks/reboundary/`)

Harness `bench.py` runs against any finished LTRquest run directory. It never writes into
the run: it works in a bench directory with copies or symlinks. Results go to
`benchmarks/reboundary/README.md` (tables) plus a TSV per experiment.

**B1 — planted obstacles (ground truth).**
- **Truth set:** ~2,000 elements with a Kmer2LTR TSD at their called ends, stratified by
  clade and age.
- **Build:** each gets a synthetic contig of 3 kb flank + element + 3 kb flank, taken from
  its genome. One obstacle is planted in one LTR, at distance δ ∈ {20, 50, 100, 200, 400}
  bp from its outer end (δ < LTR/2):
  - substitution patch (30 bp at 30%);
  - deletion (10 / 50 / 200 bp);
  - insertion (50 / 300 / 1,000 / 5,000 bp; random or a real segment from another family);
  - none (bare truncation).
- **The call is cut at the obstacle,** as a finder would: the outer end and the homologous
  inner end of the other LTR.
- **Leave-one-out:** the truth copies of a family are excluded from its references.
- **Score:**
  - recovered outer end exact / ±1 / ±5 bp;
  - over-extension (beyond truth);
  - wrong (> 5 bp off);
  - left unchanged (fail-safe);
  - broken down by obstacle type, size and δ.

**B2 — the Poa run.** Post-hoc on a copy of the run. Measured:
- candidates and accepted;
- TSD at called / new / null ends, among copies without a TSD when called;
- TG..CA at called / new (reported, not a criterion);
- K2P shift;
- reject reasons;
- per clade and age.

Family QC is chosen on a random half of each family's copies and evaluated on the other
half, so TSD-based QC cannot inflate the reported gain.

**B3 — false changes.**
- B1 truth elements with no obstacle and no truncation: fraction altered.
- B2 copies with a TSD at their called ends: fraction altered, and TSD kept.

**Tuning order.** Proposals are cached (`--dump-proposals`), so later stages re-run cheaply.
1. `method` × `references` on B1 + B2 (6 combinations).
2. Placement: `min_identity`, `anchor_len`, `anchor_bits`, `max_indel` on B1.
3. `credit` on cached proposals.
4. Family QC (`max_ratio`, `qc_min_n`, untested policy) on B2 cross-validated.

**Decision rule.** Maximise B1 exact-or-±1 recovery, subject to:
- B3 false-change ≤ 0.5%;
- B2 cross-validated TSD gain not lower than the spike's (74% vs a 1.2% null).

Ties go to the simpler or faster setting. Chosen defaults are written into the CLI and
`benchmarks/reboundary/README.md`.

**Resources:** wall time and peak RSS (`/usr/bin/time -v`) for B2. Target on Poa (212k
elements, 4 genomes): ≤1 h at 32 threads, ≤16 GB.

## 10. Testing

All tests are pytest with synthetic fixtures and no real data. `tests/reboundary_fixtures.py`
builds a ~80 kb genome holding one family of 14 copies of a synthetic LTR-RT:
- LTR 400 bp, internal 2 kb, TG..CA, 5 bp TSDs, ~2% divergence.
- Three truncated calls, one each for a 60 bp deletion, a 30 bp 30% patch, and a 1.5 kb
  insertion.
- One nested pair; one neighbour placed to trigger each conflict rule; one old copy with a
  decayed TSD.

The fixture writes 33-column clean tables and matching FASTAs (with IUPAC-masked inners),
exactly as LTRquest does.

Unit tests:
- `consensus` recovers the planted LTR termini exactly (and each method);
- end polishing;
- core placement;
- anchor placement across the 1.5 kb insertion;
- never-trim;
- family QC incl. the untested branch;
- each conflict rule;
- FASTA splice in both orientations;
- host re-masking (mirrored for `-`);
- nest_status re-keying;
- sidecar/map;
- gff3 `--reboundary-map` lookups;
- post-hoc backup/re-run/restore;
- idempotence.

Kmer2LTR-dependent tests:
- skip when Kmer2LTR's API is not importable (`kmer2ltr.api`), matching how
  `conftest.py` skips on missing data;
- one integration test runs the whole CLI on the fixture.

Regression: with `--reboundary` off, the annotation-stage outputs are identical
(driver-level test using the existing end-to-end marker `slow`).

Python ≥3.9 syntax (`from __future__ import annotations`), ruff line length 100.

## 11. Implementation phases

1. `kmer2ltr.api` + `ltr_model.py` (consensus, placement core, QC) + fixtures and tests.
2. `reboundary.py`: proposals, conflicts, arbitration, rewrite, sidecar; CLI post-hoc mode
   without plots.
3. Benchmark harness (B1, B2, B3). First results with `consensus`.
4. Terminal anchor, `subfamily`, `nearest`, `references=tsd`. Benchmark each and tune (§9).
5. gff3 map and attributes; post-hoc regeneration (gff3 + plots) and restore.
6. `ltrquest.sh` integration + hygiene; Nextflow process.
7. Docs (outputs.md, README, CHANGELOG); apply to the Poa run post hoc (with a memo).

## 12. Open knobs (resolved by §9, recorded in the benchmark README)

| knob | candidates |
|---|---|
| `method` | consensus, subfamily, nearest |
| `references` | modal, tsd |
| `subfamily_jaccard` (LTR 11-mer Jaccard to join a cluster) | 0.3, 0.5, 0.7 |
| `min_copies` | 5, 10, 20 |
| `min_identity` | 0.7, 0.8 |
| `anchor_len` / `anchor_bits` | 60 / 80 / 100 bp; bit threshold from B1 ROC |
| `max_indel` | 1 kb, 5 kb, 20 kb |
| `credit` | 0, 50, 200, 1000, model |
| `max_ratio` | 1.10, 1.15, 1.25 |
| `qc_min_n` / untested policy | 3, 5, 10 / accept, skip |
