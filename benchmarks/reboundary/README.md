# Re-boundarying benchmarks

How well `--reboundary` finds the true ends of truncated LTR-RT calls, and how to
measure it again. How the stage works:
[docs/outputs.md §8](../../docs/outputs.md#8-re-boundarying---reboundary).

## Two independent checks

The new ends come only from templates, so two signals the stage never uses can check
them:

- **Target-site duplication (TSD).** An insertion duplicates a few host bases
  (usually 5) on both sides, so a correct end lands on a TSD and a wrong one rarely
  does. Templates carry TSDs at their own loci, which are other insertions, and
  nothing moves an end toward a TSD.
- **A curated library.** Where a curated LTR library can be placed on an element, its
  terminus shows where the end should be. This check is not part of `bench.py`.

## Results on rice (*O. sativa*)

One genome, 6,709 elements, 32 threads. This design and the pre-release one
(per-family LTR models, outward moves only) ran on the same tables.

| | pre-release | this design |
|---|---|---|
| elements changed | 560 | 1,726: 1,705 moved, 21 split calls merged |
| kinds of move | outward only | 1,203 outward, 300 trimmed, 202 one of each |
| curated terminus at the new end | 98.6% of 576 ends | 99.5% of 1,965 ends |
| exact 5-bp TSD at the moved end (10–100 bp away) | 46.0% (0.17%) | 50.2% (0.14%) |
| TSD at the new ends (right end moved 1 kb) | 55.5% (0.5%) | 50.1% (1.1%) |
| K2P change where added bases have no partner (median) | +0.090 | 0.000 |
| elements gaining more than 0.05 K2P | 103 | 18 |
| calls with a TSD that were moved | 22 | 0 |
| run time | 203 s | 35 s |

- Where the two designs put the same end in different places and the curated library
  can judge (156 ends), it backs this design's end 148 times and the pre-release
  one's 6 times.
- The pre-release design aligned bases without a partner as substitutions, which
  inflated K2P and age. The 18 elements that still gain more than 0.05 pair genuine
  but more diverged homology beyond a short called LTR pair.
- A second pass over its own output changes 112 more elements, and the curated
  library agrees with only 88% of those ends (of 42 it can place). So the stage runs
  once, and `--posthoc` refuses tables a pipeline run already re-boundaried.

Why candidates were not moved, at the defaults:

| reason | n |
|---|---|
| gate_identity | 563 |
| no_support | 335 |
| overlaps_element | 122 |
| host_exceeded | 7 |
| engulfs_element | 2 |
| length_filter | 2 |
| kmer2ltr_moved_ends | 1 |

**B1 (planted truncations, leave-one-out)**, on the 489 truth elements whose
TSD-backed ends the curated library confirms (4,588 planted obstacles):

| | pre-release | this design |
|---|---|---|
| exact or within 1 bp | 55.9% | 62.4% |
| over-extended | 0.15% | 0.15% |
| wrong | 0.44% | 0.48% |
| untouched controls changed | 0 of 489 | 5 of 489 |

Both designs used the same truth set: exact TSD, stranded, both LTRs ≥ 150 bp, at
most 20 per family, from families of ≥ 10 copies. `bench.py b1` as shipped draws truth
from families of any size, which scores lower on rice.

B1's truth is "the called ends carry an exact TSD". For 9 of the 514 truth elements
the curated library places the ends elsewhere (AT/GA-repeat flanks, chance or
duplicated TSDs), and this design's templates agree with the library there. Scored
against TSD truth alone, B1 therefore overstates this design's errors: 60.8% exact or
within 1 bp, 0.90% over-extended, 0.73% wrong and 12 of 514 controls changed
(pre-release: 54.7%, 0.67%, 0.46%, 4 of 514). Where a curated library exists, score
B1 on the truth elements it confirms.

## The three benchmarks

- **B1, planted obstacles.** Truth elements are calls with an exact TSD
  (`tsd_offset` `0,0`), a strand and both LTRs ≥ 150 bp, from labelled families, at
  most 20 per family (`--n`, default 1000). Each is copied onto a synthetic contig
  with 3 kb of its real flanks, and one obstacle is planted 20–400 bp from an outer
  end, in that end's LTR or at the matching place in its partner: none (a plain
  truncation), a 30 bp diverged patch, a 10–200 bp deletion, a 50 bp to 5 kb random
  insertion, or 300 bp or 1 kb of another family's element. The call is then cut
  where a finder would stop. Every truth element is excluded from the templates, so
  none templates its own copy. Each call is scored against its true ends: exact,
  within 1 bp, within 5 bp, over-extended (more than 5 bp past a true end), wrong, or
  unchanged.
- **B2, a real run.** The run's `_clean_` tables are copied under `<out>/run/` and
  re-bounded there by the full stage, including conflicts and merges; the run itself
  is only read. It reports moved, merged and trimmed counts and every rejection
  reason; TSD at the called, new and null ends (`tsd_gain`, `tsd_null`); the
  positional TSD spike (`spike_at_0` vs `spike_background`: an exact 5-bp TSD at the
  moved end vs 10–100 bp away); TG..CA termini before and after; move length; K2P
  shift; run time and peak memory; and TSD gain by clade and by age.
- **B3, false changes.** B1's untouched controls (a truth element called at its true
  ends) that were changed anyway (`false_change`). On a real run, B2 also reports the
  share of TSD-bearing calls that were moved (`b3_real_changed`), which should be
  zero, since those calls are never targets.

## Run it

```bash
cd <a finished LTRquest run directory>          # made without --reboundary
export PATH=/path/to/ltrquest/env/bin:$PATH      # blastn and makeblastdb on it
export PYTHONPATH=<repo>/src                     # or pip install . from this checkout
B=<repo>/benchmarks/reboundary/bench.py
P="A_LTRs B_LTRs"; G="A.fa B.fa"; T=<directory holding a Kmer2LTR/ checkout>

python $B b1 --run . --prefix $P --genome $G --tools-dir $T --out bench/b1_default --n 1000
python $B b2 --run . --prefix $P --genome $G --tools-dir $T --out bench/b2_default
python $B collect --out bench
```

Point `--run` at a run made without `--reboundary` (or one put back with
`--restore`); on re-boundaried tables, B2 would measure a second pass. `--tools-dir`
needs a Kmer2LTR that `ltrquest-reboundary` accepts (see the CHANGELOG). Each mode
writes `summary.json` into its `--out`. B1 adds `b1.tsv` (one row per contig) and
`synthetic.fa`; B2 adds the re-bounded copy under `run/` and every proposal in
`proposals.tsv`. `collect` tabulates every `b1_*` and `b2_*` directory under its
`--out`, naming each row by what follows the prefix.

`--set KEY=VALUE` (repeatable) sets any field of `reboundary.Settings` or
`ltr_place.Params`:

| key | flag | default | what it sets |
|---|---|---|---|
| `n_templates` | `--templates` | 5 | nearest templates placed per element |
| `min_support` | `--min-support` | 2 | agreeing placements a moved end needs |
| `agree` | | 2 | bp within which placements agree |
| `min_identity` | `--min-identity` | 0.8 | identity a placement's outer 30 bp need |
| `min_whole_identity` | | 0.6 | identity the whole template LTR needs |
| `min_ext` | `--min-ext` | 1 | the smallest outward move, bp |
| `max_trim` | `--max-trim` | 10 | the largest inward move, bp; 0 never trims |
| `anchor` | `--no-anchor` | true | search past a large indel |
| `anchor_len` | `--anchor-len` | 30 | template bases searched for there |
| `max_indel` | `--max-indel` | 5000 | how far past the call, bp |
| `max_templates_per_family` | | 100 | templates per family in the search; 0 = no cap |
| `merge_bp` | | 5 | how near a split call's end must be to the new end, bp |
| `blastn` | `--blastn` | `blastn` | the BLASTN executable; `makeblastdb` is taken from beside it, else `PATH` |
| `blast_task` | | `blastn` | its `-task` |

`threads` and `mutation_rate` are `Settings` fields too.

`tune.sh 1|2|3` sweeps the settings most likely to matter, running B1 (`--n 500`) and
B2 for each: `n_templates` × `min_support`, then `min_identity`, then `max_trim`, each
stage taking the previous one's choice as `KNOBS`.

## What these numbers do not show

- **TG..CA is not independent evidence.** A moved end is a template's end, so a
  template whose end sits on TG..CA passes the motif on. B2 reports it (`tgca_called`,
  `tgca_new`), but the TSD is the test that counts: it is read from the host flanks at
  the target's own locus, a different insertion from any template's.
- **B1 runs no conflict rules or merges.** Each synthetic contig holds one call, so a
  move that a neighbouring element would block on a real genome goes ahead in B1. Its
  false-change rate is, if anything, pessimistic.
- **B1's controls share its truth.** A control is a truth element called at its
  TSD-backed ends; where those ends are wrong (the 9 above), a correct move counts as
  a false change.
- **The curated check is not in the harness.** It needs a curated LTR library for the
  genome, placed on it independently of LTRquest.
