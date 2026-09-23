# Re-boundarying benchmarks

How well `--reboundary` recovers truncated LTR-RT ends, and how to re-measure it.

## The problem it is measuring

LTRharvest and LTR_FINDER stop extending at the first obstacle between an element's
two LTRs, so some calls are short at one or both ends. `--reboundary` extends them to
the ends the rest of the family agrees on. The question is whether those new ends are
right.

The honest check is the **target-site duplication**: a transposon insertion duplicates
5 bp of host sequence at the insertion point, so a correct boundary lands on a TSD and
a wrong one does not. Nothing in the algorithm looks at TSDs when placing an end, so
TSD recovery is independent evidence.

## How it performs

Four *Poa* genomes, 208,018 elements, at the default settings:

```
64,181 candidates -> 20,925 extended
TSD at the new boundary   0.565    (displaced-flank null 0.010)
  cross-validated         0.535    (null 0.009, n=13,170)
TG..CA termini            0.317 -> 0.610
extension length          median 323 bp, q90 1,269 bp
K2P shift                 median +0.0018
cost                      ~26 min, 4 GB  (~80 s with a cached family phase)
```

TSD recovery falls with element age, which is what makes the result believable — TSDs
decay by mutation, and nothing in the method knows an element's age:

| K2P | < 0.005 | 0.005–0.02 | 0.02–0.05 | ≥ 0.05 |
|---|---|---|---|---|
| TSD recovered | 0.666 | 0.676 | 0.569 | 0.338 |

Against planted obstacles with known truth (B1, n=3000), it recovers the true end
exactly or within 1 bp **40%** of the time, over-extends 2%, and gets it wrong 1.6%.
Recovery depends on what it has to cross: a 30 bp diverged patch 52%, a 50 bp deletion
47%, a 300 bp insertion 40%, a 5 kb insertion 26%. Large insertions are where it
struggles.

Why elements are rejected, at the defaults:

| reason | n |
|---|---|
| gate_identity | 31,816 |
| family_qc_failed | 3,804 |
| overlaps_element | 3,741 |
| kmer2ltr_not_pass | 3,339 |
| host_exceeded | 274 |
| engulfs_element | 149 |
| kmer2ltr_reverted | 98 |
| length_filter | 29 |
| duplicates_host | 6 |

## The three benchmarks

- **B1 — planted obstacles.** Well-called elements are copied onto synthetic contigs
  with 3 kb of their real flanks, one obstacle is planted in one LTR (a 30 bp diverged
  patch, a 10–200 bp deletion, or a 50 bp–5 kb insertion), and the call is cut where
  the finders would stop. Each truth element is left out of its own family's model.
  Scored against the known true end.
- **B2 — a real run.** The clean tables are copied and re-bounded. Measures TSD
  recovery against a null with the flank displaced 1 kb, plus K2P shifts, time, memory.
  The `cv` variant picks family QC on half of each family's candidates and scores it on
  the other half, so TSD-based QC cannot inflate its own validation.
- **B3 — false changes.** Elements that should not have moved, but did.

## Run it

```bash
cd <a finished LTRquest run directory>
export PATH=/path/to/ltrquest/env/bin:$PATH
B=<repo>/benchmarks/reboundary/bench.py
P="A_LTRs B_LTRs"; G="A.fa B.fa"; T=<kmer2ltr checkout>

python $B b1 --run . --prefix $P --genome $G --tools-dir $T --out bench/b1 --n 1000
python $B b2 --run . --prefix $P --genome $G --tools-dir $T --out bench/b2 --cache bench/c.pkl
python $B collect --out bench
```

`--set KEY=VALUE` reaches any field of `reboundary.Settings` or `ltr_place.Params`.
**Always pass `--cache`** — the family phase is ~25 min, the rest is ~1 min — and give
each `method`/`references` combination its own cache file.

`tune.sh 1|2|3|4a|4b` re-runs the four-stage sweep that set the defaults below.

## Defaults, and why

| knob | value | why |
|---|---|---|
| `method` | nearest | beats `consensus` on recovery, 0.404 vs 0.394 at n=3000 (p=0.03) |
| `references` | modal | ties with `tsd`; keeps TSD out of the metric that validates it |
| `min_identity` | 0.8 | 0.7 buys +1 point of recovery and doubles the wrong rate |
| `anchor_len`, `max_indel` | 30, 5000 | no measurable effect within noise |
| `credit` | 5000 | Kmer2LTR's veto mostly overturns *correct* extensions; see below |
| `max_ratio`, `min_copies` | 1.15, 10 | inert below those values |
| `qc_min_n`, `untested` | 5, accept | all settings within noise |

`credit` is the one worth understanding. It is how much external evidence Kmer2LTR is
told the family model represents. At `credit=0` it reverts 8,698 proposals; at 5000,
98. Recovery rises and the B1 error rate *falls* as the veto is withdrawn, so those
reversions were mostly wrong. Kmer2LTR still does the re-scoring and owns every column
it always did — only its veto turns out to be redundant here.

## What these numbers do not show

- **TG..CA is partly self-fulfilling.** A template is cut from a reference's own
  *called* span, so a reference whose original call happened to land on TG..CA passes
  the motif to elements placed against it, within the ~2 bp matching window. The
  improvement is real but not independent. TSD recovery is not affected: a TSD is read
  from the host flanks at the target's own locus — a different insertion event than any
  element that built the model.
- **B1's false-change rate is pessimistic.** Its synthetic contigs have no neighbouring
  annotations, so the conflict rules that block a kb-scale jump on real data cannot fire
  there; and a control whose original call was already short counts as a false change
  even when the extension is correct.
- **`cv TSD` is blind to anything after arbitration.** It is computed from the family
  phase's proposals, so it is identical at every value of `credit`. Read `tsd_gain` and
  the B1 columns for those knobs instead.
