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
(212,360 elements). Kmer2LTR aa25f46. Tables are pasted from `bench.py collect`.
