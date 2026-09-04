# LTRquest output reference

The full format of everything an LTRquest run writes. For what the tool is and
how to install it, see the [README](../README.md); for the Nextflow pipeline,
see [nextflow.md](nextflow.md).

---

## 1. Install

The container needs no root and no dependencies:

```bash
apptainer pull ltrquest.sif docker://ghcr.io/cwb14/ltrquest:1.0.1
```

`./ltrquest.sif` then takes the same flags as `ltrquest`, so read every
`ltrquest …` below as `./ltrquest.sif …`. Or install the tools onto your `PATH`
instead, one line at a time:

```bash
mamba env create -f environment.yml
mamba activate ltrquest
pip install .
```

See the [README](../README.md#install) for Docker and pip.

## 2. Quick start (test data)

Verify your installation against the *Arabidopsis* chromosome that ships with
the repo. It should finish in a few minutes.

```bash
mkdir athal_test && cd athal_test

ltrquest \
  --genome   ../tests/data/Athal_tair10_chr2.fa.gz \
  --proteins ../tests/data/Athal.pep.gz \
  --threads  20 --max-rounds 1
```

## 3. Usage

Most runs only need three flags:

| Flag | Meaning |
|---|---|
| `--genome` | Genome FASTA (`.fa`/`.fasta`, plain or `.gz`). **Required.** Accepts several, space-separated — see [§5](#5-multiple-genomes-shared-family-names). |
| `--proteins` | **Can come from any related species** — it does not have to match the genome. Recommended. |
| `--threads` | CPU threads (default 20). |

Occasionally useful:

| Flag | Meaning |
|---|---|
| `--max-rounds 1` | Run a single detection round. Use this if you don't care about nested elements. |
| `--out_prefix` | Output prefix (default: `<genome>_LTRs`). With several genomes it names the shared family namespace instead — see [§5](#5-multiple-genomes-shared-family-names). |
| `--terminate_count` | Stop iterating when a round finds fewer than this many elements (default 100). |
| `--run-sdust` | Drop candidates made mostly of low-complexity sequence, early. Off by default. |

## 4. Primary outputs: `depth<N>_clean_ltr.{tsv,fa}`

Elements are bucketed by nesting depth. `N` =
how many layers of LTR-RT are nested *inside* the element:

```
depth0:  [===============]                     no LTR-RT inside ("un-nested")
depth1:  [=====[depth0]=====]                  one LTR-RT nested inside
depth2:  [===[depth1 [depth0] ]===]            two layers nested inside
...
```

For the test run you'd get, e.g.:

```
Athal_tair10_chr2_LTRs_depth0_clean_ltr.tsv   ← table of un-nested elements
Athal_tair10_chr2_LTRs_depth0_clean_ltr.fa    ← their sequences
Athal_tair10_chr2_LTRs_depth1_clean_ltr.tsv   ← single-nested elements (ie, one LTR-RT nested inside)
Athal_tair10_chr2_LTRs_depth1_clean_ltr.fa
...
```

### 4.1 `depth<N>_clean_ltr.tsv` format

Tab-separated, one row per element, 33 named columns (`#`-prefixed header
line). Columns 1–29 are Kmer2LTR's own output, carried through unchanged:
Kmer2LTR runs ahead of classification, so its boundary calls, divergence
figures and TSDs are the round's single source of truth. LTRquest appends the
last four — `strand`, `family`, `domains`, `nest_status` — once
classification and clustering have run.

| # | Column | Meaning |
|---|---|---|
| 1 | `seq_id` | Element id: `chrom:start-end#Order/Superfamily/Clade`, e.g. `Chr2:102001-110500#LTR/Gypsy/Reina` |
| 2 | `seq_len` | Length of the record as supplied |
| 3 | `status` | `pass`, `weak_pair`, `k2p_undefined`, `no_pair`, `too_short`, or `all_ambiguous`. Only `pass` rows reach LTRquest's output |
| 4 | `ltr5_start` | 1-based inclusive, relative to the record. Always `1` here — see below |
| 5 | `ltr5_end` | 1-based inclusive, relative to the record: last bp of the 5′ LTR |
| 6 | `ltr3_start` | 1-based inclusive, relative to the record: first bp of the 3′ LTR |
| 7 | `ltr3_end` | 1-based inclusive, relative to the record. Always equals `seq_len` here — see below |
| 8 | `ltr5_len` | Length of the 5′ LTR |
| 9 | `ltr3_len` | Length of the 3′ LTR |
| 10 | `flank5_len` | bp of over-extension the detector added on the 5′ side. Always `0` here, by construction — see below |
| 11 | `flank3_len` | bp of over-extension the detector added on the 3′ side. Always `0` here, by construction — see below |
| 12 | `aln_len` | Length of the LTR-pair alignment, **including** gap columns |
| 13 | `n_sites` | Ungapped, unambiguous columns — the denominator for `identity`, `p_dist` and `k2p` |
| 14 | `n_ts` | Transitions |
| 15 | `n_tv` | Transversions |
| 16 | `n_gapcols` | Alignment columns with a gap on either side |
| 17 | `identity` | `n_match / n_sites` |
| 18 | `p_dist` | `(n_ts + n_tv) / n_sites` |
| 19 | `k2p` | Kimura 2-parameter distance. `NA` when saturated, never clamped |
| 20 | `k2p_se` | Standard error of `k2p` |
| 21 | `bitscore` | Global alignment score of the two LTRs under a fixed generic model |
| 22 | `flank_margin_bits` | How decisively the boundary was called. `NA` when the discovery core already reached both termini |
| 23 | `cigar` | Extended CIGAR of the LTR-pair alignment; query = 5′ LTR |
| 24 | `motif` | The two terminal dinucleotides, lowercased, e.g. `tg...ca`. Read off the called boundary, never used to find it |
| 25 | `k2p_time` | `round(k2p / (2 × mutation_rate))`, years since insertion; set by `--mutation-rate` (default `3e-8`) |
| 26 | `orientation` | Whether *this record* is stored reverse-complemented relative to its own header locus. Not biology — see `strand` below |
| 27 | `tsd` | Target-site duplication at the called boundary. `.` = searched and absent; `NA` = could not be searched |
| 28 | `tsd_offset` | `d5,d3` — how far each boundary had to move for `tsd` to appear, positive meaning *into* the element. `NA` when `tsd` is `.` |
| 29 | `tsd_input` | The same measurement at the record's termini as originally supplied, before re-bounding |
| 30 | `strand` | The element's biological orientation, resolved by LTRquest from classification, domain order and pass-2 homology: `+`, `-`, or `.` if none of those resolved it |
| 31 | `family` | Family id from the pooled consensus-LTR clustering, `<prefix>_fam00001` … — see section 5 |
| 32 | `domains` | Protein domains with genomic coordinates: `GENE\|clade@start-end;...`, or `.` |
| 33 | `nest_status` | `nest-outer:chrom:s-e` (that element is inside me) and/or `nest-inner:chrom:s-e` (I am inside that element), `;`-joined, or `.` if un-nested |

Three things this table will burn you on if you skim it:

- **`.` and `NA` are not the same value.** `tsd` (and `tsd_input`, the same
  measurement taken at the record's original termini) write `.` when the
  search ran and found no duplication — a real negative result — and `NA`
  when the search had no chance to run at all, for instance because the
  element's contig could not be matched back to the reference genome. Treat
  the two as one "missing" bucket and you will misjudge every TSD call:
  presence-vs-absence and ran-vs-could-not-run are different findings, and
  only one of them says anything about the boundary.
- **`strand` and `orientation` answer different questions.** `orientation`
  is a storage fact: whether this FASTA record happens to be the reverse
  complement of the genome at its own header coordinates, decided purely by
  Kmer2LTR comparing the record against `--genome`. `strand` is a biological
  claim: which way the element is transcribed, decided later by LTRquest from
  classification, domain order and pass-2 homology. The two are independent —
  an element stored forward (`orientation=+`) can carry either strand, and
  vice versa.
- **Every row satisfies `ltr5_start == 1` and `ltr3_end == seq_len`.** That is
  what "Kmer2LTR-bounded" means: the record has already been cut down to the
  called LTR pair, so both flanks report `0` and there is nothing upstream or
  downstream of the element left in the record. This is not merely intended —
  `assert_bounded` checks it against every row at the end of each round and
  raises before the round can finish if a single one fails.

### 4.2 `depth<N>_clean_ltr.fa` format — nested regions are masked

Each record is the **full-length LTR-RT**, with every
nested inner element **hard-masked** by an IUPAC letter chosen by the *inner
element's own depth*:

| Inner element's depth | Mask letter |
|---|---|
| 0 | `N` |
| 1 | `R` |
| 2 | `D` |
| 3+ | `Y`, `S`, `W`, `K`, `M`, `B`, `H` |

So a `depth1` record looks like `ACGT...NNNN...ACGT` (its depth0 insert masked
with `N`), and a `depth2` record looks like `ACGT...RRR..NNN..RRR...ACGT`
(its depth1 child masked `R`, whose own depth0 child is masked `N`). 

**To remove the masked inner sequence** and keep only the outer element's sequence:

```bash
awk '/^>/{if(s)print s; s=""; print; next}{gsub(/[^ACGTacgt]/,""); s=s $0}END{if(s)print s}' \
  Athal_tair10_chr2_LTRs_depth1_clean_ltr.fa > depth1_removed_nested.fa
```

`depth0` records contain no masked nest-ins, so this is only needed for depth ≥ 1.

## 5. Multiple genomes: shared family names

Pass several genomes to give closely-related species **one family vocabulary**,
which is what makes between-species family comparisons meaningful:

```bash
ltrquest \
  --genome   ../tests/data/Athal_tair10_chr2.fa.gz Alyrata_chr2.fa \
  --proteins ../tests/data/Athal.pep.gz \
  --threads  40
```

Each genome is detected on its own, then **all** detected elements are pooled
for a single clustering pass, so `family` means the same thing everywhere:

```
Athal_tair10_chr2_LTRs_depth0_clean_ltr.tsv     family = merged_fam00001
Alyrata_chr2_LTRs_depth0_clean_ltr.tsv          family = merged_fam00001
Athal_tair10_chr2_LTRs_all_depth_LTR_cleaned.gff3
Alyrata_chr2_LTRs_all_depth_LTR_cleaned.gff3
merged_all_ltr.consensus_id0.75_cluster.tsv     ← the one shared table
```

False-positive families are also called over the pooled set rather than one
species at a time, so a repeat that looks convincing in isolation but wrong
across species gets caught.

Naming: every genome keeps its own `<basename>_LTRs` prefix. `--out_prefix`
names only the shared pool — `--out_prefix Arabidopsis` gives
`Arabidopsis_fam00001` and `Arabidopsis_all_ltr.*`, default `merged`. With a
single genome the two are the same thing, so nothing changes.

One `--proteins` file serves every genome.

> **Sequence IDs must be unique across genomes.** Two files both calling a
> chromosome `Chr2` are rejected before any work starts: pooled clustering keys
> elements on `chrom:start-end`, so a shared ID would cross-assign families
> *and* cross-purge real elements between species. Rename first:
>
> ```bash
> awk '/^>/{sub(/^>/,">Aly_")}1' Alyrata.fa > Alyrata.renamed.fa
> ```

## 6. GFF3 annotation

Two files, pooled across all depths and built from the FP-purged
`_clean_ltr.tsv` set (falling back to the raw set, with a warning, if the FP
stage never ran):

| Output | Contents |
|---|---|
| `<prefix>_all_depth_LTR_cleaned.gff3` | The LTR-RTs |
| `<prefix>_all_depth_protein_LTR_cleaned.gff3` | The same, plus every miniprot protein alignment. Only written when the run had `--proteins` |

Each element is one block:

```
chr  LTRquest  LTR_retrotransposon   13031  17307  .  -  .  ID=…_LTRRT_00001;Name=chr:13031-17307;
                                                            classification=LTR/Gypsy/Tekay;superfamily=Gypsy;clade=Tekay;
                                                            family=…_fam00001;family_size=5;depth=0;
                                                            K2P_d=0.026305;K2P_T=438418;strand_source=tesorter
chr  LTRquest  long_terminal_repeat  13031  14048  .  -  .  ID=…_LTRRT_00001.lLTR;Parent=…_LTRRT_00001
chr  LTRquest  long_terminal_repeat  16290  17307  .  -  .  ID=…_LTRRT_00001.rLTR;Parent=…_LTRRT_00001
chr  LTRquest  protein_match         14482  14712  .  -  .  ID=…_LTRRT_00001.CHD.1;Parent=…;Name=CHD;clade=Tekay
###
```

Notes:

- Nested elements are **flat top-level features**, not children of their host:
  GFF3 `Parent` means part-of, and a nested LTR-RT is not a part of the element
  it landed in. Use the `depth` and `nest_status` attributes instead.
- Blocks are coordinate-sorted, and lines within a block are too, but
  **overlapping blocks are not interleaved** — so a nested element's block
  follows its host's in full.
- GFF3 attribute names are their own vocabulary, independent of the element
  table's column names: `K2P_d`/`K2P_T` above read from `k2p`/`k2p_time` in
  section 4.1. Keeping the attribute names stable regardless of what the
  table's own columns are called is what lets an existing GFF3 consumer keep
  working unmodified.

## 7. Plots (`<prefix>_plots/`)

| Output | What it shows |
|---|---|
| `struct/<clade>_average.pdf` | Average element structure per classification clade with bootstrap 95% CI whiskers on feature positions |
| `struct/<clade>_individual.pdf` | Every element of the clade drawn individually, all domains shown |
| `struct/all_elements.pdf` | All elements on one page |
| `<prefix>_summary.pdf` | Multi-page summary |
| `<prefix>_TEGV.html` | Self-contained interactive genome browser (open in web browser) |

## 8. Benchmarks

On a simulated genome (PrinTE) with 4,468 true intact LTR-RTs, 20 threads,
scored at ≥90% reciprocal overlap:

| Tool | Runtime | TP | FP | FN | F1 |
|---|---|---|---|---|---|
| **LTRquest** | **5:18** | 4,328 | 170 | 140 | **0.966** |
| EDTA raw LTR module (LTR_retriever) | 20:28 | 2,085 | 0 | 2,383 | 0.64 |
