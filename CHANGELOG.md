# Changelog

All notable changes to LTRquest are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Fixes the three ways a fresh install could fail on a machine that is not a
developer workstation. All three were found by installing from scratch — pulling
the published image and creating the conda environment — rather than by testing
the source tree.

Upgrades the bundled Kmer2LTR from its original k-mer/MAFFT implementation to
the rewritten package, and reorders the pipeline around what that rewrite
makes possible. Kmer2LTR now runs *ahead of* classification instead of after
it, so its boundary calls, K2P divergence figures and TSDs are the round's
single source of truth rather than figures LTRquest re-derived and could
disagree with. The new tool's CLI and output schema are not compatible with
the old one's, which forces the schema, flag and stage changes below.

### Added

- `--mutation-rate`, forwarded to Kmer2LTR's `-u` for both the per-round
  table and the pooled family-clustering pass. Insertion age used to come
  from a fixed μ = 3×10⁻⁸ baked into the old Kmer2LTR script with no flag to
  change it; now `k2p_time` is `NA` unless you supply a rate that actually
  matches your organism.

### Changed

- `ltrquest-plots`' `--min_ltr_aln` is renamed `--min_ltr_len`, and
  `flagged_false_positives.tsv`'s `aln_len_raw` column is renamed
  `ltr3_len_raw`. Both hold a 3' LTR length, not an alignment length; the
  names say what they measure.
- `ensure_trfmod` and `ensure_sdust` look in `--tools-dir`, then on `PATH`, then
  clone and build — the order `ensure_tools` has used since 1.0.1. `sdust` is
  packaged (`bioconda::sdust`) and is now in `environment.yml`, so `--run-sdust`
  no longer needs a compiler.
- Failures to find or build a helper name the conda package that provides it and
  the compiler toolchain that would build it, instead of only the tool.
- LTRquest's own `depth<N>_clean_ltr.tsv` is now Kmer2LTR's 33-column schema
  (see [docs/outputs.md](docs/outputs.md)), replacing the old 21-column one
  built around `LTR_len`, `raw_d`/`raw_T`, `JC69_d`/`JC69_T` and
  `left_trim`/`right_trim` — none of which the new boundary model or K2P
  estimator produces. `K2P_d`/`K2P_T` survive only as the GFF3 attribute
  names they always were; the table columns underneath them are now
  `k2p`/`k2p_time`.
- `flag_fp_families.py` moves out of the Kmer2LTR clone and into LTRquest as
  `ltrquest.flag_fp` / `ltrquest-flagfp`. Flagging false-positive families
  reads Kmer2LTR's clustering output; it was never a Kmer2LTR concern, only
  ever colocated with it because that happened to be where the script lived.

### Removed

- Six flags that had nothing left to forward to once Kmer2LTR's internals
  changed underneath them: `--wfa-align` (the new tool always aligns with
  WFA, so there is no MAFFT default left to opt out of), `--kmer2ltr-max-win-
  overdisp` and `--kmer2ltr-min-retained-fraction` (tuning for the old k-mer
  window search; the new discovery stage is exact Smith-Waterman with no
  overdispersion or retained-fraction knobs), `--kmer2ltr-domains` (the old
  script's own domain calls; domains have only ever come from TEsorter2 in
  this pipeline), `--tsd-min-len` (the new Kmer2LTR accepts a TSD on an
  exact-match, two-distinct-base rule with no length floor to set), and
  `--tesorter-use-ret` (chose whether classification saw the ret span or the
  internal span of a candidate; structurally inert now that Kmer2LTR needs
  both LTRs to measure divergence, so classification must always receive the
  full-length element and the internals-only path can no longer feed it).
- `ltrquest-reconcile --scn`. It fed a shared-LTR lookup keyed on each round's
  pre-trim SCN locus, but Kmer2LTR renames every element to its post-trim
  locus before reconcile ever sees it, so the flag could no longer reach the
  elements it was meant to describe. Boundaries are read from each round's
  own `--tsv` table instead (`ltr_bounds_from_table`), which carries the
  right key and the right coordinates on the same row.

### Fixed

- **`./ltrquest.sif --genome x.fa` did not work.** The image declared
  `ENTRYPOINT []` with only a `CMD`, and Apptainer builds a SIF's runscript from
  those two: with no ENTRYPOINT it *replaces* the command with the user's
  arguments, so the first flag was run as a program —
  `FATAL: "--help": executable file not found in $PATH`. The image now has an
  entry point ([`bin/entrypoint.sh`](bin/entrypoint.sh)) that runs a first
  argument naming a command and treats anything else as `ltrquest` arguments.
  `apptainer exec IMG ltrquest …` and `docker run IMG ltrquest …` are unchanged
  — `exec` never reaches an entry point, and an explicit command still wins.
- **A working directory reached through a symlink left the run in `$HOME`.**
  Apptainer mounts your home and working directories and nothing else. Where
  home points elsewhere — `/home/you/data -> /scratch/you`, the usual cluster
  layout — it adds no mount for the target, cannot follow the link from inside,
  and falls back to `$HOME` with a warning that scrolls past. Every relative
  path in the command then resolved somewhere else, and the only symptom was
  `ERROR: Genome not found: your_genome.fa`. Two changes: a missing input inside
  a container now says so and names the fix, and
  [`bin/ltrquest-container`](bin/ltrquest-container) — shipped inside the image
  — works the mounts out from the command it is given.
- **The conda environment could not build the helpers it needs.**
  `environment.yml` had no `git`, `make` or compiler, but TRF-mod is compiled on
  first use and runs by default, so a clean environment died on `cc: not found`
  the first time it saw a genome. This was invisible on any machine with system
  build tools. The environment now carries them, and TRF-mod's build is retried
  with `CC=$CC` because its makefile hardcodes `CC=gcc`, which a conda toolchain
  does not provide.
- **GFF3 `long_terminal_repeat` sub-features landed on the wrong span for
  minus-strand elements.** This predates the rewrite: `gff3.py` applies
  `start + ltr5_end - 1` as a forward-genomic offset, but the pipeline used to
  measure those coordinates against a library where minus-strand records had
  already been reverse-complemented for classification, so the two
  sub-features were mis-sized whenever the LTRs differed in length or the
  flanks were asymmetric. It is fixed as a consequence of Kmer2LTR now
  running on the original forward-genomic sequence, ahead of classification,
  rather than of anything `gff3.py` itself had to change.
- **Same-round nest masking painted the wrong bases on reverse-complemented
  elements.** `nest_status` intervals are always forward-genomic, but on the
  ~45% of elements whose library record Kmer2LTR reverse-complemented for
  family clustering, the interval was never mirrored before being painted
  onto the flipped record. `mask_same_round_inners_in_fa` now mirrors it for
  every record `bounded_fasta_oriented` flipped.
- **A truncated or extended element could be recorded as a nested insertion.**
  `_ltrs_shared` tells a truncation/extension variant from a genuine nested
  insertion by checking whether the two elements' LTRs overlap, and returns
  `False` ("distinct LTRs", i.e. a real nest) whenever either element's key is
  absent from the LTR-boundary map it is given. Reconcile's map came from
  `--scn`, keyed on each round's pre-trim locus, and was checked against
  elements Kmer2LTR had already renamed to their post-trim locus, so a missing
  key read as nesting instead of as a failed lookup -- disabling the
  truncation guard for 38% of elements (79/206 survivors on a full
  Arabidopsis chr2 run). `dedup_kmer2ltr_tsv` made the same call from a
  boundary map re-keyed through the rename rather than read from the
  post-trim table. Both now derive their map with `ltr_bounds_from_table`,
  which reads each element's `ltr5_start`/`ltr5_end`/`ltr3_start`/`ltr3_end`
  off its own row, so the key and the bounds it is checked against can no
  longer fall in different frames.

## [1.0.1] - 2026-08-24

Fixes two bugs that made the published 1.0.0 container unusable for a real run.
Both were found by pulling the released image and running it, rather than by
testing the source tree.

### Fixed

- **The container could not run a genome at all.** `ensure_tools` cloned and
  compiled minimap2 and miniprot from source *before* it would consider the
  copies already on `PATH`. The runtime image ships no `git` — deliberately, so
  a run never reaches the network — so detection died on
  `FileNotFoundError: 'git'` before it processed a single base. Resolution order
  is now: a prebuilt copy in `--tools-dir`, then `PATH`, then clone and build.
  Conda and pip users benefit too: they no longer pay a source build for tools
  their environment already provides.
- **`-profile singularity` could not pull its image.** Every Nextflow module
  pointed Singularity and Apptainer at `oras://…:1.0.0-singularity`, an artifact
  the release workflow never publishes, so the pull 404'd. Both engines convert
  an OCI image on the fly, so all engines now use the one published image.
- A missing `git` now raises a message naming the tool, the conda command that
  installs it, and the `--tools-dir` alternative, instead of a bare traceback
  from six frames inside `subprocess`.

### Added

- The container smoke test now resolves the helper binaries offline, which is
  the step `--help` never reached and the reason 1.0.0 shipped broken.
- Tests pinning the helper-resolution order, and guards against the package
  version drifting out of sync with `nextflow.config` or the module container
  tags, and against any module referencing an unpublished ORAS artifact.

### Changed

- README leads with Apptainer/Singularity rather than Docker, and gives
  copy-pasteable commands for the bundled test data on every install route.
  The audience is biologists on shared clusters, who cannot install Docker and
  do not have root.

## [1.0.0] - 2026-08-21

First release as a standalone tool. LTRquest was extracted from
[synLTR](https://github.com/cwb14/synLTR) `module2` with `git filter-repo`, so
the development history of every file it carries came with it.

### Added

- `ltrquest` console entry point, plus one per pipeline stage
  (`ltrquest-detect`, `-mask`, `-reconcile`, `-annotate`, `-gff3`, `-plots`).
- Installable Python package (`pip install .`), with the pipeline driver shipped
  as package data rather than a loose script.
- Nextflow DSL2 pipeline (`main.nf`) following nf-core conventions, with `stub:`
  blocks on every process so the whole DAG can be validated without the external
  toolchain.
- Container image with the helper repositories (Kmer2LTR, TEsorter2, TRF-mod,
  sdust) pre-built at pinned commits, so a run never reaches the network.
- `LTRQUEST_TOOLS_DIR` to point a run at a pre-built helper-tool directory.
- pytest suite covering the masking, reconciliation, annotation and GFF3 stages,
  and the installed package's wiring.
- GitHub Actions running the tests, lint, a Nextflow stub run, and a container
  build on every push.

### Changed

- Filenames collapsed to module names: `ltrharvest5.py` → `ltrquest.detect`,
  `mask_ltr.py` → `ltrquest.mask`, `reconcile_nests.py` → `ltrquest.reconcile`,
  `ltr_annotate.py` → `ltrquest.annotate`, `ltr_tsv_to_gff3.py` →
  `ltrquest.gff3`, `ltrharvest_wrapper2.sh` → the `ltrquest` driver.
- `--ltrharvest5-args` / `--ltrharvest5-args-from-round` renamed to
  `--detect-args` / `--detect-args-from-round`.
- `--script_path` removed. Stages are resolved as modules of the installed
  package, so there is nothing left to point at.
- GFF3 column 2 now reads `LTRquest` rather than `synLTR`.

### Removed

- The superseded detector and phasing revisions, the orphaned utilities, and the
  `solo_LTR` pipeline, none of which the entry point reached. They remain in the
  synLTR history.
