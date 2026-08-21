# Changelog

All notable changes to LTRquest are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

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
