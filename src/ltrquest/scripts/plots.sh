#!/usr/bin/env bash
# plots.sh — post-completion plotting driver for the ltrquest pipeline.
#
# Discovers the final annotation files an ltrquest run left behind (preferring
# the FP-purged *_depth<N>_clean_ltr.* set), derives the small inputs the
# plotters need (chrom sizes, combined TSV), and runs:
#   1. ltrquest.plot_struct   — per-family structure PDFs, drawn both with
#                               nested insertions in place and excised
#   2. ltrquest.plot_summary  — multi-page summary PDF (K2P, karyotype, sizes)
#   3. ltrquest.tegv          — interactive HTML viewer
# Each plotter is independent: a failure prints [WARN] and the others still
# run. Exit code = number of failed plotters (0 = all good). Missing/invalid
# inputs die immediately with a clear message.
#
# Also runnable standalone on any existing ltrquest output directory.
set -euo pipefail

# Python steps run as modules of the installed package, so this behaves the same
# from `pip install .`, from a source checkout with PYTHONPATH=src, and in the
# container. The `ltrquest` console script exports LTRQUEST_PYTHON as the
# interpreter that owns the package; running this file directly falls back to
# whatever python3 is on PATH.
PY="${LTRQUEST_PYTHON:-python3}"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[ltrquest-plots] $*"; }
warn() { echo "[WARN] $*" >&2; }
vlog() { [[ "$VERBOSE" == true ]] && echo "  \$ $*" || true; }

abspath() ( cd "$(dirname "$1")" && printf '%s/%s\n' "$PWD" "$(basename "$1")" )

usage() {
  cat <<'EOF'
Usage: plots.sh --prefix PREFIX --genome GENOME[.gz] [options]

Required:
  --prefix PREFIX      Output prefix used by the ltrquest run
                       (e.g. Athal_tair10_chr2_LTRs).
  --genome FASTA       Genome FASTA (plain or .gz). Streamed once for
                       chromosome sizes; never loaded into memory.

Options:
  --indir DIR          Directory holding the wrapper outputs (default: .).
  --out-dir DIR        Plot output directory (default: <indir>/<prefix>_plots).
  --nesting MODE       Structure PDFs: expanded | collapsed | both
                       (default: both -- writes <family>_individual_expanded.pdf
                       and <family>_individual_collapsed.pdf; any other mode
                       writes a single <family>_individual.pdf).
  --k2p-mode MODE      Summary PDF K2P panels: hist | kde (default: hist,
                       i.e. binned density bars rather than a smoothed curve).
  --k2p-bin-width F    Bin width for --k2p-mode hist (default: 0.0005).
  --family-k2p-xmax F  Lock family-page K2P x-axes to [0, F]. Default: the
                       99th percentile of the pooled K2P.
  -v, --verbose        Echo each plotter command before running it.
  -h, --help           This help.
EOF
}

PREFIX="" GENOME="" INDIR="." OUT_DIR="" VERBOSE=false
NESTING="both" K2P_MODE="hist" K2P_BIN_WIDTH="" FAMILY_K2P_XMAX=""
[[ $# -eq 0 ]] && { usage; exit 1; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="${2:-}"; shift 2;;
    --genome) GENOME="${2:-}"; shift 2;;
    --indir) INDIR="${2:-}"; shift 2;;
    --out-dir) OUT_DIR="${2:-}"; shift 2;;
    --nesting) NESTING="${2:-}"; shift 2;;
    --k2p-mode) K2P_MODE="${2:-}"; shift 2;;
    --k2p-bin-width) K2P_BIN_WIDTH="${2:-}"; shift 2;;
    --family-k2p-xmax) FAMILY_K2P_XMAX="${2:-}"; shift 2;;
    -v|--verbose) VERBOSE=true; shift;;
    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -n "$PREFIX" ]] || die "--prefix is required"
case "$NESTING" in expanded|collapsed|both) ;; *) die "--nesting must be expanded, collapsed or both (got: $NESTING)";; esac
case "$K2P_MODE" in hist|kde) ;; *) die "--k2p-mode must be hist or kde (got: $K2P_MODE)";; esac
[[ -n "$GENOME" ]] || die "--genome is required"
[[ -f "$GENOME" ]] || die "Genome not found: $GENOME"
[[ -d "$INDIR" ]] || die "Input directory not found: $INDIR"
"$PY" -c 'import ltrquest' 2>/dev/null \
  || die "the ltrquest package is not importable by '$PY'; pip install ltrquest, or point LTRQUEST_PYTHON at an interpreter that has it"
INDIR="$(abspath "$INDIR")"
GENOME="$(abspath "$GENOME")"
[[ -z "$OUT_DIR" ]] && OUT_DIR="${INDIR}/${PREFIX}_plots"
mkdir -p "$OUT_DIR"
OUT_DIR="$(abspath "$OUT_DIR")"

# ---- Discover annotation set: clean preferred, raw fallback -------------
# Depth files are probed numerically (wrapper hard-caps rounds at 10; 20 is
# comfortable headroom) so ordering is depth0, depth1, ... regardless of glob
# lexicography.
depth_tsvs=() variant=""
for v in clean raw; do
  found=()
  for ((d = 0; d <= 20; d++)); do
    if [[ "$v" == clean ]]; then f="${INDIR}/${PREFIX}_depth${d}_clean_ltr.tsv"
    else f="${INDIR}/${PREFIX}_depth${d}_ltr.tsv"; fi
    [[ -s "$f" ]] && found+=("$f")
  done
  if ((${#found[@]} > 0)); then depth_tsvs=("${found[@]}"); variant="$v"; break; fi
done
((${#depth_tsvs[@]} > 0)) \
  || die "No ${PREFIX}_depth<N>[_clean]_ltr.tsv files found in ${INDIR} — is this a completed ltrquest output directory?"
[[ "$variant" == raw ]] \
  && warn "no *_clean_ltr.tsv files found; falling back to raw depth files (FP families not purged)"
log "annotation set (${variant}): $(basename -a "${depth_tsvs[@]}" | tr '\n' ' ')"

# Matching FASTAs (TEGV needs them; a missing one drops that depth from TEGV only).
depth_fas=()
for t in "${depth_tsvs[@]}"; do
  fa="${t%.tsv}.fa"
  if [[ -s "$fa" ]]; then depth_fas+=("$fa")
  else warn "no FASTA next to $(basename "$t"); TEGV will omit that depth"; fi
done

# ---- Derived inputs -----------------------------------------------------
log "writing ${PREFIX}.chrom.sizes (streaming genome pass)"
zcat -f -- "$GENOME" | awk '
  /^>/ { if (name != "") print name "\t" len; name = substr($1, 2); len = 0; next }
  { gsub(/[ \t\r]/, ""); len += length($0) }
  END { if (name != "") print name "\t" len }
' > "${OUT_DIR}/${PREFIX}.chrom.sizes"
[[ -s "${OUT_DIR}/${PREFIX}.chrom.sizes" ]] || die "Genome produced no sequences: $GENOME"

combined="${OUT_DIR}/${PREFIX}_combined_ltr.tsv"
head -n 1 "${depth_tsvs[0]}" > "$combined"
for t in "${depth_tsvs[@]}"; do grep -v '^#' "$t" >> "$combined" || true; done

# ---- Run the three plotters (independent; warn-on-fail) -----------------
FAILED=0
run_plotter() {
  local name="$1"; shift
  log "running ${name}..."
  if "$@"; then log "${name}: OK"
  else warn "${name} failed; continuing with the remaining plotters."; FAILED=$((FAILED + 1)); fi
}

plot_struct() {
  local args=(--prefix "$PREFIX" --indir "$INDIR" --out_dir "${OUT_DIR}/struct" \
              --depth all --nesting "$NESTING")
  vlog "$PY" -m ltrquest.plot_struct "${args[*]}"
  "$PY" -m ltrquest.plot_struct "${args[@]}"
}

plot_summary() {
  local args=(--species "$PREFIX" --fai-suffix .chrom.sizes \
              --aln-suffix _combined_ltr.tsv --outpdf "${PREFIX}_summary.pdf" \
              --k2p-mode "$K2P_MODE")
  [[ -n "$K2P_BIN_WIDTH" ]] && args+=(--k2p-bin-width "$K2P_BIN_WIDTH")
  [[ -n "$FAMILY_K2P_XMAX" ]] && args+=(--family-k2p-xmax "$FAMILY_K2P_XMAX")
  local cluster="${INDIR}/${PREFIX}_all_ltr.consensus_id0.75_cluster.tsv"
  if [[ -s "$cluster" ]]; then
    args+=(--mmseqs-tsv "$cluster" --family-min-count 10)
  else
    log "no consensus cluster TSV; skipping per-family pages"
  fi
  vlog "(cd ${OUT_DIR} &&) $PY -m ltrquest.plot_summary ${args[*]}"
  (cd "$OUT_DIR" && "$PY" -m ltrquest.plot_summary "${args[@]}")
}

plot_tegv() {
  local args=(--te-fastas "${depth_fas[@]}" --ltr-divergence "${depth_tsvs[@]}" \
              --output "${OUT_DIR}/${PREFIX}_TEGV.html")
  shopt -s nullglob
  local genic=("${INDIR}/${PREFIX}"_r1.work/*.genic.gff)
  shopt -u nullglob
  if [[ ${#genic[@]} -eq 1 && -s "${genic[0]}" ]]; then
    log "gene track: ${genic[0]##*/} (miniprot alignments; run with --proteins to produce this)"
    args+=(--gff3 "${genic[0]}")
  else
    log "no single r1 genic.gff found; TEGV runs without a gene track"
  fi
  vlog "$PY" -m ltrquest.tegv "${args[*]}"
  "$PY" -m ltrquest.tegv "${args[@]}"
}

run_plotter "plot_struct" plot_struct
run_plotter "plot_summary (summary PDF)" plot_summary
if ((${#depth_fas[@]} > 0)); then
  run_plotter "TEGV" plot_tegv
else
  warn "no depth FASTAs available; skipping TEGV"
  FAILED=$((FAILED + 1))
fi

log "done: ${OUT_DIR} (${FAILED} plotter failure(s))"
exit "$FAILED"
