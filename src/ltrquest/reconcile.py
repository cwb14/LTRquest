#!/usr/bin/env python3
"""Reconcile nested-LTR-RT calls across multiple rounds of ltrquest.detect
into depth-bucketed libraries.

Each surviving element is assigned an inward-chain depth:
    chain_inward(x) = 0 if no LTR-RT is strictly inside x, else
                      1 + max(chain_inward(direct children of x)).

Direct children: y is a direct child of x iff x strictly contains y AND
no other z in the pool satisfies x strictly-contains z strictly-contains y.

Outputs {out_prefix}_depth{N}_ltr.tsv and {out_prefix}_depth{N}_ltr.fa
for every observed depth N (shadows the raw per-round files, which are
left untouched).

Before containment is resolved, a host that a tenant in one of its LTRs split
into a truncated call and a free-looking neighbour is rebuilt around that tenant
(ltrquest.rejoin, needs --genome); {out_prefix}_rejoined.tsv lists each one.

Usage:
  python -m ltrquest.reconcile \
      --out-prefix mafft_update \
      --tsv mafft_update_r1_ltr.tsv mafft_update_r2_ltr.tsv ... \
      --fa  mafft_update_r1_ltr.fa  mafft_update_r2_ltr.fa  ... \
      --genome genome.fa --tools-dir tools
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

from . import overlaps as ov
from .annotate import REJOIN_SIDECAR_SUFFIX, write_rejoin_sidecar
from .detect import (
    _is_contained,
    _ltrs_shared,
    iter_fasta,
    ltr_bounds_from_table,
)
from .table import Columns, parse_header

COORD_RE = re.compile(r"^([^:]+):(\d+)-(\d+)")

# Depth -> IUPAC mask char. Must match ltrquest's IUPAC_SEQ so a
# cross-round N from a round-1 element (already embedded in the extracted outer
# sequence) lines up with depth0's char here, depth1 -> R, depth2 -> D, etc.
# V is reserved for the wrapper's far-character and is NOT in this list.
IUPAC_DEPTH_SEQ = ("N", "R", "D", "Y", "S", "W", "K", "M", "B", "H")

# The columns resolve_conflicts ranks calls on (ltrquest.overlaps). None is
# fatal to reading a table, but all of them defaulting is.
ROUND_TSV_FIELDS = ("tsd", "tsd_offset", "motif", "bitscore")


def _check_round_columns(path: str, cols: Columns) -> None:
    """Confirm a round table's header names the columns conflicts are ranked on,
    once per file and before any row of it is used.

    None present means there is no header to key off of at all -- every record
    would rank the same, which is not a failure the caller can see: the survivor
    of each conflict is then whichever record sorted first. That file is refused
    outright. Some present is still a usable table; it is read, but the defaulted
    fields are named so they are on the record.
    """
    present = [name for name in ROUND_TSV_FIELDS if name in cols]
    if not present:
        raise ValueError(
            f"{path}: no header, or a header naming none of "
            f"{', '.join(ROUND_TSV_FIELDS)}; this is not an LTRquest "
            f"element table"
        )
    missing = [name for name in ROUND_TSV_FIELDS if name not in cols]
    if missing:
        print(f"[reconcile] WARNING: {path}: header has no "
              f"{', '.join(missing)}; missing column(s) count as no evidence "
              f"for every row", file=sys.stderr)


def parse_tsv(path: str, round_idx: int) -> Tuple[Optional[str], List[dict]]:
    """Return (header_line_or_None, list_of_record_dicts).

    `rank` is the call's structural evidence (ltrquest.overlaps.evidence), read by
    name rather than position, and a file that names none of its columns is
    refused rather than ranked on defaults.
    """
    header: Optional[str] = None
    cols_spec = Columns.of([])
    checked = False
    recs: List[dict] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if header is None:
                    header = line
                    cols_spec = Columns.from_line(line)
                    _check_round_columns(path, cols_spec)
                    checked = True
                continue
            if not checked:
                _check_round_columns(path, cols_spec)
                checked = True
            cols = line.split("\t")
            if not cols:
                continue
            m = COORD_RE.match(cols[0])
            if not m:
                continue
            chrom, s, e = m.group(1), int(m.group(2)), int(m.group(3))
            key = f"{chrom}:{s}-{e}"
            rank = ov.evidence(cols_spec.get(cols, "tsd"), cols_spec.get(cols, "tsd_offset"),
                               cols_spec.get(cols, "motif"), cols_spec.get(cols, "bitscore"),
                               ov.n_domains(cols_spec.get(cols, "domains")))
            recs.append({
                "key": key,
                "chrom": chrom,
                "s": s,
                "e": e,
                "round": round_idx,
                "col1": cols[0],
                "line": line,
                "cols": cols,
                "rank": rank,
            })
    return header, recs


def build_all_in(survivors: List[dict],
                 ltr_bounds: Dict[str, Tuple[int, int, int, int]],
                 tenant_in_ltr: Optional[Dict[str, Tuple[int, int]]] = None
                 ) -> Dict[str, List[str]]:
    """For each survivor, list all survivors strictly contained by it
    (distinct LTRs required). Returns dict key -> list of descendant keys.

    `tenant_in_ltr` maps a host rebuilt around a tenant in one of its LTRs to that
    tenant's span: the tenant and whatever it holds sit inside the host's LTR
    interval without sharing that LTR, so they are nested all the same.
    """
    tenant_in_ltr = tenant_in_ltr or {}
    all_in: Dict[str, List[str]] = defaultdict(list)
    by_chr: Dict[str, List[dict]] = defaultdict(list)
    for e in survivors:
        by_chr[e["chrom"]].append(e)

    for _chrom, arr in by_chr.items():
        # Sort by start asc, then end desc (outer-first on ties)
        arr.sort(key=lambda x: (x["s"], -x["e"]))
        n = len(arr)
        for i in range(n):
            x = arr[i]
            xe = x["e"]
            for j in range(i + 1, n):
                y = arr[j]
                if y["s"] > xe:
                    break  # no further overlap possible
                if y["e"] > xe:
                    continue  # partial overlap, not containment
                if _is_contained((y["s"], y["e"]), (x["s"], x["e"])) != "a_in_b":
                    continue
                span = tenant_in_ltr.get(x["key"])
                in_tenant = span is not None and span[0] <= y["s"] and y["e"] <= span[1]
                if not in_tenant and _ltrs_shared(x["key"], y["key"], ltr_bounds):
                    continue
                all_in[x["key"]].append(y["key"])
    return all_in


def build_direct_children(all_in: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """y is a direct child of x iff y in all_in[x] and no z in all_in[x]
    has y in all_in[z]."""
    children: Dict[str, List[str]] = defaultdict(list)
    for xk, descendants in all_in.items():
        for y in descendants:
            is_direct = True
            for z in descendants:
                if z == y:
                    continue
                if y in all_in.get(z, ()):
                    is_direct = False
                    break
            if is_direct:
                children[xk].append(y)
    return children


def compute_chain_inward(keys: List[str],
                         children: Dict[str, List[str]]) -> Dict[str, int]:
    """Longest chain from each node going inward (downward)."""
    memo: Dict[str, int] = {}

    def walk(k: str) -> int:
        if k in memo:
            return memo[k]
        kids = children.get(k, [])
        if not kids:
            memo[k] = 0
            return 0
        m = 1 + max(walk(c) for c in kids)
        memo[k] = m
        return m

    for k in keys:
        walk(k)
    return memo


def build_updated_nest_status(key: str,
                              all_in: Dict[str, List[str]]) -> str:
    """Pairwise nest_status string: all strict containments for `key`
    across the entire pool (both inners and outers)."""
    rels: List[Tuple[str, str]] = []
    for y in all_in.get(key, []):
        rels.append(("nest-outer", y))
    for xk, descendants in all_in.items():
        if key in descendants:
            rels.append(("nest-inner", xk))
    if not rels:
        return "."
    seen = set()
    uniq: List[Tuple[str, str]] = []
    for r in rels:
        if r not in seen:
            seen.add(r)
            uniq.append(r)
    return ";".join(f"{role}:{k}" for role, k in uniq)


_ORIENTATION_UNREADABLE_WARNED = False


def _outer_is_revcomped(outer_rec: dict, table_cols: Optional[Columns]) -> bool:
    """True if `outer_rec` names a row whose `orientation` column reads `-`.

    `restate_orientation_to_match_library` keeps this column in step with
    `bounded_fasta_oriented`'s flips, so `-` is a direct readout of what the
    library actually stored rather than something re-derived here that could
    drift out of sync with it.

    A row with no usable value (no header, no `orientation` column, or a row
    too short to hold it) is reported once and treated as forward: guessing
    `-` for a row that cannot say so risks mirroring content the library
    never flipped, which is worse than leaving an already-rare case unmarked.
    """
    global _ORIENTATION_UNREADABLE_WARNED
    row = outer_rec.get("cols")
    if table_cols is not None and row is not None and "orientation" in table_cols:
        orient = table_cols.get(row, "orientation")
        if orient in ("+", "-"):
            return orient == "-"
    if not _ORIENTATION_UNREADABLE_WARNED:
        print("[reconcile] WARNING: no usable 'orientation' column for a "
              "nest-outer; its depth mask will not be mirrored even if the "
              "library stores it reverse-complemented", file=sys.stderr)
        _ORIENTATION_UNREADABLE_WARNED = True
    return False


def apply_depth_masking(outer_seq: str,
                        outer_rec: dict,
                        direct_children: Dict[str, List[str]],
                        rec_by_key: Dict[str, dict],
                        depth_map: Dict[str, int],
                        table_cols: Optional[Columns] = None) -> str:
    """Return outer_seq with every descendant region overwritten by the IUPAC
    char corresponding to that descendant's depth.

    The traversal is parent-first / grandchild-after so deeper descendants
    overwrite the shallower marks placed by their ancestors, producing the
    depth-indexed pattern:
        depth-1 outer  -> ...N... (direct child at depth 0)
        depth-2 outer  -> ...R..N..R... (direct child at depth 1 with its own
                                         depth-0 grandchild inside)
    """
    outer_s = outer_rec["s"]
    outer_chrom = outer_rec["chrom"]
    outer_len = len(outer_seq)
    chars = list(outer_seq)
    is_revcomped = _outer_is_revcomped(outer_rec, table_cols)

    def paint(parent_key: str) -> None:
        for child_key in direct_children.get(parent_key, []):
            if child_key == parent_key:
                continue
            child = rec_by_key.get(child_key)
            if child is None or child["chrom"] != outer_chrom:
                continue
            cd = depth_map.get(child_key, 0)
            if 0 <= cd < len(IUPAC_DEPTH_SEQ):
                ch = IUPAC_DEPTH_SEQ[cd]
            else:
                ch = "X"
            # 1-based inclusive coords -> 0-based half-open relative to outer
            rel_s = max(0, child["s"] - outer_s)
            rel_e = min(outer_len, child["e"] - outer_s + 1)
            if rel_e <= rel_s:
                continue
            if is_revcomped:
                # outer_seq is the library's own record: minus-strand elements
                # are stored reverse-complemented (bounded_fasta_oriented),
                # while child coords above are always forward-genomic, so the
                # interval has to be mirrored within the record before it is
                # painted -- see mask_same_round_inners_in_fa, which mirrors
                # for the same reason on the round-local twin of this mask.
                rel_s, rel_e = outer_len - rel_e, outer_len - rel_s
            for p in range(rel_s, rel_e):
                chars[p] = ch
            paint(child_key)

    paint(outer_rec["key"])
    return "".join(chars)


def resolve_conflicts(pool: List[dict],
                      ltr_bounds: Dict[str, Tuple[int, int, int, int]]
                      ) -> Tuple[List[dict], int]:
    """Pooled over the rounds, any two calls that overlap without one being properly
    nested in the other keep only the stronger (ltrquest.overlaps), as detection does
    within a round. Returns (survivors in pool order, n dropped).
    """
    calls = [ov.Call(r["key"], r["chrom"], r["s"], r["e"], ltr_bounds.get(r["key"]), r["rank"])
             for r in pool]
    kept, dropped = ov.resolve(calls)
    kept_ids = {id(c) for c in kept}
    return [r for r, c in zip(pool, calls) if id(c) in kept_ids], len(dropped)


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        description="Reconcile nested LTR-RT calls across rounds into "
                    "depth-bucketed libraries."
    )
    ap.add_argument("--tsv", nargs="+", required=True,
                    help="Per-round _ltr.tsv files (in round order).")
    ap.add_argument("--fa", nargs="+", required=True,
                    help="Per-round _ltr.fa files matching --tsv.")
    ap.add_argument("--out-prefix", required=True,
                    help="Output prefix for _depth{N}_ltr.{tsv,fa} files.")
    ap.add_argument("--genome", default=None,
                    help="The genome the rounds were run on (.fa or .fa.gz). Needed to "
                         "rebuild a host a tenant in one of its LTRs split; without it "
                         "such pairs stay two adjacent calls.")
    ap.add_argument("--tools-dir", default="./tools",
                    help="Where a Kmer2LTR checkout lives if it is not installed "
                         "(its TSD rule is read at rebuilt ends). Default: ./tools")
    args = ap.parse_args(argv)

    if len(args.tsv) != len(args.fa):
        ap.error("--tsv and --fa must have the same number of entries")

    # 1. Load pool
    pool: List[dict] = []
    header: Optional[str] = None
    for r_idx, tsv in enumerate(args.tsv, start=1):
        h, recs = parse_tsv(tsv, r_idx)
        if header is None and h is not None:
            header = h
        pool.extend(recs)
    print(f"[reconcile] loaded {len(pool)} records from {len(args.tsv)} rounds",
          file=sys.stderr)

    # Step 6 below rewrites the last column in place, on the assumption that
    # it is nest_status; a schema that moved nest_status elsewhere would make
    # that a silent corruption instead of a failure here.
    table_cols: Optional[Columns] = None
    if header is not None:
        names = parse_header(header)
        if not names or names[-1] != "nest_status":
            raise ValueError(
                f"expected the element table header to end in 'nest_status', "
                f"got: {header!r}"
            )
        table_cols = Columns.of(names)

    # 2. Load LTR boundaries (union across rounds). Read from each round's own
    # table rather than a separate file, so a boundary and the key it is
    # filed under always come from the same row -- see ltr_bounds_from_table.
    ltr_bounds: Dict[str, Tuple[int, int, int, int]] = {}
    for tsv in args.tsv:
        lb = ltr_bounds_from_table(tsv)
        for k, v in lb.items():
            ltr_bounds.setdefault(k, v)
    print(f"[reconcile] LTR boundaries loaded: {len(ltr_bounds)} keys",
          file=sys.stderr)
    if len(ltr_bounds) < len(pool):
        print(f"[reconcile] WARNING: LTR boundary map covers {len(ltr_bounds)} "
              f"of {len(pool)} pooled element(s); {len(pool) - len(ltr_bounds)} "
              f"have no usable LTR coordinates and will be treated as sharing "
              f"no LTRs with anything", file=sys.stderr)

    # 3. Calls from different rounds that overlap without proper nesting
    survivors, n_dropped = resolve_conflicts(pool, ltr_bounds)
    print(f"[reconcile] dropped {n_dropped} call(s) overlapping a stronger call from "
          f"another round without being nested in it", file=sys.stderr)

    # Load FASTA records keyed by full header (col1 == "chrom:s-e#class"); a
    # rebuilt host's record is built from its stored one in step 3b.
    fa_seqs: Dict[str, str] = {}
    for fa in args.fa:
        for h, seq in iter_fasta(fa):
            if h not in fa_seqs or len(seq) > len(fa_seqs[h]):
                fa_seqs[h] = seq

    # 3b. Rebuild hosts that a tenant in one of their LTRs split in two.
    # Always written, a header alone when nothing was rebuilt, so the strand
    # loaders never read one left by an earlier run.
    rejoins = []
    sidecar_rows: List[List[str]] = []
    if args.genome and table_cols is not None:
        # Imported here: it aligns (parasail, pyfaidx), and reconciling without
        # --genome needs nothing beyond the standard library.
        from . import kmer2ltr, rejoin
        rejoins = rejoin.run(survivors, ltr_bounds, fa_seqs, table_cols, args.genome,
                             kmer2ltr.genome_module(args.tools_dir),
                             os.path.dirname(os.path.abspath(args.out_prefix)),
                             lambda msg: print(f"[reconcile] {msg}", file=sys.stderr))
        sidecar_rows = rejoin.sidecar_rows(rejoins, table_cols)
    elif survivors:
        print("[reconcile] WARNING: no --genome; a host split by a tenant in one of "
              "its LTRs stays two adjacent calls", file=sys.stderr)
    write_rejoin_sidecar(args.out_prefix + REJOIN_SIDECAR_SUFFIX, sidecar_rows)

    # 4. Containment graph
    tenant_in_ltr = {rj.key: (rj.tenant["s"], rj.tenant["e"]) for rj in rejoins}
    all_in = build_all_in(survivors, ltr_bounds, tenant_in_ltr)
    children = build_direct_children(all_in)

    # 5. Inward-chain depth
    keys = [e["key"] for e in survivors]
    depth = compute_chain_inward(keys, children)

    rec_by_key: Dict[str, dict] = {e["key"]: e for e in survivors}

    # 6. Bucket survivors; rewrite nest_status col with cross-round view
    buckets: Dict[int, List[dict]] = defaultdict(list)
    for e in survivors:
        d = depth.get(e["key"], 0)
        new_cols = list(e["cols"])
        new_cols[-1] = build_updated_nest_status(e["key"], all_in)
        e["line_out"] = "\t".join(new_cols)
        buckets[d].append(e)

    # 7. Write outputs
    out_prefix = args.out_prefix
    any_written = False
    n_repainted = 0
    for d in sorted(buckets):
        recs = buckets[d]
        tsv_out = f"{out_prefix}_depth{d}_ltr.tsv"
        fa_out = f"{out_prefix}_depth{d}_ltr.fa"
        with open(tsv_out, "w") as tout:
            if header:
                tout.write(header + "\n")
            for r in recs:
                tout.write(r["line_out"] + "\n")
        with open(fa_out, "w") as fout:
            for r in recs:
                seq = fa_seqs.get(r["col1"])
                if seq is None:
                    continue
                if d > 0:
                    seq = apply_depth_masking(seq, r, children, rec_by_key, depth,
                                              table_cols)
                    n_repainted += 1
                fout.write(f">{r['col1']}\n")
                for i in range(0, len(seq), 60):
                    fout.write(seq[i:i + 60] + "\n")
        any_written = True
        print(f"[reconcile] depth{d}: {len(recs)} -> {tsv_out}, {fa_out}",
              file=sys.stderr)
    if n_repainted:
        print(f"[reconcile] depth-indexed IUPAC masking applied to {n_repainted} "
              f"outer sequence(s)", file=sys.stderr)
    if not any_written:
        print("[reconcile] no records to write", file=sys.stderr)

    # Summary line for logs
    total = sum(len(v) for v in buckets.values())
    parts = [f"depth{d}={len(buckets[d])}" for d in sorted(buckets)]
    print(f"[reconcile] total={total} ({', '.join(parts)})", file=sys.stderr)


if __name__ == "__main__":
    main()
