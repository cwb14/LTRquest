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

Usage:
  python -m ltrquest.reconcile \
      --out-prefix mafft_update \
      --tsv mafft_update_r1_ltr.tsv mafft_update_r2_ltr.tsv ... \
      --fa  mafft_update_r1_ltr.fa  mafft_update_r2_ltr.fa  ...
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .detect import (
    _is_contained,
    _ltrs_shared,
    _sub_dedup_shared_ltr_group,
    iter_fasta,
    ltr_bounds_from_table,
)
from .table import Columns, as_float, as_int, parse_header

COORD_RE = re.compile(r"^([^:]+):(\d+)-(\d+)")

# Depth -> IUPAC mask char. Must match ltrquest's IUPAC_SEQ so a
# cross-round N from a round-1 element (already embedded in the extracted outer
# sequence) lines up with depth0's char here, depth1 -> R, depth2 -> D, etc.
# V is reserved for the wrapper's far-character and is NOT in this list.
IUPAC_DEPTH_SEQ = ("N", "R", "D", "Y", "S", "W", "K", "M", "B", "H")

# The columns cross_round_dedup scores survivors on. Neither is fatal to
# reading a table, but both defaulting is.
ROUND_TSV_FIELDS = ("p_dist", "aln_len")


def _check_round_columns(path: str, cols: Columns) -> None:
    """Confirm a round table's header names the columns dedup scores on, once
    per file and before any row of it is used.

    Neither present means there is no header to key off of at all -- every
    record would come back with p_dist 0 and aln_len 0, which is not a failure
    the caller can see: the score `aln * (1 - p)` is then uniformly zero and
    the survivor of each duplicate group is whichever record sorted first.
    That file is refused outright. One of the two present is still a usable
    table; it is read, but the defaulted field is named so it is on the record.
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
              f"{', '.join(missing)}; missing column(s) default to zero for "
              f"every row", file=sys.stderr)


def parse_tsv(path: str, round_idx: int) -> Tuple[Optional[str], List[dict]]:
    """Return (header_line_or_None, list_of_record_dicts).

    `p` (p_dist) and `aln` (aln_len) feed _sub_dedup_shared_ltr_group's
    matching-bases score, so they are read by name rather than position, and a
    file that names neither is refused rather than scored on defaults.
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
            p = as_float(cols_spec.get(cols, "p_dist"), default=0.0)
            aln = as_int(cols_spec.get(cols, "aln_len"), default=0)
            recs.append({
                "key": key,
                "chrom": chrom,
                "s": s,
                "e": e,
                "round": round_idx,
                "col1": cols[0],
                "line": line,
                "cols": cols,
                "p": p,
                "aln": aln,
            })
    return header, recs


def build_all_in(survivors: List[dict],
                 ltr_bounds: Dict[str, Tuple[int, int, int, int]]
                 ) -> Dict[str, List[str]]:
    """For each survivor, list all survivors strictly contained by it
    (distinct LTRs required). Returns dict key -> list of descendant keys.
    """
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
                if _ltrs_shared(x["key"], y["key"], ltr_bounds):
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


def cross_round_dedup(pool: List[dict],
                      ltr_bounds: Dict[str, Tuple[int, int, int, int]]
                      ) -> Tuple[List[dict], int]:
    """Union-find across the pool: any two records (from different rounds)
    that share LTR boundaries are merged; keep the best per group via
    _sub_dedup_shared_ltr_group. Returns (survivors, n_merged_groups).
    """
    n = len(pool)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    by_chr: Dict[str, List[int]] = defaultdict(list)
    for i, e in enumerate(pool):
        by_chr[e["chrom"]].append(i)

    for _chrom, idxs in by_chr.items():
        m = len(idxs)
        for ii in range(m):
            i = idxs[ii]
            for jj in range(ii + 1, m):
                j = idxs[jj]
                if pool[i]["round"] == pool[j]["round"]:
                    continue
                if _ltrs_shared(pool[i]["key"], pool[j]["key"], ltr_bounds):
                    union(i, j)

    groups: Dict[int, List[int]] = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)

    survivors: List[dict] = []
    n_merged = 0
    for members in groups.values():
        if len(members) == 1:
            survivors.append(pool[members[0]])
        else:
            n_merged += 1
            recs = [pool[m] for m in members]
            best = _sub_dedup_shared_ltr_group(recs)
            survivors.append(best)
    return survivors, n_merged


def main() -> None:
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
    args = ap.parse_args()

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

    # 3. Cross-round dedup (shared-LTR collapse)
    survivors, n_merged = cross_round_dedup(pool, ltr_bounds)
    if n_merged:
        print(f"[reconcile] cross-round shared-LTR merges: {n_merged} groups",
              file=sys.stderr)
    else:
        print("[reconcile] no cross-round shared-LTR merges", file=sys.stderr)

    # 4. Containment graph
    all_in = build_all_in(survivors, ltr_bounds)
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

    # 7. Load FASTA records keyed by full header (col1 == "chrom:s-e#class")
    fa_seqs: Dict[str, str] = {}
    for fa in args.fa:
        for h, seq in iter_fasta(fa):
            if h not in fa_seqs or len(seq) > len(fa_seqs[h]):
                fa_seqs[h] = seq

    # 8. Write outputs
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
