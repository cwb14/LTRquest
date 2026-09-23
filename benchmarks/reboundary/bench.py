#!/usr/bin/env python3
"""Benchmarks for ltrquest.reboundary (design spec section 9).

  b1       planted obstacles with known true ends (B1) plus untouched controls (B3)
  b2       a finished run re-bounded in a scratch copy (B2, and B3 on TSD-bearing calls)
  collect  one table from every <out>/*/summary.json

  python bench.py b1 --run RUN --prefix P [P ...] --genome G [G ...] --tools-dir T --out DIR
                     [--n 1000] [--seed 11] [--threads 32] [--set key=value ...]
  python bench.py b2 --run RUN --prefix P [P ...] --genome G [G ...] --tools-dir T --out DIR
                     [--threads 32] [--cache FILE] [--set key=value ...]
  python bench.py collect --out DIR

--set takes any reboundary.Settings field (method, references, min_copies, credit,
max_ratio, qc_min_n, untested, subfamily_jaccard) or ltr_place.Params field
(min_identity, min_whole_identity, min_ext, anchor, anchor_len, max_indel, n_models).
The run directory is only read: b2 works on copies under <out>/run/.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import random
import resource
import shutil
import time
import zlib
from collections import Counter, defaultdict
from dataclasses import asdict, fields, replace
from multiprocessing import Pool
from typing import Dict, List, Optional, Sequence

from ltrquest import reboundary as rb
from ltrquest.ltr_model import Genomes, Member, ratio_ok, tsd_enrichment
from ltrquest.ltr_place import Params, propose
from ltrquest.reboundary_io import load_clean_tables, members_from

OBSTACLES = [("none", 0), ("patch", 30), ("del", 10), ("del", 50), ("del", 200),
             ("ins", 50), ("ins", 300), ("ins", 1000), ("ins", 5000),
             ("inste", 300), ("inste", 1000)]     # inste: a real segment of another family
DELTAS = (20, 50, 100, 200, 400)
FLANK = 3000
BENCH = "bench"
NOT_FOUND = ("", ".", "NA")


def settings(pairs: Sequence[str], threads: int, mafft: str) -> rb.Settings:
    s = rb.Settings(threads=threads, mafft=mafft)
    place_names = {f.name for f in fields(Params)}
    top_names = {f.name for f in fields(rb.Settings)} - {"place"}
    place = {}
    for kv in pairs:
        key, value = kv.split("=", 1)
        if key in place_names:
            cur = getattr(Params(), key)
            place[key] = ((value.lower() in ("1", "true", "yes")) if isinstance(cur, bool)
                          else type(cur)(value))
        elif key in top_names:
            s = replace(s, **{key: type(getattr(s, key))(value)})
        else:
            raise SystemExit(f"bench: unknown setting {key!r}")
    return replace(s, place=replace(s.place, **place))


def recorded(s: rb.Settings) -> Dict:
    """The settings as they actually ran: `place_params` overrides `combine` for `nearest`.

    `asdict(s)` records the placement parameters before that override, so a
    `nearest` run used to be published as `combine: "best"` when every
    placement was combined by median.
    """
    return asdict(replace(s, place=rb.place_params(s)))


def frac(a: int, n: int) -> Optional[float]:
    return round(a / n, 4) if n else None


def md(header: Sequence[str], rows: Sequence[Sequence]) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join("" if v is None else str(v) for v in r) + " |" for r in rows]
    return "\n".join(out)


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path) as fh:
        head = fh.readline().rstrip("\n").lstrip("#").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def write_tsv(path: str, rows: Sequence[Dict]) -> None:
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def found(tsd: str) -> bool:
    return tsd not in NOT_FOUND


# ---------------------------------------------------------------- B1 / B3 (planted)
def pick_truth(members: List[Member], s: rb.Settings, n: int, seed: int) -> List[Member]:
    """TSD-bearing calls (independent evidence the ends are right), <= 20 per family."""
    fams: Dict[str, List[Member]] = defaultdict(list)
    for m in members:
        fams[m.family].append(m)
    pool: List[Member] = []
    for fam, ms in sorted(fams.items()):
        if fam in NOT_FOUND or len(ms) < s.min_copies:
            continue
        ok = sorted((m for m in ms if m.has_tsd and m.stranded
                     and min(m.len_left, m.len_right) >= 150), key=lambda m: m.uid)
        random.Random(zlib.crc32(fam.encode()) ^ seed).shuffle(ok)
        pool += ok[:20]
    random.Random(seed).shuffle(pool)
    return pool[:n]


def scenarios(n_truth: int, seed: int):
    rng = random.Random(seed)
    for i in range(n_truth):
        for kind, size in OBSTACLES:
            yield (i, kind, size, rng.choice(DELTAS), rng.choice(("left", "right")),
                   rng.choice(("same", "other")))
        yield i, "control", 0, 0, "left", "same"


def plant(g: Genomes, t: Member, kind: str, size: int, delta: int, end: str, copy: str,
          rng: random.Random, cid: str, donors: Sequence[Member] = ()):
    """One synthetic contig: the truth element plus 3 kb flanks, an obstacle planted in one LTR
    `delta` bp from an outer end, and the call cut at the obstacle the way a finder stops.
    `copy`: 'same' = the LTR whose outer end is truncated, 'other' = its partner.
    `inste` inserts `size` bp of another family's element (a TE-like insertion) instead of
    random sequence. Returns (contig, called Member, (true left, true right)) or None."""
    geo = "ins" if kind == "inste" else kind
    lo = max(1, t.start - FLANK)
    hi = min(g.length(t.prefix, t.chrom), t.end + FLANK)
    seq, _ = g.fetch(t.prefix, t.chrom, lo, hi)
    L0, L1, R0, R1 = (x - lo + 1 for x in (t.start, t.l1, t.r0, t.end))
    ltr = min(L1 - L0 + 1, R1 - R0 + 1)
    width = {"none": 0, "patch": 30, "del": size, "ins": 0, "control": 0}[geo]
    if kind != "control" and delta + width + 30 >= ltr // 2:
        return None
    olen = size if geo in ("del", "ins") else 0
    p = None
    if kind != "control":
        if end == "left":
            p = (L0 if copy == "same" else R0) + delta
        else:
            last = R1 if copy == "same" else L1
            p = last - delta - width + 1 if geo in ("del", "patch") else last - delta + 1
    if geo == "del":
        seq = seq[:p - 1] + seq[p - 1 + olen:]
    elif kind == "ins":
        seq = seq[:p - 1] + "".join(rng.choice("ACGT") for _ in range(olen)) + seq[p - 1:]
    elif kind == "inste":
        pool = [d for d in donors if d.family != t.family and d.end - d.start + 1 >= olen + 200]
        if not pool:
            return None
        d = rng.choice(pool)
        seg, _ = g.fetch(d.prefix, d.chrom, d.start + 100, d.start + 99 + olen)
        if len(seg) != olen:
            return None
        seq = seq[:p - 1] + seg + seq[p - 1:]
    elif kind == "patch":
        chars = list(seq)
        for i in rng.sample(range(p - 1, p + 29), 9):
            chars[i] = rng.choice([b for b in "ACGT" if b != chars[i]])
        seq = "".join(chars)

    def moved(x: int) -> int:
        if p is None or x < p:
            return x
        return x - olen if geo == "del" else (x + olen if geo == "ins" else x)

    tL0, tL1, tR0, tR1 = (moved(x) for x in (L0, L1, R0, R1))
    if kind == "control":
        start, stop, l1, r0 = tL0, tR1, tL1, tR0
    else:
        same = copy == "same"
        off = {"none": delta, "patch": delta + 30,
               "del": delta if same else delta + olen,
               "ins": delta + olen if same else delta}[geo]
        if end == "left":
            start, stop, l1 = tL0 + off, tR1, tL1
            r0 = tR1 - (l1 - start)
        else:
            start, stop, r0 = tL0, tR1 - off, tR0
            l1 = tL0 + (stop - r0)
    suffix = "#" + t.name.split("#", 1)[1] if "#" in t.name else ""
    m = Member(prefix=BENCH, name=f"{cid}:{start}-{stop}{suffix}", chrom=cid, start=start,
               end=stop, l1=l1, r0=r0, strand=t.strand, orientation="+", family=t.family,
               depth=0, k2p=t.k2p, tsd=".", nest_status=".")
    return seq, m, (tL0, tR1)


def outcome(m: Member, truth, left: int, right: int, kind: str) -> str:
    tl, tr = truth
    if kind == "control":
        return "kept" if (left, right) == (tl, tr) else "false_change"
    if (left, right) == (m.start, m.end):
        return "unchanged"
    if left < tl - 5 or right > tr + 5:
        return "over"
    d = max(abs(left - tl), abs(right - tr))
    return "exact" if d == 0 else "within1" if d <= 1 else "within5" if d <= 5 else "wrong"


def _b1_family(job):
    family, real, synth, s, exclude = job
    models, _modal, status = rb.build_models(family, real, rb._G, s, exclude)
    kept = [mo for mo in models if ratio_ok(mo, s.max_ratio)]
    if models and not kept:
        status = "ratio_failed"
    place = rb.place_params(s)
    return [(m.uid, propose(m, kept, rb._G, place) if kept else None, status) for m in synth]


def summarize_b1(rows: List[Dict]) -> Dict:
    groups = defaultdict(list)
    for r in rows:
        groups[(r["kind"], r["size"])].append(r["outcome"])
    table = []
    for (kind, size), outs in sorted(groups.items()):
        c, n = Counter(outs), len(outs)
        table.append([kind, size, n] + [frac(c[k], n) for k in
                                        ("exact", "within1", "within5", "over", "wrong",
                                         "unchanged", "false_change")])
    obst = Counter(r["outcome"] for r in rows if r["kind"] != "control")
    n_obst = sum(obst.values())
    ctrl = Counter(r["outcome"] for r in rows if r["kind"] == "control")
    return {"exact_or_1": frac(obst["exact"] + obst["within1"], n_obst),
            "exact": frac(obst["exact"], n_obst), "over": frac(obst["over"], n_obst),
            "wrong": frac(obst["wrong"], n_obst), "unchanged": frac(obst["unchanged"], n_obst),
            "false_change": frac(ctrl["false_change"], sum(ctrl.values())), "table": table}


def b1(args) -> None:
    s = settings(args.set, args.threads, args.mafft)
    os.makedirs(args.out, exist_ok=True)
    paths = dict(zip(args.prefix, args.genome))
    members = [m for p in args.prefix for m in members_from(p, load_clean_tables(args.run, p))]
    truth = pick_truth(members, s, args.n, args.seed)
    g = Genomes(paths)
    rng = random.Random(args.seed)
    donors = random.Random(args.seed + 1).sample(
        [m for m in members if m.end - m.start + 1 >= 6000],
        min(2000, sum(1 for m in members if m.end - m.start + 1 >= 6000)))
    synth, contigs = [], []
    for j, (i, kind, size, delta, end, copy) in enumerate(scenarios(len(truth), args.seed)):
        got = plant(g, truth[i], kind, size, delta, end, copy, rng, f"b{j}", donors)
        if got is not None:
            seq, m, tr = got
            contigs.append((m.chrom, seq))
            synth.append((m, tr, (kind, size, delta, end, copy), truth[i]))
    fasta = os.path.join(args.out, "synthetic.fa")
    with open(fasta, "w") as fh:
        for name, seq in contigs:
            fh.write(f">{name}\n")
            for k in range(0, len(seq), 80):
                fh.write(seq[k:k + 80] + "\n")
    if os.path.exists(fasta + ".fai"):
        os.remove(fasta + ".fai")
    paths[BENCH] = fasta
    gb = Genomes({BENCH: fasta})
    gb.length(BENCH, contigs[0][0])        # build the .fai here, not in racing workers
    real = defaultdict(list)
    for m in members:
        real[m.family].append(m)
    exclude = defaultdict(set)
    for t in truth:
        exclude[t.family].add(t.key)
    by_family = defaultdict(list)
    for m, _, _, _ in synth:
        by_family[m.family].append(m)
    jobs = [(f, real[f], ms, s, frozenset(exclude[f])) for f, ms in sorted(by_family.items())]
    t0 = time.time()
    with Pool(s.threads, initializer=rb._init, initargs=(paths, args.tools_dir)) as pool:
        placed = {}
        for chunk in pool.imap_unordered(_b1_family, jobs):
            for uid, prop, status in chunk:
                placed[uid] = (prop, status)
        arb = []
        for m, _, _, _ in synth:
            prop, _ = placed.get(m.uid, (None, "no_model"))
            if prop is not None and prop.gate_ok:
                rec = gb.fetch(BENCH, m.chrom, m.start, m.end)[0]
                arb.append((prop, rec, rb.credit_for(prop, s), s.mutation_rate, s.place.min_ext))
        verdicts = {v.uid: v for v in pool.map(rb.arbitrate, arb, chunksize=8)}
    rows = []
    for m, tr, (kind, size, delta, end, copy), t in synth:
        prop, status = placed.get(m.uid, (None, "no_model"))
        v = verdicts.get(m.uid)
        left, right = (v.accepted.left, v.accepted.right) if v and v.accepted else (m.start, m.end)
        rows.append({"kind": kind, "size": size, "delta": delta, "end": end, "copy": copy,
                     "family": t.family, "model": status,
                     "proposal": "none" if prop is None else ("ok" if prop.gate_ok else "gated"),
                     "verdict": v.status if v else ".", "truth_left": tr[0], "truth_right": tr[1],
                     "called_left": m.start, "called_right": m.end, "final_left": left,
                     "final_right": right, "outcome": outcome(m, tr, left, right, kind)})
    write_tsv(os.path.join(args.out, "b1.tsv"), rows)
    summary = summarize_b1(rows)
    summary.update(settings=recorded(s), n_truth=len(truth), n_contigs=len(synth),
                   wall_s=round(time.time() - t0))
    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(md(["kind", "size", "n", "exact", "±1", "±5", "over", "wrong", "unchanged",
              "false change"], summary["table"]))
    print(f"\nB1 exact-or-±1 {summary['exact_or_1']}  over {summary['over']}  "
          f"wrong {summary['wrong']}  | B3 false change {summary['false_change']}  "
          f"({len(truth)} truth elements, {summary['wall_s']} s)")


# ---------------------------------------------------------------- B2 / B3 (real run)
def cross_validated(props: List[Dict[str, str]], s: rb.Settings) -> Dict:
    """Family QC chosen on one half of each family's candidates, TSD gain measured on the other."""
    halves = defaultdict(lambda: ([], []))
    for r in props:
        if (r["gate_ok"] != "1" or r["called_tsd"] == "1"
                or r["family_status"] in ("ratio_failed", "too_few_references", "no_ltr_span")):
            continue
        halves[r["family"]][zlib.crc32(r["uid"].encode()) % 2].append(r)
    train = [x for a, _ in halves.values() for x in a]
    p0 = (sum(found(x["tsd_null"]) for x in train) / len(train)) if train else 0.0
    n = hits = null = 0
    for a, b in halves.values():
        st = tsd_enrichment(sum(found(x["tsd_new"]) for x in a), len(a), p0, s.qc_min_n)
        if st == "pass" or (st == "untested" and s.untested == "accept"):
            n += len(b)
            hits += sum(found(x["tsd_new"]) for x in b)
            null += sum(found(x["tsd_null"]) for x in b)
    return {"cv_n": n, "cv_tsd": frac(hits, n), "cv_null": frac(null, n), "cv_p0": round(p0, 4)}


def b2(args) -> None:
    s = settings(args.set, args.threads, args.mafft)
    work = os.path.join(args.out, "run")
    os.makedirs(work, exist_ok=True)
    for p in args.prefix:
        for path in (glob.glob(os.path.join(args.run, f"{p}_depth*_clean_ltr.*"))
                     + glob.glob(os.path.join(args.run, f"{p}.detect.json"))):
            shutil.copy2(path, os.path.join(work, os.path.basename(path)))
    tsd_calls = sum(m.has_tsd for p in args.prefix
                    for m in members_from(p, load_clean_tables(args.run, p)))
    dump = os.path.join(args.out, "proposals.tsv")
    t0 = time.time()
    counts = rb.run(work, args.prefix, args.genome, s, args.tools_dir, dump=dump,
                    cache=args.cache)
    wall = time.time() - t0
    rss_kb = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    side = [r for p in args.prefix for r in read_tsv(os.path.join(work, p + "_reboundary.tsv"))]
    ext = [r for r in side if r["decision"] == "extended"]
    fresh = [r for r in ext if not found(r["tsd_called"])]
    had = [r for r in ext if found(r["tsd_called"])]
    exts = sorted(max(int(r["ext5"]), int(r["ext3"])) for r in ext)
    dk = sorted(float(r["k2p_new"]) - float(r["k2p_called"]) for r in ext
                if r["k2p_new"] not in NOT_FOUND and r["k2p_called"] not in NOT_FOUND)
    q = lambda xs, f: xs[min(len(xs) - 1, int(f * len(xs)))] if xs else None  # noqa: E731
    summary = {"settings": recorded(s), "candidates": len(side), "extended": len(ext),
               "counts": counts, "n_fresh": len(fresh),
               "tsd_gain": frac(sum(found(r["tsd_new"]) for r in fresh), len(fresh)),
               "tsd_null": frac(sum(found(r["tsd_null"]) for r in fresh), len(fresh)),
               "b3_real_changed": frac(len(had), tsd_calls),
               "b3_real_tsd_kept": frac(sum(found(r["tsd_new"]) for r in had), len(had)),
               "ext_median": q(exts, 0.5), "ext_q90": q(exts, 0.9),
               "dk2p_median": q(dk, 0.5), "dk2p_q90": q(dk, 0.9),
               "wall_s": round(wall), "peak_rss_gb": round(rss_kb / 1e6, 2)}
    summary.update(cross_validated(read_tsv(dump), s))
    motif_called = {r[0]: t.cols.get(r, "motif") for p in args.prefix
                    for t in load_clean_tables(args.run, p) for r in t.rows}
    motif_new = {r[0]: t.cols.get(r, "motif") for p in args.prefix
                 for t in load_clean_tables(work, p) for r in t.rows}
    summary["tgca_called"] = frac(sum(motif_called.get(r["old_seq_id"]) == "tg...ca" for r in ext),
                                  len(ext))
    summary["tgca_new"] = frac(sum(motif_new.get(r["new_seq_id"]) == "tg...ca" for r in ext),
                               len(ext))

    def breakdown(label):
        groups = defaultdict(list)
        for r in side:
            groups[label(r)].append(r)
        out = {}
        for k, rs in sorted(groups.items()):
            e = [r for r in rs if r["decision"] == "extended"]
            f = [r for r in e if not found(r["tsd_called"])]
            out[k] = {"candidates": len(rs), "extended": len(e),
                      "tsd_gain": frac(sum(found(r["tsd_new"]) for r in f), len(f))}
        return out

    def age(r):
        k = r["k2p_called"]
        if k in NOT_FOUND:
            return "NA"
        k = float(k)
        return ("<0.005" if k < 0.005 else "0.005-0.02" if k < 0.02
                else "0.02-0.05" if k < 0.05 else ">=0.05")

    summary["by_clade"] = breakdown(lambda r: r["old_seq_id"].rsplit("/", 1)[-1]
                                    if "/" in r["old_seq_id"] else "unknown")
    summary["by_age"] = breakdown(age)
    with open(os.path.join(args.out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(md(["metric", "value"], [[k, v] for k, v in summary.items()
                                   if k not in ("settings", "by_clade", "by_age")]))
    for title in ("by_clade", "by_age"):
        print(f"\n{title}\n")
        print(md(["group", "candidates", "extended", "TSD gain"],
                 [[k, v["candidates"], v["extended"], v["tsd_gain"]]
                  for k, v in summary[title].items()]))


# ---------------------------------------------------------------- collect
def collect(args) -> None:
    rows1, rows2 = [], []
    for path in sorted(glob.glob(os.path.join(args.out, "*", "summary.json"))):
        name = os.path.basename(os.path.dirname(path))
        with open(path) as fh:
            s = json.load(fh)
        if name.startswith("b1_"):
            rows1.append([name[3:], s["exact_or_1"], s["exact"], s["over"], s["wrong"],
                          s["unchanged"], s["false_change"], s["wall_s"]])
        elif name.startswith("b2_"):
            rows2.append([name[3:], s["extended"], s["tsd_gain"], s["tsd_null"], s["cv_tsd"],
                          s["cv_null"], s["tgca_called"], s["tgca_new"], s["b3_real_changed"],
                          s["b3_real_tsd_kept"], s["dk2p_q90"], s["wall_s"], s["peak_rss_gb"]])
    if rows1:
        print("B1 planted obstacles (B3 = false change on untouched controls)\n")
        print(md(["setting", "exact/±1", "exact", "over", "wrong", "unchanged", "B3 false",
                  "s"], rows1))
    if rows2:
        print("\nB2 real run (cv = family QC cross-validated)\n")
        print(md(["setting", "extended", "TSD gain", "null", "cv TSD", "cv null", "TGCA called",
                  "TGCA new", "B3 changed", "B3 TSD kept", "dK2P q90", "s", "RSS GB"], rows2))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=("b1", "b2", "collect"))
    ap.add_argument("--run", help="finished LTRquest run directory (read only)")
    ap.add_argument("--prefix", nargs="+")
    ap.add_argument("--genome", nargs="+")
    ap.add_argument("--tools-dir", help="Kmer2LTR checkout at the pinned commit")
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=1000, help="b1: truth elements (default 1000)")
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--mafft", default="mafft")
    ap.add_argument("--cache", default=None, help="b2: reuse a family phase across sweeps")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    args = ap.parse_args()
    if args.mode == "collect":
        collect(args)
        return
    if not (args.run and args.prefix and args.genome and args.tools_dir):
        raise SystemExit("bench: b1/b2 need --run, --prefix, --genome and --tools-dir")
    if len(args.prefix) != len(args.genome):
        raise SystemExit("bench: --prefix and --genome must pair up one to one")
    (b1 if args.mode == "b1" else b2)(args)


if __name__ == "__main__":
    main()
