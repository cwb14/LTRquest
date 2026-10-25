"""Step 8e: bounded calls moved onto canonical ends, re-measured by Kmer2LTR."""

from __future__ import annotations

import os
import random

import pytest

from ltrquest import detect
from ltrquest.kmer2ltr import COLUMNS

TOOLS = os.environ.get("LTRQUEST_TOOLS_DIR", ".")
RNG = random.Random(11)


def rand(n):
    return "".join(RNG.choice("ACGT") for _ in range(n))


LTR = "TG" + rand(396) + "CA"
INT = rand(1999) + "C"          # ends C, as the base before A's TSD copy: see `wide`
LTR_B = "TG" + rand(396) + "CA"


def build():
    """chr1: A (5-bp TSD AATGC), B (4-bp TSD TAAA), D (no TSD). Returns (seq, spans)."""
    parts, spans = [rand(1000)], {}

    def add(name, tsd, ltr, right_tsd=None):
        parts.append(tsd)
        start = sum(map(len, parts)) + 1
        parts.append(ltr + INT + ltr)
        spans[name] = (start, start + len(ltr + INT + ltr) - 1)
        parts.append(right_tsd if right_tsd is not None else tsd)
        parts.append(rand(800))

    add("A", "AATGC", LTR)
    add("B", "TAAA", LTR_B, right_tsd="TAAAT")    # the T after it makes TAAAT 1 bp in
    add("D", "GGCAT", LTR, right_tsd="CTTAG")
    return "".join(parts), spans


def row(seq_id, length, l5e, l3s, cigar, motif, tsd=".", off="NA"):
    v = dict.fromkeys(COLUMNS, "NA")
    v.update(seq_id=seq_id, seq_len=str(length), status="pass", ltr5_start="1",
             ltr5_end=str(l5e), ltr3_start=str(l3s), ltr3_end=str(length),
             ltr5_len=str(l5e), ltr3_len=str(length - l3s + 1), flank5_len="0",
             flank3_len="0", aln_len="400", bitscore="800", cigar=cigar, motif=motif,
             orientation="+", tsd=tsd, tsd_offset=off, tsd_input=tsd)
    return "\t".join(v[c] for c in COLUMNS)


@pytest.fixture
def snapped(tmp_path, k2l_api):
    seq, spans = build()
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1\n" + seq + "\n")
    (a0, a1), (b0, b1), (d0, d1) = spans["A"], spans["B"], spans["D"]
    n = a1 - a0 + 1
    calls = {
        # A called 1 bp wide on the left: its pairing kept the base, which matches
        # the last internal base, so Kmer2LTR found the TSD only 1 bp in.
        "wide": (f"chr1:{a0 - 1}-{a1}", row(f"chr1:{a0 - 1}-{a1}", n + 1, 401, n - 399, "401=",
                                            "ct...ca", "AATGC", "1,0")),
        # A called again, 2 bp short on the right: it moves onto the same ends.
        "short": (f"chr1:{a0}-{a1 - 2}", row(f"chr1:{a0}-{a1 - 2}", n - 2, 400, n - 399,
                                             "398=2I", "tg...gt")),
        # B is canonical as called, but with a 4-bp TSD Kmer2LTR does not look for.
        "b": (f"chr1:{b0}-{b1}", row(f"chr1:{b0}-{b1}", b1 - b0 + 1, 400, b1 - b0 - 398,
                                     "400=", "tg...ca", "TAAAT", "1,0")),
        "d": (f"chr1:{d0}-{d1}", row(f"chr1:{d0}-{d1}", d1 - d0 + 1, 400, d1 - d0 - 398,
                                     "400=", "tg...ca")),
    }
    tsv = tmp_path / "k.tsv"
    tsv.write_text("#" + "\t".join(COLUMNS) + "\n" + "\n".join(r for _, r in calls.values())
                   + "\n")
    fa = tmp_path / "b.fa"
    with open(fa, "w") as fh:
        for name, _ in calls.values():
            s, e = (int(x) for x in name.split(":")[1].split("-"))
            fh.write(f">{name}\n{seq[s - 1:e]}\n")
    moved, removed, counts = detect.snap_to_canonical_ends(
        str(tsv), str(fa), str(genome), TOOLS, 1, 3e-8)
    rows = {r["seq_id"]: r for r in detect.k2l.read_rows(str(tsv))}
    recs = dict(detect.iter_fasta(str(fa)))
    return seq, spans, calls, moved, removed, counts, rows, recs


def test_a_call_one_bp_wide_moves_onto_its_canonical_ends(snapped):
    seq, spans, calls, moved, _removed, _counts, rows, recs = snapped
    a0, a1 = spans["A"]
    new = f"chr1:{a0}-{a1}"
    assert moved == {calls["wide"][0]: new}
    r = rows[new]
    assert (r["tsd"], r["tsd_offset"], r["motif"]) == ("AATGC", "0,0", "tg...ca")
    assert (r["ltr5_len"], r["ltr3_len"], r["seq_len"]) == ("400", "400", str(a1 - a0 + 1))
    assert r["tsd_input"] == "AATGC"
    assert recs[new] == seq[a0 - 1:a1]
    assert calls["wide"][0] not in rows and calls["wide"][0] not in recs


def test_a_second_call_reaching_the_same_ends_is_removed(snapped):
    _seq, _spans, calls, _moved, removed, counts, rows, recs = snapped
    assert removed == {calls["short"][0]}
    assert calls["short"][0] not in rows and calls["short"][0] not in recs
    assert counts["duplicate"] == 1 and counts["moved"] == 1


def test_canonical_ends_keep_their_place_and_report_their_own_tsd(snapped):
    _seq, _spans, calls, _moved, _removed, counts, rows, _recs = snapped
    r = rows[calls["b"][0]]
    assert (r["tsd"], r["tsd_offset"]) == ("TAAA", "0,0")
    assert counts["canonical"] == 1


def test_a_call_with_no_canonical_ends_nearby_is_untouched(snapped):
    _seq, _spans, calls, _moved, _removed, _counts, rows, _recs = snapped
    assert "\t".join(rows[calls["d"][0]][c] for c in COLUMNS) == calls["d"][1]
