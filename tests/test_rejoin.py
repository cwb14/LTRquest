"""Rejoining a host whose LTR pair a tenant split.

A tenant that inserts at position p of a host's LTR leaves `[1..p] TENANT
[p+1..L]`. LTRharvest and LTR_FINDER pair the two `[1..p]` fragments, so the host
comes out truncated at the tenant and the tenant looks free beside it. The
reconciler rebuilds the host when the bases beyond the tenant continue the host's
own LTR, starting at the junction. These tests plant both orientations of that
split in a synthetic genome, next to an unrelated abutting pair, an isolated
element and a split whose remainder another call already claims.
"""

from __future__ import annotations

import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import pytest

from ltrquest import kmer2ltr
from ltrquest import reconcile
from ltrquest.annotate import REJOIN_SIDECAR_SUFFIX, load_pass2_links, load_tesorter_strands
from ltrquest.detect import DETECT_COLUMNS, iter_fasta, revcomp

LTR = 400          # host LTR length
P = 250            # insertion offset inside the host LTR
TENANT_TSD = 5     # the tenant duplicates the 5 host bases before its insertion site
CHROM = "chr1"


@pytest.fixture(scope="module")
def tools_dir() -> str:
    """A directory Kmer2LTR's `genome` module imports from; skip without one."""
    path = os.environ.get("LTRQUEST_TOOLS_DIR", ".")
    try:
        kmer2ltr.genome_module(path)
    except RuntimeError as exc:
        pytest.skip(f"Kmer2LTR unavailable: {exc}")
    return path


@dataclass
class Call:
    """One element as a detection round reports it: genomic span and LTR inner ends."""
    s: int
    e: int
    l1: int
    r0: int
    orientation: str = "+"
    tsd: str = "."

    @property
    def key(self) -> str:
        return f"{CHROM}:{self.s}-{self.e}"


class Builder:
    def __init__(self, seed: int):
        self.r = random.Random(seed)
        self.seq: List[str] = []
        self.pos = 0                      # bases written so far

    def rand(self, n: int) -> str:
        return "".join(self.r.choice("ACGT") for _ in range(n))

    def ltr(self, n: int) -> str:
        return "TG" + self.rand(n - 4) + "CA"

    def mutate(self, s: str, rate: float = 0.02, keep: int = 12) -> str:
        """Substitutions away from both ends, so the planted termini stay exact."""
        out = list(s)
        for i in range(keep, len(s) - keep):
            if self.r.random() < rate:
                out[i] = self.r.choice([b for b in "ACGT" if b != s[i]])
        return "".join(out)

    def tsd(self) -> str:
        while True:
            t = self.rand(5)
            if len(set(t)) > 2:
                return t

    def add(self, s: str) -> int:
        """Append `s`; return its first base, 1-based."""
        first = self.pos + 1
        self.seq.append(s)
        self.pos += len(s)
        return first

    def element(self, ltr_len: int = 300, internal: int = 1000) -> Call:
        """A complete element with its own TSD."""
        t = self.tsd()
        self.add(t)
        left = self.ltr(ltr_len)
        s = self.add(left)
        r0 = self.add(self.rand(internal)) + internal
        self.add(self.mutate(left))
        e = self.pos
        self.add(t)
        return Call(s, e, s + ltr_len - 1, r0, tsd=t)

    def tenant(self) -> str:
        left = self.ltr(300)
        return left + self.rand(1000) + self.mutate(left)

    @property
    def genome(self) -> str:
        return "".join(self.seq)


@dataclass
class Split:
    """A planted host with a tenant in one LTR: what the callers report, and the truth."""
    host_called: Call
    tenant: Call
    true_s: int
    true_e: int
    true_l1: int
    true_r0: int
    host_tsd: str


def plant_right(b: Builder, internal: int = 1500) -> Split:
    """Tenant in the host's genomic-right LTR; the callers pair the [1..P] fragments."""
    host_tsd = b.tsd()
    b.add(host_tsd)
    left = b.ltr(LTR)
    s = b.add(left)
    b.add(b.rand(internal))
    right = b.mutate(left)
    r0 = b.add(right[:P])
    ten = b.tenant()
    ts = b.add(ten)
    b.add(right[P - TENANT_TSD:P] + right[P:])
    e = b.pos
    b.add(host_tsd)
    called = Call(s, r0 + P - 1, s + P - 1, r0)
    tenant = Call(ts, ts + len(ten) - 1, ts + 299, ts + len(ten) - 300,
                  tsd=right[P - TENANT_TSD:P])
    return Split(called, tenant, s, e, s + LTR - 1, r0, host_tsd)


def plant_left(b: Builder, internal: int = 1500, orientation: str = "+") -> Split:
    """Tenant in the host's genomic-left LTR; the callers pair the far-side fragments."""
    host_tsd = b.tsd()
    b.add(host_tsd)
    left = b.ltr(LTR)
    s = b.add(left[:P])
    ten = b.tenant()
    ts = b.add(ten)
    after = b.add(left[P - TENANT_TSD:P] + left[P:])
    b.add(b.rand(internal))
    right = b.mutate(left)
    r0 = b.add(right)
    e = b.pos
    b.add(host_tsd)
    frag = LTR - P + TENANT_TSD                 # the fragment both called LTRs share
    called = Call(after, e, after + frag - 1, e - frag + 1, orientation=orientation)
    tenant = Call(ts, ts + len(ten) - 1, ts + 299, ts + len(ten) - 300,
                  tsd=left[P - TENANT_TSD:P])
    return Split(called, tenant, s, e, after + frag - 1, r0, host_tsd)


def row(call: Call, genome: str, cls: str = "LTR/Gypsy/Tekay") -> List[str]:
    """A Kmer2LTR-bounded round-table row for `call`, every column filled."""
    n = call.e - call.s + 1
    l5 = call.l1 - call.s + 1
    l3s = call.r0 - call.s + 1
    seq = genome[call.s - 1:call.e]
    vals: Dict[str, str] = {
        "seq_id": f"{call.key}#{cls}", "seq_len": str(n), "status": "pass",
        "ltr5_start": "1", "ltr5_end": str(l5), "ltr3_start": str(l3s), "ltr3_end": str(n),
        "ltr5_len": str(l5), "ltr3_len": str(n - l3s + 1), "flank5_len": "0",
        "flank3_len": "0", "aln_len": str(l5), "n_sites": str(l5), "n_ts": "2", "n_tv": "1",
        "n_gapcols": "0", "identity": "0.99", "p_dist": "0.01", "k2p": "0.0101",
        "k2p_se": "0.003", "bitscore": "500", "flank_margin_bits": "NA",
        "cigar": f"{l5}=", "motif": f"{seq[:2]}...{seq[-2:]}".lower(), "k2p_time": "168333",
        "orientation": call.orientation, "tsd": call.tsd,
        "tsd_offset": "0,0" if call.tsd != "." else "NA", "tsd_input": call.tsd,
        "domains": ".", "nest_status": ".",
    }
    return [vals[c] for c in DETECT_COLUMNS]


def write_round(tmp: Path, genome: str, calls: List[Call]) -> tuple:
    gpath = tmp / "g.fa"
    gpath.write_text(f">{CHROM}\n" + "\n".join(genome[i:i + 60]
                                               for i in range(0, len(genome), 60)) + "\n")
    tsv, fa = tmp / "x_r1_ltr.tsv", tmp / "x_r1_ltr.fa"
    rows = [row(c, genome) for c in calls]
    tsv.write_text("#" + "\t".join(DETECT_COLUMNS) + "\n"
                   + "".join("\t".join(r) + "\n" for r in rows))
    with fa.open("w") as fh:
        for c, r in zip(calls, rows):
            seq = genome[c.s - 1:c.e]
            fh.write(f">{r[0]}\n{seq if c.orientation == '+' else revcomp(seq)}\n")
    return gpath, tsv, fa


def run_reconcile(tmp: Path, tsv: Path, fa: Path, genome=None, tools=None) -> Path:
    out = tmp / "out" / "x"
    out.parent.mkdir(exist_ok=True)
    argv = ["--out-prefix", str(out), "--tsv", str(tsv), "--fa", str(fa)]
    if genome is not None:
        argv += ["--genome", str(genome)]
    if tools is not None:
        argv += ["--tools-dir", tools]
    reconcile.main(argv)
    return out


def load_outputs(out: Path):
    rows, seqs = {}, {}
    for t in sorted(out.parent.glob(out.name + "_depth*_ltr.tsv")):
        depth = int(t.name.split("_depth")[1].split("_")[0])
        names = None
        for line in t.read_text().splitlines():
            if line.startswith("#"):
                names = line.lstrip("#").split("\t")
                continue
            cols = dict(zip(names, line.split("\t")))
            cols["depth"] = depth
            rows[cols["seq_id"].split("#")[0]] = cols
    for f in out.parent.glob(out.name + "_depth*_ltr.fa"):
        for h, s in iter_fasta(str(f)):
            seqs[h.split("#")[0]] = s
    return rows, seqs


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    b = Builder(seed=7)
    b.add(b.rand(3000))
    right = plant_right(b)
    b.add(b.rand(3000))
    left = plant_left(b, orientation="-")
    b.add(b.rand(3000))
    d1 = b.element()
    d2 = b.element()                  # abuts d1 across two TSDs: an unrelated neighbour
    b.add(b.rand(3000))
    lone = b.element()
    b.add(b.rand(3000))
    genome = b.genome
    calls = [right.host_called, right.tenant, left.host_called, left.tenant, d1, d2, lone]
    tmp = tmp_path_factory.mktemp("rejoin")
    return dict(genome=genome, right=right, left=left, decoy=(d1, d2), lone=lone,
                files=write_round(tmp, genome, calls), tmp=tmp)


@pytest.fixture(scope="module")
def rejoined(planted, tools_dir):
    gpath, tsv, fa = planted["files"]
    out = run_reconcile(planted["tmp"], tsv, fa, genome=gpath, tools=tools_dir)
    return out, *load_outputs(out)


def test_a_tenant_in_the_right_ltr_rebuilds_the_host_to_its_true_ends(planted, rejoined):
    _, rows, _ = rejoined
    sp = planted["right"]
    key = f"{CHROM}:{sp.true_s}-{sp.true_e}"
    assert key in rows, sorted(rows)
    host = rows[key]
    assert sp.host_called.key not in rows
    assert host["depth"] == 1
    assert f"nest-outer:{sp.tenant.key}" in host["nest_status"]
    assert rows[sp.tenant.key]["nest_status"] == f"nest-inner:{key}"
    assert int(host["ltr5_end"]) == sp.true_l1 - sp.true_s + 1
    assert int(host["ltr3_start"]) == sp.true_r0 - sp.true_s + 1
    assert int(host["ltr3_end"]) == int(host["seq_len"]) == sp.true_e - sp.true_s + 1
    # LTR bases only: the tenant is not part of the LTR it sits in, its TSD copy is
    assert int(host["ltr5_len"]) == LTR
    assert int(host["ltr3_len"]) == LTR + TENANT_TSD


def test_a_tenant_in_the_left_ltr_rebuilds_the_host_on_the_minus_strand(planted, rejoined):
    _, rows, _ = rejoined
    sp = planted["left"]
    key = f"{CHROM}:{sp.true_s}-{sp.true_e}"
    assert key in rows, sorted(rows)
    host = rows[key]
    assert host["depth"] == 1 and host["orientation"] == "-"
    assert f"nest-outer:{sp.tenant.key}" in host["nest_status"]
    assert int(host["ltr5_len"]) == LTR + TENANT_TSD
    assert int(host["ltr3_len"]) == LTR
    assert int(host["ltr3_start"]) == sp.true_r0 - sp.true_s + 1


def test_the_rebuilt_ends_carry_the_hosts_own_tsd(planted, rejoined):
    _, rows, _ = rejoined
    for side in ("right", "left"):
        sp = planted[side]
        host = rows[f"{CHROM}:{sp.true_s}-{sp.true_e}"]
        assert host["tsd"] == host["tsd_input"] == sp.host_tsd
        assert host["tsd_offset"] == "0,0"
        g = planted["genome"]
        ends = g[sp.true_s - 1:sp.true_s + 1], g[sp.true_e - 2:sp.true_e]
        assert host["motif"] == "{}...{}".format(*ends).lower()


def test_the_divergence_measured_on_the_called_pair_is_kept(planted, rejoined):
    _, rows, _ = rejoined
    sp = planted["right"]
    host = rows[f"{CHROM}:{sp.true_s}-{sp.true_e}"]
    assert (host["k2p"], host["aln_len"], host["cigar"]) == ("0.0101", str(P), f"{P}=")


def test_the_host_record_is_the_genome_with_its_tenant_painted(planted, rejoined):
    _, _, seqs = rejoined
    g = planted["genome"]
    for side in ("right", "left"):
        sp = planted[side]
        rec = seqs[f"{CHROM}:{sp.true_s}-{sp.true_e}"]
        want = list(g[sp.true_s - 1:sp.true_e])
        for i in range(sp.tenant.s - sp.true_s, sp.tenant.e - sp.true_s + 1):
            want[i] = reconcile.IUPAC_DEPTH_SEQ[0]
        want = "".join(want)
        if sp.host_called.orientation == "-":
            want = revcomp(want)
        assert rec == want


def test_unrelated_neighbours_and_a_lone_element_are_untouched(planted, rejoined):
    _, rows, _ = rejoined
    for c in (*planted["decoy"], planted["lone"]):
        assert rows[c.key]["depth"] == 0
        assert rows[c.key]["nest_status"] == "."


def test_every_table_row_has_its_fasta_record(rejoined):
    _, rows, seqs = rejoined
    assert set(rows) == set(seqs)
    for k, r in rows.items():
        assert len(seqs[k]) == int(r["seq_len"])


def test_the_sidecar_maps_each_old_host_to_its_new_key(planted, rejoined):
    out, _, _ = rejoined
    side = Path(str(out) + REJOIN_SIDECAR_SUFFIX)
    lines = side.read_text().splitlines()
    assert lines[0].startswith("#old_key\tnew_key\ttenant\tltr")
    got = {tuple(line.split("\t")[:4]) for line in lines[1:]}
    r, left = planted["right"], planted["left"]
    assert got == {
        (r.host_called.key, f"{CHROM}:{r.true_s}-{r.true_e}", r.tenant.key, "right"),
        (left.host_called.key, f"{CHROM}:{left.true_s}-{left.true_e}", left.tenant.key, "left"),
    }


def test_without_a_genome_nothing_is_rejoined_and_the_sidecar_is_empty(planted, tmp_path, capsys):
    _, tsv, fa = planted["files"]
    out = run_reconcile(tmp_path, tsv, fa)
    rows, _ = load_outputs(out)
    assert planted["right"].host_called.key in rows
    assert all(r["depth"] == 0 for r in rows.values())
    assert "--genome" in capsys.readouterr().err
    lines = Path(str(out) + REJOIN_SIDECAR_SUFFIX).read_text().splitlines()
    assert len(lines) == 1 and lines[0].startswith("#old_key")


def test_a_remainder_another_call_already_claims_is_left_alone(tmp_path, tools_dir):
    b = Builder(seed=11)
    b.add(b.rand(3000))
    sp = plant_right(b)
    b.add(b.rand(3000))
    genome = b.genome
    # a call sitting on the host's LTR remainder beyond the tenant
    squatter = Call(sp.tenant.e + 20, sp.true_e + 1500, sp.tenant.e + 120, sp.true_e + 1400)
    gpath, tsv, fa = write_round(tmp_path, genome, [sp.host_called, sp.tenant, squatter])
    out = run_reconcile(tmp_path, tsv, fa, genome=gpath, tools=tools_dir)
    rows, _ = load_outputs(out)
    assert sp.host_called.key in rows
    assert rows[sp.tenant.key]["nest_status"] == "."


def test_strand_loaders_follow_a_rejoined_host_to_its_new_key(tmp_path):
    work = tmp_path / "x_r1.work" / "minimap2_pass2"
    work.mkdir(parents=True)
    (tmp_path / "x_r1.work" / "x_r1.bounded.fa.cls.tsv").write_text(
        "#TE\tOrder\tSuperfamily\tClade\tComplete\tStrand\tDomains\n"
        "chr1:100-500\tLTR\tGypsy\tTekay\tno\t-\tRT|Tekay\n")
    (work / "pass2.tsv").write_text("chr1:100-500\tpass\t0.9\t0.9\t0.9\tchr2:1-900\n")
    (work / "pass2.paf").write_text(
        "chr1:100-500\t401\t0\t400\t-\tchr2:1-900\t900\t0\t400\t390\t400\t60\n")
    assert "chr1:100-9000" not in load_tesorter_strands("x", str(tmp_path))
    (tmp_path / f"x{REJOIN_SIDECAR_SUFFIX}").write_text(
        "#old_key\tnew_key\ttenant\tltr\tremainder_bp\tidentity\tevalue\ttsd\n"
        "chr1:100-500\tchr1:100-9000\tchr1:501-8500\tright\t500\t0.97\t1e-80\tACGTA\n")
    strands = load_tesorter_strands("x", str(tmp_path))
    assert strands["chr1:100-9000"] == strands["chr1:100-500"] == "-"
    target_of, orientation = load_pass2_links("x", str(tmp_path))
    assert target_of["chr1:100-9000"] == "chr2:1-900"
    assert orientation[("chr1:100-9000", "chr2:1-900")] == "-"


def test_a_genome_that_lacks_the_calls_sequences_is_refused(planted, tmp_path, tools_dir):
    _, tsv, fa = planted["files"]
    other = tmp_path / "other.fa"
    other.write_text(">chrX\n" + "ACGT" * 100 + "\n")
    with pytest.raises(RuntimeError, match="no sequence named chr1"):
        run_reconcile(tmp_path, tsv, fa, genome=other, tools=tools_dir)


def test_a_continuation_longer_than_any_ltr_is_a_tandem_copy_and_is_refused(tmp_path, tools_dir):
    # A host whose internal region is duplicated in tandem right after its tenant:
    # past the tenant, its own sequence "continues" for far longer than any LTR
    # (as the next unit of an array does), which is no LTR remainder.
    b = Builder(seed=13)
    b.add(b.rand(3000))
    host_tsd = b.tsd()
    b.add(host_tsd)
    left = b.ltr(LTR)
    s = b.add(left)
    internal = b.rand(9000)
    b.add(internal)
    right = b.mutate(left)
    r0 = b.add(right[:P])
    ten = b.tenant()
    ts = b.add(ten)
    b.add(right[P - TENANT_TSD:P] + right[P:] + b.mutate(internal, 0.01))
    b.add(b.rand(3000))
    host = Call(s, r0 + P - 1, s + P - 1, r0)
    tenant = Call(ts, ts + len(ten) - 1, ts + 299, ts + len(ten) - 300)
    gpath, tsv, fa = write_round(tmp_path, b.genome, [host, tenant])
    out = run_reconcile(tmp_path, tsv, fa, genome=gpath, tools=tools_dir)
    rows, _ = load_outputs(out)
    assert host.key in rows
    assert rows[tenant.key]["nest_status"] == "."
