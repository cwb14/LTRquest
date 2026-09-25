"""Tandem-array calls: purged from the _clean_ tables alongside FP families.

A call cut from a tandem array (rDNA, a satellite, a tandem segmental
duplication) has "LTRs" that are two slices of one longer direct repeat, so the
repeat carries on past both of its termini. A genuine LTR-RT's LTR homology
stops at its termini: what lies outside is unrelated target-site sequence.
"""

import gzip
import random

import pytest

from ltrquest import flag_fp

LTR = 400
HEADER = "#seq_id\tseq_len\tstatus\tltr5_start\tltr5_end\tltr3_start\tltr3_end\tnest_status\n"


def rand(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(seq, rate, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT".replace(b, "")) if rng.random() < rate else b for b in seq)


def array_candidate(unit_len=3000, copies=5, flank=5000, seed=1):
    """A tandem array of `copies` units, and the call spanning unit 2's first
    LTR-length slice to unit 3's: 1-based (s5, e5, s3, e3)."""
    unit = rand(unit_len, seed)
    seq = rand(flank, seed + 1) + unit * copies + rand(flank, seed + 2)
    s5 = flank + unit_len + 1
    s3 = s5 + unit_len
    return seq, (s5, s5 + LTR - 1, s3, s3 + LTR - 1)


def element_candidate(internal=3000, flank=5000, seed=10):
    """One LTR-RT inserted into unrelated sequence, 5 bp TSD either side."""
    ltr, tsd = rand(LTR, seed), rand(5, seed + 1)
    seq = (rand(flank, seed + 2) + tsd + ltr + rand(internal, seed + 3)
           + mutate(ltr, 0.02, seed + 4) + tsd + rand(flank, seed + 5))
    s5 = flank + 6
    s3 = s5 + LTR + internal
    return seq, (s5, s5 + LTR - 1, s3, s3 + LTR - 1)


def row(chrom, s5, e5, s3, e3, cls="LTR/unknown/unknown", nest="."):
    """One depth-table row: Kmer2LTR-bounded, so ltr5_start == 1, ltr3_end == seq_len."""
    n = e3 - s5 + 1
    return (f"{chrom}:{s5}-{e3}#{cls}\t{n}\tpass\t1\t{e5 - s5 + 1}\t{s3 - s5 + 1}\t{n}"
            f"\t{nest}\n")


# ---------------------------------------------------------------- the score


def test_a_unit_inside_a_tandem_array_continues_on_both_sides():
    seq, c = array_candidate()
    up, dn = flag_fp.tandem_flank_scores(seq, *c)
    assert up > 0.9 and dn > 0.9


def test_a_solitary_ltr_rt_does_not_continue_on_either_side():
    seq, c = element_candidate()
    up, dn = flag_fp.tandem_flank_scores(seq, *c)
    assert up < 0.1 and dn < 0.1


def test_two_elements_sharing_a_middle_ltr_are_each_one_sided():
    # LTR-int-LTR-int-LTR: a genuine tandem insertion. Each element continues
    # into its partner on one side only, so neither is an array.
    ltr, body, flank = rand(LTR, 20), rand(3000, 21), 5000
    seq = rand(flank, 22) + ltr + body + ltr + mutate(body, 0.02, 23) + ltr + rand(flank, 24)
    s1 = flank + 1
    s2 = s1 + LTR + 3000
    s3 = s2 + LTR + 3000
    first = flag_fp.tandem_flank_scores(seq, s1, s1 + LTR - 1, s2, s2 + LTR - 1)
    second = flag_fp.tandem_flank_scores(seq, s2, s2 + LTR - 1, s3, s3 + LTR - 1)
    # 2% divergence leaves 0.98**15 ~ 74% of 15-mers intact
    assert first[0] < 0.1 and first[1] > 0.5
    assert second[0] > 0.5 and second[1] < 0.1


def test_units_that_differ_in_length_next_to_the_ltr_still_continue():
    # rDNA spacers vary in subrepeat copy number between units, so the flank
    # sits hundreds of bp out of register with the matching internal sequence.
    core, sub, ltr = rand(2500, 30), rand(300, 31), rand(LTR, 32)
    units = [core + sub * n + ltr for n in (1, 3, 2, 1, 3)]
    flank = rand(5000, 33)
    seq = flank + "".join(units) + rand(5000, 34)
    s5 = len(flank) + len(units[0]) - LTR + 1          # unit 1's trailing LTR
    s3 = s5 + len(units[1])                            # unit 2's trailing LTR
    up, dn = flag_fp.tandem_flank_scores(seq, s5, s5 + LTR - 1, s3, s3 + LTR - 1)
    assert up > 0.5 and dn > 0.5


def test_masked_flanks_leave_nothing_to_judge():
    seq, c = array_candidate()
    s5, _, _, e3 = c
    seq = "N" * (s5 - 1) + seq[s5 - 1:e3] + "N" * (len(seq) - e3)
    assert flag_fp.tandem_flank_scores(seq, *c) == (0.0, 0.0)


def test_a_call_at_the_sequence_ends_is_scored_without_error():
    seq, (s5, e5, s3, e3) = array_candidate(flank=0)
    shift = s5 - 1
    seq = seq[shift:e3]                                # cut flush to both LTRs
    up, dn = flag_fp.tandem_flank_scores(seq, 1, e5 - shift, s3 - shift, e3 - shift)
    assert up == 0.0 and dn == 0.0


def test_soft_masked_sequence_scores_like_upper_case():
    seq, c = array_candidate()
    assert flag_fp.tandem_flank_scores(seq.lower(), *c) == flag_fp.tandem_flank_scores(seq, *c)


def test_no_internal_region_means_no_verdict():
    seq, (s5, e5, _, _) = array_candidate()
    assert flag_fp.tandem_flank_scores(seq, s5, e5, e5 + 1, e5 + LTR) == (0.0, 0.0)


# ---------------------------------------------------------------- the scan


@pytest.fixture
def two_contigs(tmp_path):
    arr, arr_c = array_candidate()
    ele, ele_c = element_candidate()
    genome = tmp_path / "g.fa.gz"                      # gzip, soft-masked, described
    with gzip.open(genome, "wt") as fh:
        fh.write(f">arr\n{arr}\n>ele desc\n{ele.lower()}\n")
    table = tmp_path / "g_depth0_ltr.tsv"
    table.write_text(HEADER + row("arr", *arr_c) + row("ele", *ele_c))
    arr_id = f"arr:{arr_c[0]}-{arr_c[3]}"
    ele_id = f"ele:{ele_c[0]}-{ele_c[3]}"
    return tmp_path, genome, table, arr_id, ele_id


def test_the_scan_flags_the_array_unit_and_not_the_element(two_contigs):
    tmp, genome, table, arr_id, _ = two_contigs
    report = tmp / "tandem.tsv"
    flagged, judged = flag_fp.tandem_array_coords([str(table)], [str(genome)], 0.5, str(report))
    assert flagged == {arr_id} and judged == 2
    verdicts = {r.split("\t")[0].split(":")[0]: r.split("\t")[-1]
                for r in report.read_text().splitlines()[1:]}
    assert verdicts == {"arr": "tandem", "ele": "keep"}


def test_a_threshold_above_every_score_flags_nothing(two_contigs):
    _, genome, table, _, _ = two_contigs
    assert flag_fp.tandem_array_coords([str(table)], [str(genome)], 1.01) == (set(), 2)


def test_a_call_on_a_sequence_in_no_genome_is_left_alone(two_contigs):
    _, genome, table, _, _ = two_contigs
    with open(table, "a") as fh:
        fh.write(row("elsewhere", 5001, 5400, 8001, 8400))
    flagged, judged = flag_fp.tandem_array_coords([str(table)], [str(genome)], 0.5)
    assert judged == 2 and all(not c.startswith("elsewhere") for c in flagged)


# ---------------------------------------------------------------- Stage B


def flag_fp_argv(tmp, table, *extra):
    ids = [line.split("\t")[0] for line in table.read_text().splitlines()[1:]]
    (tmp / "cons.tsv").write_text("".join(f"{ids[0]}\t{i}\n" for i in ids))
    (tmp / "int.tsv").write_text("".join(f"{ids[0]}\t{i}\n" for i in ids))
    (tmp / "ltrs.fa").write_text(">rep\nACGT\n")
    return ["--consensus-cluster", str(tmp / "cons.tsv"), "--internal-cluster", str(tmp / "int.tsv"),
            "--ltr-fasta", str(tmp / "ltrs.fa"), "--domains-tsv", str(table),
            "--no-plot", "-o", str(tmp / "out"), *extra]


def test_stage_b_purges_the_array_unit_and_scrubs_references_to_it(two_contigs, capsys):
    tmp, genome, table, arr_id, ele_id = two_contigs
    # have the element name the array unit as nested in it, to see the token go
    lines = table.read_text().splitlines(keepends=True)
    lines[2] = lines[2].rsplit("\t", 1)[0] + f"\tnest-outer:{arr_id}\n"
    table.write_text("".join(lines))

    assert flag_fp.main(flag_fp_argv(tmp, table, "--tandem-genome", str(genome))) == 0
    clean = (tmp / "g_depth0_clean_ltr.tsv").read_text()
    assert arr_id not in clean and ele_id in clean and "nest-outer" not in clean
    err = capsys.readouterr().err
    assert "FP fraction: 0/2" in err                   # a purged array is not an FP family
    assert (tmp / "out.tandem.tsv").is_file()


def test_without_a_tandem_genome_nothing_is_purged_for_being_an_array(two_contigs):
    tmp, _, table, arr_id, _ = two_contigs
    assert flag_fp.main(flag_fp_argv(tmp, table)) == 0
    assert arr_id in (tmp / "g_depth0_clean_ltr.tsv").read_text()
    assert not (tmp / "out.tandem.tsv").exists()


def test_a_tandem_genome_needs_the_depth_tables(two_contigs):
    tmp, genome, table, _, _ = two_contigs
    argv = [a for a in flag_fp_argv(tmp, table, "--tandem-genome", str(genome))
            if a not in ("--domains-tsv", str(table))]
    with pytest.raises(SystemExit) as excinfo:
        flag_fp.main(argv)
    assert "--domains-tsv" in str(excinfo.value)
