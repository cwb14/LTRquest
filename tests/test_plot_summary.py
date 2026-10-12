"""Size summaries: lengths with nested elements taken out."""

from __future__ import annotations

from ltrquest.plot_summary import parse_alignment_file


def _table(tmp_path, nest_status):
    # element chr1:1001-9000; left LTR 1001-1400, right LTR 3001-9000
    path = tmp_path / "t.tsv"
    path.write_text("#seq_id\tltr5_len\tltr5_end\tltr3_start\tk2p\tnest_status\n"
                    f"chr1:1001-9000#LTR/Gypsy/Tekay\t400\t400\t2001\t0.01\t{nest_status}\n")
    return parse_alignment_file(str(path))[0]


def test_a_tenant_in_the_internal_region_is_taken_out_of_the_internal_length(tmp_path):
    rec = _table(tmp_path, "nest-outer:chr1:1501-2500")
    # spans count end - start, as nested_span_bp does
    assert rec["internal_len"] == 1600 - (2500 - 1501)


def test_a_tenant_in_an_ltr_is_not_taken_out_of_the_internal_length(tmp_path):
    rec = _table(tmp_path, "nest-outer:chr1:3251-8600")
    assert rec["internal_len"] == 1600
    assert rec["full_len"] == (9000 - 1001) - (8600 - 3251)
