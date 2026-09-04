import pytest

from ltrquest.table import Columns, as_float, as_int, is_missing, parse_header

HEADER = "#seq_id\tseq_len\tstatus\tk2p\ttsd\n"


def test_parse_header_strips_hash_and_newline():
    assert parse_header(HEADER) == ["seq_id", "seq_len", "status", "k2p", "tsd"]


def test_parse_header_rejects_a_data_line():
    assert parse_header("chr1:1-2\t100\tpass\t0.01\tAAAAT\n") == []


def test_index_and_membership():
    cols = Columns.from_line(HEADER)
    assert cols.index("status") == 2
    assert cols.index("absent") is None
    assert "k2p" in cols
    assert "absent" not in cols


def test_require_names_what_is_available():
    cols = Columns.from_line(HEADER)
    with pytest.raises(KeyError) as excinfo:
        cols.require("K2P_d")
    assert "seq_id" in str(excinfo.value)


def test_get_returns_default_for_short_rows():
    cols = Columns.from_line(HEADER)
    assert cols.get(["chr1:1-2", "100"], "tsd") == "."
    assert cols.get(["chr1:1-2", "100"], "tsd", default="") == ""


def test_get_returns_default_for_unknown_column():
    cols = Columns.from_line(HEADER)
    assert cols.get(["a", "b", "c", "d", "e"], "nope") == "."


def test_missing_covers_both_sentinels():
    assert is_missing("NA") and is_missing(".") and is_missing("")
    assert not is_missing("0")


def test_as_int_and_as_float_tolerate_sentinels():
    assert as_int("NA") is None
    assert as_int(".", default=0) == 0
    assert as_int("12") == 12
    assert as_int("12.0") == 12
    assert as_float("NA") is None
    assert as_float("0.5") == 0.5
    assert as_float("not-a-number", default=-1.0) == -1.0
