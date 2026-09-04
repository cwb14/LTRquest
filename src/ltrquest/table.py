"""Header-driven access to LTRquest's element tables.

Every element table carries a '#'-prefixed header, so a column can be found by
name. Reading by position instead is what makes a schema change a silent
corruption rather than an error, and the schema does change: the columns come
from Kmer2LTR, whose output has grown before and will again.

Two sentinels travel through these tables and mean different things. 'NA' is a
measurement that could not be made; '.' is a search that ran and found nothing.
Code that only needs to know whether a value is usable can ask `is_missing`,
but code that reports to the user should keep them apart.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

MISSING = frozenset({"", ".", "NA"})


def is_missing(value: Optional[str]) -> bool:
    return value is None or value.strip() in MISSING


def parse_header(line: str) -> list:
    """Column names from a header line, or [] if this is not one."""
    if not line.startswith("#"):
        return []
    return line[1:].rstrip("\n").split("\t")


@dataclass(frozen=True)
class Columns:
    names: list

    @classmethod
    def of(cls, names: Sequence[str]) -> "Columns":
        return cls(list(names))

    @classmethod
    def from_line(cls, line: str) -> "Columns":
        return cls(parse_header(line))

    def __contains__(self, name: str) -> bool:
        return name in self.names

    def index(self, name: str) -> Optional[int]:
        try:
            return self.names.index(name)
        except ValueError:
            return None

    def require(self, name: str) -> int:
        i = self.index(name)
        if i is None:
            raise KeyError(f"no column {name!r}; table has: {', '.join(self.names)}")
        return i

    def get(self, row: Sequence[str], name: str, default: str = ".") -> str:
        i = self.index(name)
        if i is None or i >= len(row):
            return default
        return row[i]


def as_int(value: Optional[str], default: Optional[int] = None) -> Optional[int]:
    if is_missing(value):
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def as_float(value: Optional[str], default: Optional[float] = None) -> Optional[float]:
    if is_missing(value):
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default
