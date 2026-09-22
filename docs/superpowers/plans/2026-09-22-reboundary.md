# Family-guided re-boundarying (`--reboundary`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** An opt-in LTRquest stage (pipeline flag `--reboundary` and a post-hoc command `ltrquest-reboundary`) that extends truncated LTR-RT calls to the ends their family's LTR model supports, with Kmer2LTR re-scoring every changed element.

**Architecture:** Per family (pooled over all genomes), build an LTR model whose termini come only from 5′/3′-LTR agreement. Place it on every member (glocal core alignment plus a terminal anchor for large indels) to propose outward-only ends. Apply per-family QC and conflict rules, then re-score each candidate with Kmer2LTR's Python API (`classify` with a family credit), so Kmer2LTR stays the single source of the 29 Kmer2LTR columns. Only the `_clean_` tables and FASTAs change (two-phase atomic commit). A per-genome sidecar doubles as the old→new key map the GFF3 writer reads.

**Tech Stack:** Python ≥3.9 (runs on 3.11 in the `ltrquest` env), parasail-python, pyfaidx, numpy, MAFFT, Kmer2LTR (pinned `aa25f46`), bash driver, Nextflow DSL2, pytest.

**Spec:** `docs/superpowers/specs/2026-09-22-reboundary-design.md` — read it before starting; this plan argues from it.

**Pre-validated:** every Python file and test in Tasks 1–10, and the driver edits in Task 11,
were assembled into a scratch copy of this repo and run against the pinned Kmer2LTR:
**548 passed, 0 failed** (baseline 463 + 85 new), including Kmer2LTR arbitration on the
synthetic fixture and the benchmark harness end to end. So the code in those tasks is known
to work as written; the TDD steps are still worth following, because a test that passes on the
first run has not shown you what it is for. Tasks 12 (Nextflow), 13 (docs), 14 and 15 were not
pre-run.

## Global Constraints

- Python syntax must stay 3.9-compatible: every new module starts with `from __future__ import annotations`; no `match`; ruff line length 100 (ruff is not installed here; keep lines ≤ 100 by hand).
- No new dependencies. Allowed: mafft, parasail-python, pywfa, pyfaidx, numpy (already in `environment.yml`, `recipe/meta.yaml`, Dockerfile).
- **Never install anything** (pip/mamba/conda). If a tool or module is missing, stop and tell the user exactly what is needed.
- Never load a whole genome into memory: genome access is pyfaidx random access only.
- Extend only, never trim. Element ends are never found from a TSD or from TG..CA; TSDs are validation (and family-level QC) only.
- Raw `<prefix>_depth<N>_ltr.{tsv,fa}` tables are detection outputs and are never modified.
- Kmer2LTR columns (the first 29 of every depth table, `kmer2ltr.COLUMNS`) come from Kmer2LTR only.
- Work on branch `reboundary` in `/data2/chris/poa_LTR/LTRquest`. Commit after every task; never push. Every commit message ends with:
  `Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>`
- Test command (used in every task; `PYT` below means exactly this):
  ```bash
  cd /data2/chris/poa_LTR/LTRquest && LTRQUEST_TOOLS_DIR=/data2/chris/poa_LTR/ltrquest_tools \
    PYTHONPATH=src /home/chris/bin/mambaforge/envs/ltrquest/bin/python -m pytest -p no:cacheprovider
  ```
  and every command that runs mafft needs `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH` (prefix it; the ltrquest env holds mafft).
  Baseline before this plan: `PYT -m "not slow"` → 463 passed, 14 deselected.
- Long jobs (> ~10 min) are launched detached: `setsid nohup <cmd> > log 2>&1 < /dev/null &`, then verify the process has its own session (`ps -o pid,ppid,sid,cmd`). Never use the Bash tool's `run_in_background`.
- After running any pipeline/benchmark on real data, write a terse `memo.md` beside the outputs (date | purpose | env | inputs | exact commands | key outputs | gotchas | `claude --resume` id).

## File structure

| path | responsibility |
|---|---|
| `src/ltrquest/kmer2ltr.py` (modify) | `api()`: import Kmer2LTR's Python API with a version guard; `_clone()` factored out of `resolve()` |
| `src/ltrquest/ltr_model.py` (new) | `Member`, `Model`, `Genomes`; alignment primitives (`glocal`); reference selection; consensus model (MAFFT + 5′/3′ agreement + end polishing); subfamily models; family QC maths |
| `src/ltrquest/ltr_place.py` (new) | placing a model on one element: core + terminal anchor → `Proposal`; obstacle report; TSD probe via Kmer2LTR; nearest templates |
| `src/ltrquest/reboundary_io.py` (new) | clean-table loading → `Member`s; FASTA streaming; record splice/paint/sanitize; nest_status re-keying; conflict index; sidecar format; two-phase commit; `rewrite()` |
| `src/ltrquest/reboundary.py` (new) | `Settings`; pooled driver `run()` (family phase, QC, conflicts, Kmer2LTR arbitration, rewrite); post-hoc backup/restore/regenerate; CLI `main()` |
| `src/ltrquest/gff3.py` (modify) | `--reboundary-map`: old-key lookups for families/strand provenance; `boundary_source`, `boundary_shift` attributes |
| `src/ltrquest/scripts/ltrquest.sh` (modify) | `--reboundary`; annotation stage restructured when on; `run_reboundary_stage`; carry-forward/promotion hygiene |
| `pyproject.toml` (modify) | console script `ltrquest-reboundary` |
| `modules/local/ltrquest/reboundary/{main.nf,environment.yml}` (new), `workflows/ltrquest.nf`, `modules/local/ltrquest/gff3/main.nf`, `nextflow.config`, `nextflow_schema.json`, `conf/modules.config` (modify) | Nextflow parity |
| `benchmarks/reboundary/bench.py`, `benchmarks/reboundary/README.md` (new) | B1 planted obstacles, B2 real run, B3 false changes |
| `tests/reboundary_fixtures.py` (new) | synthetic one-contig family with truncated calls, written exactly as LTRquest writes a run |
| `tests/test_reboundary_fixtures.py`, `tests/test_ltr_model.py`, `tests/test_ltr_place.py`, `tests/test_reboundary_io.py`, `tests/test_reboundary.py` (new); `tests/conftest.py`, `tests/test_kmer2ltr.py`, `tests/test_gff3.py`, `tests/test_packaging.py` (modify) | tests |
| `docs/outputs.md`, `README.md`, `CHANGELOG.md` (modify) | docs |

Deviation from the spec's file table, on purpose: the spec lists one `ltr_model.py` and one `reboundary.py`; the plan splits them into `ltr_model.py` + `ltr_place.py` and `reboundary.py` + `reboundary_io.py` so each file has one job and stays readable. One knob changed after design: `anchor_len` candidates are 30/50/80 bp (not 60/80/100) — a local alignment of the model's outer bases can only pin the true end if the insertion lies beyond the anchor, and 30 bp is the TERM identity window already used by the gates.

---

### Task 0: Kmer2LTR checkout for tests and benchmarks

The tests and benchmarks import Kmer2LTR's Python API at the exact commit the pipeline pins. This is a source checkout (what LTRquest itself does at runtime), not an install.

**Files:** none in the repo.

- [ ] **Step 1: Clone Kmer2LTR at the pinned commit**

```bash
mkdir -p /data2/chris/poa_LTR/ltrquest_tools
git clone -q https://github.com/cwb14/Kmer2LTR.git /data2/chris/poa_LTR/ltrquest_tools/Kmer2LTR
git -C /data2/chris/poa_LTR/ltrquest_tools/Kmer2LTR checkout -q aa25f46262a32826c1fabf4a01977c032cd692dd
git -C /data2/chris/poa_LTR/ltrquest_tools/Kmer2LTR log -1 --format='%h %s'
```
Expected: `aa25f46 Take the last tied minimum as a diagonal segment's start`

- [ ] **Step 2: Confirm the checkout imports in the ltrquest env**

```bash
PYTHONPATH=/data2/chris/poa_LTR/ltrquest_tools/Kmer2LTR/src \
  /home/chris/bin/mambaforge/envs/ltrquest/bin/python -c "from kmer2ltr import align, genome, runner; print(genome.PAD, genome.TSD_K)"
```
Expected: `32 (6, 5)`

- [ ] **Step 3: Confirm the baseline suite**

Run: `PYT -m "not slow" -q 2>&1 | tail -1`
Expected: `463 passed, 14 deselected, 2 warnings`

No commit (nothing in the repo changed).

---

### Task 1: `kmer2ltr.api()` — Kmer2LTR's Python API behind a version guard

**Files:**
- Modify: `src/ltrquest/kmer2ltr.py`
- Modify: `tests/conftest.py`
- Test: `tests/test_kmer2ltr.py`

**Interfaces:**
- Consumes: nothing new.
- Produces:
  - `kmer2ltr.Api` (NamedTuple): `classify(seq_id, seq, **kw) -> Result`, `format_row(result) -> str` (29 tab-separated fields), `Window(up, head, tail, down)`, `Options()`, `orient(seq, window) -> Context|None`, `annotate(result, seq, ctx, options) -> Result`, `find_tsd(ctx, b5, b3, ks, shifts) -> (motif, d5, d3)|None`, `PAD: int`, `PROBE: int`, `TSD_K: tuple`, `TSD_SHIFTS: tuple`.
  - `kmer2ltr.api(tools_dir, clone: bool = True) -> Api` (cached per process).
  - pytest fixtures `k2l_api` (skips when unavailable) and `mafft` (skips when mafft is not on PATH).

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_kmer2ltr.py`:

```python
import random
import sys


def test_api_refuses_without_a_checkout_when_cloning_is_off(tmp_path, monkeypatch):
    monkeypatch.setattr(kmer2ltr, "_API", None)
    monkeypatch.setitem(sys.modules, "kmer2ltr", None)   # make `import kmer2ltr` fail
    with pytest.raises(RuntimeError, match="not importable"):
        kmer2ltr.api(tmp_path, clone=False)


def test_api_exposes_the_pieces_reboundary_needs(k2l_api):
    assert callable(k2l_api.classify) and callable(k2l_api.find_tsd)
    assert k2l_api.PAD == 32 and k2l_api.PROBE == 40
    assert 5 in k2l_api.TSD_K and 0 in k2l_api.TSD_SHIFTS


def test_api_classify_formats_to_the_twentynine_columns(k2l_api):
    rng = random.Random(1)
    ltr = "TG" + "".join(rng.choice("ACGT") for _ in range(396)) + "CA"
    internal = "".join(rng.choice("ACGT") for _ in range(1500))
    seq = ltr + internal + ltr
    res = k2l_api.classify(f"chr1:1-{len(seq)}", seq, period_rule="outermost",
                           mutation_rate=3e-8, tsd_credit=0.0)
    row = k2l_api.format_row(res).split("\t")
    assert len(row) == len(kmer2ltr.COLUMNS) == 29
    assert row[2] == "pass" and row[3] == "1" and row[6] == str(len(seq))


def test_api_find_tsd_reads_a_duplication(k2l_api):
    ctx = "GGGGGACGTC" + "T" * 20 + "ACGTCGGGGG"
    hit = k2l_api.find_tsd(ctx, 10, 30, k2l_api.TSD_K, k2l_api.TSD_SHIFTS)
    assert hit is not None and hit[0].endswith("ACGTC")
```

Append to `tests/conftest.py`:

```python
import os
import shutil


@pytest.fixture(scope="session")
def k2l_api():
    """Kmer2LTR's Python API, from the environment or $LTRQUEST_TOOLS_DIR; skip if neither."""
    from ltrquest import kmer2ltr
    try:
        return kmer2ltr.api(os.environ.get("LTRQUEST_TOOLS_DIR", "."), clone=False)
    except RuntimeError as exc:
        pytest.skip(f"Kmer2LTR API unavailable: {exc}")


@pytest.fixture(scope="session")
def mafft() -> str:
    path = shutil.which("mafft")
    if path is None:
        pytest.skip("mafft is not on PATH")
    return path
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `PYT tests/test_kmer2ltr.py -q`
Expected: FAIL — `AttributeError: module 'ltrquest.kmer2ltr' has no attribute '_API'` / `... no attribute 'api'`.

- [ ] **Step 3: Implement `_clone` and `api`**

In `src/ltrquest/kmer2ltr.py`, add `import importlib` and `import inspect` to the imports, and change `from typing import Iterator, Optional, Sequence` to `from typing import Any, Iterator, NamedTuple, Optional, Sequence`.

Replace the clone block inside `resolve()`:

```python
    tools_dir.mkdir(parents=True, exist_ok=True)
    target = tools_dir / "Kmer2LTR"
    if not target.exists():
        try:
            subprocess.run(["git", "clone", "--depth", "1", REPO_URL, str(target)],
                           check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                ...
            ) from e
```

with a call `_clone(tools_dir)` and move the block, unchanged in behaviour and message, into:

```python
def _clone(tools_dir: Path) -> None:
    """git-clone Kmer2LTR into tools_dir/Kmer2LTR unless something is already there."""
    tools_dir = Path(tools_dir)
    tools_dir.mkdir(parents=True, exist_ok=True)
    target = tools_dir / "Kmer2LTR"
    if target.exists():
        return
    try:
        subprocess.run(["git", "clone", "--depth", "1", REPO_URL, str(target)], check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"failed to clone {REPO_URL} into {target}: {e}. Compute nodes "
            f"on this cluster often have no outbound network access -- "
            f"clone it on a login node instead, or install Kmer2LTR "
            f"(pip install kmer2ltr) so it is on PATH. If {target} exists "
            f"from this failed attempt, delete it first: a leftover "
            f"directory there stops the next run from retrying the clone."
        ) from e
```

so `resolve()` reads:

```python
    tools_dir = Path(tools_dir)
    src_dir = tools_dir / "Kmer2LTR" / "src"
    if _importable(src_dir):
        return _module_argv(src_dir)

    _clone(tools_dir)
    if _importable(src_dir):
        return _module_argv(src_dir)
```

Append:

```python
class Api(NamedTuple):
    """The pieces of Kmer2LTR's Python API that per-record stages call.

    The CLI is the right interface for a batch; ltrquest.reboundary needs a
    per-record `tsd_credit`, which only the Python API takes.
    """
    classify: Any
    format_row: Any
    Window: Any
    Options: Any
    orient: Any
    annotate: Any
    find_tsd: Any
    PAD: int
    PROBE: int
    TSD_K: tuple
    TSD_SHIFTS: tuple


_API: Optional[Api] = None


def api(tools_dir, clone: bool = True) -> Api:
    """Kmer2LTR's Python API: an installed package first, else a checkout in tools_dir.

    Refuses a Kmer2LTR whose table or classify signature differs from the one
    this LTRquest was written against -- the same guard `assert_schema` applies
    to the CLI's output, moved to where the skew would enter.
    """
    global _API
    if _API is not None:
        return _API
    try:
        importlib.import_module("kmer2ltr")
    except ImportError:
        src_dir = Path(tools_dir) / "Kmer2LTR" / "src"
        if not _importable(src_dir) and clone:
            _clone(Path(tools_dir))
        if not _importable(src_dir):
            raise RuntimeError(
                f"Kmer2LTR is not importable and {src_dir} holds no checkout. Install "
                f"it (pip install kmer2ltr) or point --tools-dir at a Kmer2LTR clone.")
        sys.path.insert(0, str(src_dir))
        sys.modules.pop("kmer2ltr", None)
    align = importlib.import_module("kmer2ltr.align")
    genome = importlib.import_module("kmer2ltr.genome")
    runner = importlib.import_module("kmer2ltr.runner")
    params = inspect.signature(align._classify).parameters
    missing = [p for p in ("tsd_credit", "period_rule", "mutation_rate") if p not in params]
    if missing or list(runner.COLUMNS) != COLUMNS:
        raise RuntimeError(
            "the Kmer2LTR found does not match this LTRquest "
            f"(classify lacks {missing or 'nothing'}; columns match: "
            f"{list(runner.COLUMNS) == COLUMNS}). Use the commit the Dockerfile pins.")
    _API = Api(align.classify, runner.format_row, genome.Window, genome.Options,
               genome.orient, genome.annotate, genome.find_tsd, genome.PAD, genome.PROBE,
               genome.TSD_K, genome.TSD_SHIFTS)
    return _API
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYT tests/test_kmer2ltr.py -q`
Expected: all pass (the three API tests run because `LTRQUEST_TOOLS_DIR` points at Task 0's checkout).

Run: `PYT -m "not slow" -q 2>&1 | tail -1`
Expected: `467 passed, 14 deselected` (463 + 4).

- [ ] **Step 5: Commit**

```bash
git add src/ltrquest/kmer2ltr.py tests/conftest.py tests/test_kmer2ltr.py
git commit -m "kmer2ltr: expose Kmer2LTR's Python API behind a version guard

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Synthetic run fixture

A tiny run written exactly as LTRquest writes one, with truncated calls whose true ends are known. Every later test uses it.

**Files:**
- Create: `tests/reboundary_fixtures.py`
- Test: `tests/test_reboundary_fixtures.py`

**Interfaces:**
- Consumes: `ltrquest.kmer2ltr.COLUMNS`, `ltrquest.detect.revcomp`.
- Produces: `reboundary_fixtures.build(root: Path, prefix="syn_LTRs", seed=7) -> Fixture`; `Fixture(genome, indir, prefix, elements, ltr)` with `.kind(name) -> Element` and `.family_members()`; `Element` fields `kind, family, strand, true_start, true_end, start, end, l1, r0, depth, tsd, k2p, nest_status`, properties `key`, `name`; constants `CHROM="chrS"`, `FAMILY="syn_fam00001"`, `LTR_LEN=400`, `TAIL`.

- [ ] **Step 1: Write the fixture module**

Create `tests/reboundary_fixtures.py`:

```python
"""A synthetic LTR-RT family with truncated calls, written as LTRquest writes a run.

One contig, `chrS`. Family `syn_fam00001` (LTR 400 bp, TG..CA, internal 2 kb) has
13 intact copies whose calls are right and six whose calls stop short at an
obstacle between their two LTRs -- the way LTRharvest / LTR_FINDER stop:

  del_left         60 bp deleted 50 bp into the left LTR; the call starts after it
  patch_right      30 bp at 30% divergence, 40-70 bp from the right end; the call ends before it
  ins_left         1.5 kb inserted 40 bp into the left LTR; the call starts after it
  del_left_minus   del_left on a minus-strand copy: the obstacle sits at the forward right end
  nested_del_left  del_left inside the internal region of a one-copy host (depth 1)
  conflict_left    del_left with an annotated neighbour overlapping the recoverable bases

Every copy has a 5 bp TSD at its true ends except `decayed`, an intact copy whose
TSD was broken. Tables carry LTRquest's 33 columns; FASTAs store minus-oriented
records reverse-complemented and mask the host's nested inner with `N` (depth 0).
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path
from typing import List

from ltrquest.detect import revcomp as revcomp_record
from ltrquest.kmer2ltr import COLUMNS

CHROM = "chrS"
FAMILY = "syn_fam00001"
LTR_LEN, INT_LEN, SPACER = 400, 2000, 1500
TAIL = ["strand", "family", "domains", "nest_status"]
_COMP = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass
class Element:
    kind: str
    family: str
    strand: str
    true_start: int
    true_end: int
    start: int
    end: int
    l1: int
    r0: int
    depth: int = 0
    tsd: str = "."
    k2p: float = 0.01
    nest_status: str = "."

    @property
    def key(self) -> str:
        return f"{CHROM}:{self.start}-{self.end}"

    @property
    def name(self) -> str:
        return f"{self.key}#LTR/Gypsy/Synth"


@dataclass
class Fixture:
    genome: Path
    indir: Path
    prefix: str
    elements: List[Element]
    ltr: str

    def kind(self, kind: str) -> Element:
        (hit,) = [e for e in self.elements if e.kind == kind]
        return hit

    def family_members(self) -> List[Element]:
        return [e for e in self.elements if e.family == FAMILY]


class _Contig:
    def __init__(self, rng: random.Random):
        self.rng = rng
        self.parts: List[str] = []
        self.pos = 0

    def rand(self, n: int) -> str:
        return "".join(self.rng.choice("ACGT") for _ in range(n))

    def add(self, seq: str) -> int:
        """Append `seq`; return its 1-based start."""
        self.parts.append(seq)
        start = self.pos + 1
        self.pos += len(seq)
        return start

    def tsd(self) -> str:
        while True:
            t = self.rand(5)
            if len(set(t)) >= 2:
                return t


def _mutate(rng: random.Random, seq: str, rate: float, keep_ends: int = 2) -> str:
    out = list(seq)
    for i in range(keep_ends, len(out) - keep_ends):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def _patch(rng: random.Random, seq: str, n: int) -> str:
    out = list(seq)
    for i in rng.sample(range(len(out)), n):
        out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def build(root: Path, prefix: str = "syn_LTRs", seed: int = 7) -> Fixture:
    rng = random.Random(seed)
    g = _Contig(rng)
    ltr0 = "TG" + g.rand(LTR_LEN - 4) + "CA"
    internal0 = g.rand(INT_LEN)
    elements: List[Element] = []

    def copy_pair():
        base = _mutate(rng, ltr0, 0.02)
        return (_mutate(rng, base, 0.003), _mutate(rng, base, 0.003),
                _mutate(rng, internal0, 0.02, keep_ends=0))

    def place(kind, left, internal, right, strand, family=FAMILY, broken_tsd=False, k2p=0.01):
        """Insert one copy between TSDs. `left`/`right` are the biological 5'/3' LTRs."""
        t = g.tsd()
        g.add(t)
        body = left + internal + right if strand == "+" else rc(left + internal + right)
        start = g.add(body)
        end = start + len(body) - 1
        g.add(t if not broken_tsd else t[:2] + ("A" if t[2] != "A" else "C") + t[3:])
        forward_left = len(left) if strand == "+" else len(right)
        forward_right = len(right) if strand == "+" else len(left)
        e = Element(kind, family, strand, start, end, start, end,
                    start + forward_left - 1, end - forward_right + 1,
                    tsd="." if broken_tsd else t, k2p=k2p)
        elements.append(e)
        return e

    for i in range(12):
        g.add(g.rand(SPACER))
        a, b, internal = copy_pair()
        place("normal", a, internal, b, "+" if i % 3 else "-", k2p=0.004 + 0.001 * i)
    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    place("decayed", a, internal, b, "+", broken_tsd=True, k2p=0.02)

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("del_left", a[:50] + a[110:], internal, b, "+", k2p=0.03)
    e.start, e.r0, e.tsd = e.true_start + 50, e.r0 + 110, "."

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("patch_right", a, internal, b[:330] + _patch(rng, b[330:360], 9) + b[360:], "+",
              k2p=0.03)
    e.end, e.l1, e.tsd = e.true_end - 70, e.l1 - 70, "."

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("ins_left", a[:40] + g.rand(1500) + a[40:], internal, b, "+", k2p=0.03)
    e.start, e.r0, e.tsd = e.true_start + 40 + 1500, e.r0 + 40, "."

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("del_left_minus", a[:50] + a[110:], internal, b, "-", k2p=0.03)
    e.end, e.l1, e.tsd = e.true_end - 50, e.l1 - 110, "."

    g.add(g.rand(SPACER))
    host_ltr = "TG" + g.rand(296) + "CA"
    host_tsd = g.tsd()
    g.add(host_tsd)
    host_start = g.add(host_ltr + g.rand(1500))
    a, b, internal = copy_pair()
    inner = place("nested_del_left", a[:50] + a[110:], internal, b, "+", k2p=0.03)
    inner.start, inner.r0, inner.tsd = inner.true_start + 50, inner.r0 + 110, "."
    g.add(g.rand(1500) + _mutate(rng, host_ltr, 0.01))
    host_end = g.pos
    g.add(host_tsd)
    host = Element("host", "syn_fam00002", "+", host_start, host_end, host_start, host_end,
                   host_start + 299, host_end - 299, depth=1, tsd=host_tsd, k2p=0.01)
    elements.append(host)
    host.nest_status = f"nest-outer:{inner.key}"
    inner.nest_status = f"nest-inner:{host.key}"

    g.add(g.rand(SPACER))
    a, b, internal = copy_pair()
    e = place("conflict_left", a[:50] + a[110:], internal, b, "+", k2p=0.03)
    e.start, e.r0, e.tsd = e.true_start + 50, e.r0 + 110, "."
    xs, xe = e.true_start - 405, e.true_start + 20
    elements.append(Element("neighbour", "syn_fam00003", "+", xs, xe, xs, xe, xs + 99, xe - 99))
    g.add(g.rand(SPACER))

    contig = "".join(g.parts)
    root.mkdir(parents=True, exist_ok=True)
    genome = root / "genome.fa"
    with open(genome, "w") as fh:
        fh.write(f">{CHROM}\n")
        for i in range(0, len(contig), 60):
            fh.write(contig[i:i + 60] + "\n")

    header = "#" + "\t".join(COLUMNS + TAIL)
    for depth in (0, 1):
        rows = [e for e in elements if e.depth == depth]
        for variant in ("_clean_ltr", "_ltr"):
            base = root / f"{prefix}_depth{depth}{variant}"
            with open(f"{base}.tsv", "w") as fh:
                fh.write(header + "\n")
                for e in rows:
                    fh.write("\t".join(_row(e)) + "\n")
            with open(f"{base}.fa", "w") as fh:
                for e in rows:
                    seq = contig[e.start - 1:e.end]
                    if e.kind == "host":
                        a0, a1 = inner.start - e.start, inner.end - e.start + 1
                        seq = seq[:a0] + "N" * (a1 - a0) + seq[a1:]
                    if e.strand == "-":
                        seq = revcomp_record(seq)
                    fh.write(f">{e.name}\n")
                    for i in range(0, len(seq), 60):
                        fh.write(seq[i:i + 60] + "\n")
    return Fixture(genome, root, prefix, elements, ltr0)


def _row(e: Element) -> List[str]:
    n = e.end - e.start + 1
    l5, l3 = e.l1 - e.start + 1, e.end - e.r0 + 1
    vals = {
        "seq_id": e.name, "seq_len": str(n), "status": "pass",
        "ltr5_start": "1", "ltr5_end": str(l5), "ltr3_start": str(e.r0 - e.start + 1),
        "ltr3_end": str(n), "ltr5_len": str(l5), "ltr3_len": str(l3),
        "flank5_len": "0", "flank3_len": "0", "aln_len": str(max(l5, l3)),
        "n_sites": str(min(l5, l3)), "n_ts": "1", "n_tv": "1", "n_gapcols": "0",
        "identity": "0.99", "p_dist": "0.01", "k2p": f"{e.k2p:g}", "k2p_se": "0.001",
        "bitscore": "500", "flank_margin_bits": "NA", "cigar": ".", "motif": "tg...ca",
        "k2p_time": str(round(e.k2p / 6e-8)), "orientation": "-" if e.strand == "-" else "+",
        "tsd": e.tsd, "tsd_offset": "NA" if e.tsd == "." else "0,0", "tsd_input": e.tsd,
        "strand": e.strand, "family": e.family, "domains": ".", "nest_status": e.nest_status,
    }
    return [vals[c] for c in COLUMNS + TAIL]
```

- [ ] **Step 2: Write the fixture's own tests**

Create `tests/test_reboundary_fixtures.py`:

```python
"""The synthetic run is only useful if its truth is what it claims."""

from __future__ import annotations

import pytest

from ltrquest.detect import revcomp as revcomp_record

from reboundary_fixtures import CHROM, FAMILY, build


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build(tmp_path_factory.mktemp("syn"))


def contig(fx):
    lines = fx.genome.read_text().splitlines()
    assert lines[0] == f">{CHROM}"
    return "".join(lines[1:])


def test_every_family_copy_but_the_decayed_one_has_a_tsd_at_its_true_ends(fx):
    seq = contig(fx)
    for e in fx.family_members():
        left, right = seq[e.true_start - 6:e.true_start - 1], seq[e.true_end:e.true_end + 5]
        assert (left == right) == (e.kind != "decayed"), e.kind


def test_calls_sit_inside_their_true_spans_and_only_truncated_ones_differ(fx):
    for e in fx.family_members():
        assert e.true_start <= e.start and e.end <= e.true_end
        truncated = e.kind not in ("normal", "decayed")
        assert ((e.start, e.end) != (e.true_start, e.true_end)) == truncated, e.kind


def test_the_intact_ltrs_start_tg_and_end_ca(fx):
    seq = contig(fx)
    e = [x for x in fx.family_members() if x.kind == "normal" and x.strand == "+"][0]
    assert seq[e.start - 1:e.start + 1] == "TG" and seq[e.end - 2:e.end] == "CA"
    assert e.l1 - e.start + 1 == 400 and e.end - e.r0 + 1 == 400


def test_fasta_records_are_the_called_spans_stored_as_the_tables_say(fx):
    seq = contig(fx)
    records = {}
    for path in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.fa"):
        name = None
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                name = line[1:]
                records[name] = ""
            else:
                records[name] += line
    for e in fx.elements:
        fwd = seq[e.start - 1:e.end]
        stored = records[e.name]
        if e.kind == "host":
            inner = fx.kind("nested_del_left")
            assert stored[inner.start - e.start:inner.end - e.start + 1] == "N" * (
                inner.end - inner.start + 1)
            continue
        assert stored == (revcomp_record(fwd) if e.strand == "-" else fwd), e.kind


def test_the_tables_carry_the_33_columns_and_the_nesting(fx):
    rows = (fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv").read_text().splitlines()
    assert rows[0].startswith("#seq_id\tseq_len\tstatus") and rows[0].endswith("nest_status")
    assert len(rows[0].split("\t")) == 33
    host = fx.kind("host")
    inner = fx.kind("nested_del_left")
    assert host.depth == 1 and inner.nest_status == f"nest-inner:{host.key}"
    assert all(e.family == FAMILY for e in fx.family_members())
```

- [ ] **Step 3: Run the tests**

Run: `PYT tests/test_reboundary_fixtures.py -q`
Expected: 5 passed. (`tests/` is on `sys.path` because pytest's rootdir inserts it; `from reboundary_fixtures import ...` works the same way `test_recover_strand.py` imports its helpers. If the import fails, add `sys.path.insert(0, str(Path(__file__).parent))` at the top of the test module.)

- [ ] **Step 4: Commit**

```bash
git add tests/reboundary_fixtures.py tests/test_reboundary_fixtures.py
git commit -m "tests: synthetic run with truncated LTR-RT calls for re-boundarying

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: `ltr_model.py` — members, genome access, alignment primitives, consensus model, QC maths

**Files:**
- Create: `src/ltrquest/ltr_model.py`
- Modify: `tests/reboundary_fixtures.py` (add `members()` helper)
- Test: `tests/test_ltr_model.py`

**Interfaces:**
- Consumes: nothing from earlier tasks except the fixture.
- Produces (exact names; later tasks rely on them):
  - constants `MATCH=2, MISMATCH=-3, GAP_OPEN=5, GAP_EXTEND=2, TERM=30, PAD=200, MIN_REFS=10, KMER=11`
  - `matrix()`, `rc(seq) -> str`
  - `Member(prefix, name, chrom, start, end, l1, r0, strand, orientation, family, depth, k2p, tsd, nest_status=".")` with properties `key`, `uid`, `len_left`, `len_right`, `stranded`, `has_tsd`
  - `Model(model_id, family, seq, modal_len, n_refs)` with field `kmers`
  - `Genomes(paths: Dict[str, str])` with `.fetch(prefix, chrom, lo, hi) -> (seq, first)` and `.length(prefix, chrom)`
  - `canonical_kmers(seq, k=KMER)`, `jaccard(a, b)`
  - `Glocal(start, end, id_first, id_last, id_all, score)`; `glocal(model, window, term=TERM) -> Glocal`
  - `oriented_ltrs(m, g, pad=0) -> (five, three)`
  - `modal_length(lengths) -> float`
  - `select_references(members, strategy, seed_key, exclude=frozenset(), n_young=40, n_random=40) -> (refs, modal)`
  - `mafft_align(seqs, mafft="mafft") -> List[str]`
  - `consensus_from_alignment(aln, sides) -> Optional[str]`
  - `polish_ends(seq, refs, g, k=8, thr=0.7, rounds=3) -> (seq, (ds, de))`
  - `consensus_model(family, refs, g, modal, mafft="mafft", model_id=None) -> Optional[Model]`
  - `ratio_ok(model, max_ratio) -> bool`, `binom_sf(k, n, p) -> float`, `tsd_enrichment(hits, n, p0, min_n, alpha=0.01) -> str`

- [ ] **Step 1: Add the fixture helper**

Append to `tests/reboundary_fixtures.py`:

```python
def members(fx: Fixture):
    """The fixture's elements as ltr_model.Member objects, as reboundary_io would load them."""
    from ltrquest.ltr_model import Member
    return [Member(prefix=fx.prefix, name=e.name, chrom=CHROM, start=e.start, end=e.end,
                   l1=e.l1, r0=e.r0, strand=e.strand, orientation="-" if e.strand == "-" else "+",
                   family=e.family, depth=e.depth, k2p=e.k2p, tsd=e.tsd,
                   nest_status=e.nest_status)
            for e in fx.elements]
```

- [ ] **Step 2: Write the failing tests**

Create `tests/test_ltr_model.py`:

```python
"""Family LTR models: termini from 5'/3' agreement, never from a TSD or a motif."""

from __future__ import annotations

import random

import pytest

from ltrquest import ltr_model as lm

from reboundary_fixtures import FAMILY, build, members


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build(tmp_path_factory.mktemp("syn_model"))


@pytest.fixture(scope="module")
def g(fx):
    return lm.Genomes({fx.prefix: str(fx.genome)})


@pytest.fixture(scope="module")
def fam(fx):
    return [m for m in members(fx) if m.family == FAMILY]


def test_genomes_fetch_clips_to_the_contig(fx, g):
    seq, first = g.fetch(fx.prefix, "chrS", -5, 10)
    assert first == 1 and len(seq) == 10
    n = g.length(fx.prefix, "chrS")
    assert g.fetch(fx.prefix, "chrS", n - 2, n + 50)[0] == g.fetch(fx.prefix, "chrS", n - 2, n)[0]


def test_genomes_survive_pickling(fx, g):
    import pickle
    g2 = pickle.loads(pickle.dumps(g))
    assert g2.fetch(fx.prefix, "chrS", 1, 20) == g.fetch(fx.prefix, "chrS", 1, 20)


def test_member_geometry(fam):
    m = [x for x in fam if x.name.startswith("chrS") and x.len_left == 400][0]
    assert m.len_right == 400 and m.key == m.name.split("#")[0] and m.uid.endswith(m.key)


def test_modal_length_picks_the_dense_bin():
    assert abs(lm.modal_length([400] * 20 + [290, 330, 360, 1900]) - 400) <= 5


def test_select_references_takes_the_modal_class_and_honours_exclusion(fam):
    refs, modal = lm.select_references(fam, "modal", FAMILY)
    assert abs(modal - 400) <= 5
    assert all(abs(r.len_left - modal) <= 0.15 * modal for r in refs)
    assert len(refs) >= lm.MIN_REFS
    keys = {r.key for r in refs}
    refs2, _ = lm.select_references(fam, "modal", FAMILY, exclude=frozenset(list(keys)[:2]))
    assert not (set(list(keys)[:2]) & {r.key for r in refs2})


def test_select_references_by_tsd_uses_only_tsd_bearing_copies(fam):
    refs, _ = lm.select_references(fam, "tsd", FAMILY)
    assert refs and all(r.has_tsd for r in refs)


def test_consensus_from_alignment_finds_the_ltr_between_disagreeing_flanks():
    rng = random.Random(3)
    ltr = "TG" + "".join(rng.choice("ACGT") for _ in range(196)) + "CA"
    rows, sides = [], []
    for _ in range(12):   # pre-aligned rows, no MAFFT: 5' rows = flank A + LTR + internal G,
        rows.append("A" * 60 + ltr + "G" * 60)   # 3' rows = internal C + LTR + flank T, so the
        sides.append(5)                           # two sides disagree everywhere but the LTR
        rows.append("C" * 60 + ltr + "T" * 60)
        sides.append(3)
    assert lm.consensus_from_alignment(rows, sides) == ltr


def test_glocal_places_a_model_exactly():
    rng = random.Random(4)
    rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    model = rand(300)
    window = rand(120) + model + rand(80)
    hit = lm.glocal(model, window)
    assert (hit.start, hit.end) == (120, 419)
    assert hit.id_first == hit.id_last == hit.id_all == 1.0


def test_consensus_model_recovers_the_planted_ltr(fx, g, fam, mafft):
    refs, modal = lm.select_references(fam, "modal", FAMILY)
    model = lm.consensus_model(FAMILY, refs, g, modal, mafft)
    assert model is not None and abs(len(model.seq) - 400) <= 2
    assert model.seq[:20] == fx.ltr[:20] and model.seq[-20:] == fx.ltr[-20:]
    assert lm.glocal(model.seq, fx.ltr).id_all >= 0.99


def test_polish_ends_restores_clipped_termini(fx, g, fam):
    refs, _ = lm.select_references(fam, "modal", FAMILY)
    seq, (ds, de) = lm.polish_ends(fx.ltr[3:-3], refs, g)
    assert (ds, de) == (3, 3) and seq[:10] == fx.ltr[:10] and seq[-10:] == fx.ltr[-10:]


def test_polish_ends_trims_overhangs(fx, g, fam):
    refs, _ = lm.select_references(fam, "modal", FAMILY)
    seq, (ds, de) = lm.polish_ends("ACG" + fx.ltr + "TTA", refs, g)
    assert (ds, de) == (-3, -3) and seq[:10] == fx.ltr[:10] and seq[-10:] == fx.ltr[-10:]


def test_ratio_ok():
    m = lm.Model("f:consensus", "f", "A" * 459, 400.0, 20)
    assert lm.ratio_ok(m, 1.15) and not lm.ratio_ok(m, 1.10)


def test_binomial_tail():
    assert lm.binom_sf(0, 5, 0.3) == 1.0
    assert abs(lm.binom_sf(5, 5, 0.5) - 1 / 32) < 1e-12


def test_tsd_enrichment_has_an_explicit_untested_branch():
    assert lm.tsd_enrichment(3, 3, 0.01, min_n=5) == "untested"
    assert lm.tsd_enrichment(6, 8, 0.01, min_n=5) == "pass"
    assert lm.tsd_enrichment(0, 8, 0.01, min_n=5) == "fail"
```

- [ ] **Step 3: Run to verify failure**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT tests/test_ltr_model.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'ltrquest.ltr_model'`.

- [ ] **Step 4: Implement `src/ltrquest/ltr_model.py`**

```python
"""Family LTR models for ltrquest.reboundary: what a family's LTR is, end to end.

A model's termini come only from sequence agreement -- where a family's 5' and 3'
LTR copies stop matching each other -- never from a target-site duplication or a
terminal motif, so the TSD stays an independent check on what reboundary does.

Coordinates are genome-forward, 1-based, inclusive. `Member.l1` and `Member.r0`
are the inner ends of the genomic-left and genomic-right LTR, which is how the
depth tables store them even for minus-strand elements.
"""

from __future__ import annotations

import math
import os
import random
import subprocess
import tempfile
import zlib
from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

import numpy as np
import parasail
import pyfaidx

MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2, -3, 5, 2
TERM = 30         # outer bases whose identity gates a moved end
PAD = 200         # flank kept around each reference LTR in the model alignment
MIN_REFS = 10     # fewer reference copies than this and no model is built
KMER = 11
_COMP = str.maketrans("ACGTN", "TGCAN")
_MATRIX = None


def matrix():
    """The parasail scoring matrix, built once per process (it does not pickle)."""
    global _MATRIX
    if _MATRIX is None:
        _MATRIX = parasail.matrix_create("ACGTN", MATCH, MISMATCH)
    return _MATRIX


def rc(seq: str) -> str:
    """Reverse complement of genomic sequence (ACGTN). Records use detect.revcomp."""
    return seq.translate(_COMP)[::-1]


@dataclass(frozen=True)
class Member:
    prefix: str
    name: str
    chrom: str
    start: int
    end: int
    l1: int
    r0: int
    strand: str
    orientation: str
    family: str
    depth: int
    k2p: Optional[float]
    tsd: str
    nest_status: str = "."

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def uid(self) -> str:
        """Unique across genomes: the same coordinates can occur in two assemblies."""
        return f"{self.prefix}\t{self.key}"

    @property
    def len_left(self) -> int:
        return self.l1 - self.start + 1

    @property
    def len_right(self) -> int:
        return self.end - self.r0 + 1

    @property
    def stranded(self) -> bool:
        return self.strand in ("+", "-")

    @property
    def has_tsd(self) -> bool:
        return self.tsd not in ("", ".", "NA")


def canonical_kmers(seq: str, k: int = KMER) -> FrozenSet[str]:
    """Strand-independent k-mer set: each k-mer or its reverse complement, whichever is smaller."""
    out = set()
    for i in range(len(seq) - k + 1):
        w = seq[i:i + k]
        if "N" in w:
            continue
        r = rc(w)
        out.add(w if w <= r else r)
    return frozenset(out)


def jaccard(a: FrozenSet[str], b: FrozenSet[str]) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


@dataclass(frozen=True)
class Model:
    model_id: str
    family: str
    seq: str            # biological 5'->3'
    modal_len: float    # the family's modal called LTR length
    n_refs: int
    kmers: FrozenSet[str] = field(default=frozenset(), repr=False, compare=False)

    def __post_init__(self):
        if not self.kmers:
            object.__setattr__(self, "kmers", canonical_kmers(self.seq))


class Genomes:
    """Random access to each genome by prefix. Nothing is ever read whole."""

    def __init__(self, paths: Dict[str, str]):
        self.paths = dict(paths)
        self._open: Dict[str, pyfaidx.Fasta] = {}

    def __getstate__(self):
        return {"paths": self.paths}

    def __setstate__(self, state):
        self.paths = state["paths"]
        self._open = {}

    def _fasta(self, prefix: str) -> pyfaidx.Fasta:
        fa = self._open.get(prefix)
        if fa is None:
            fa = pyfaidx.Fasta(self.paths[prefix], as_raw=True, sequence_always_upper=True)
            self._open[prefix] = fa
        return fa

    def length(self, prefix: str, chrom: str) -> int:
        return len(self._fasta(prefix)[chrom])

    def fetch(self, prefix: str, chrom: str, lo: int, hi: int) -> Tuple[str, int]:
        """Bases lo..hi (1-based, inclusive), clipped to the sequence. Returns (seq, first)."""
        n = self.length(prefix, chrom)
        a, b = max(1, lo), min(n, hi)
        if a > b:
            return "", a
        return self._fasta(prefix)[chrom][a - 1:b], a


@dataclass(frozen=True)
class Glocal:
    start: int        # 0-based first window base the model covers
    end: int          # 0-based last window base
    id_first: float   # identity over the model's first `term` bases
    id_last: float
    id_all: float
    score: int
    head: Tuple[int, ...] = ()   # window index of model bases 0..edge-1 (-1: aligned to a gap)
    tail: Tuple[int, ...] = ()   # same for the last bases, outermost first


def _identity(q: str, t: str) -> float:
    n = sum(1 for a in q if a != "-")
    m = sum(1 for a, b in zip(q, t) if a == b and a != "-")
    return m / n if n else 0.0


def glocal(model: str, window: str, term: int = TERM, edge: int = 8) -> Glocal:
    """Align `model` end to end inside `window`; only the window's ends are free."""
    r = parasail.sg_dx_trace_scan_sat(model, window, GAP_OPEN, GAP_EXTEND, matrix())
    q, t = r.traceback.query, r.traceback.ref
    i0 = 0
    while i0 < len(q) and q[i0] == "-":
        i0 += 1
    i1 = len(q) - 1
    while i1 >= 0 and q[i1] == "-":
        i1 -= 1
    start = i0 - t[:i0].count("-")
    end = start + sum(1 for c in t[i0:i1 + 1] if c != "-") - 1
    q, t = q[i0:i1 + 1], t[i0:i1 + 1]
    k = j = 0
    while j < len(q) and k < term:
        k += q[j] != "-"
        j += 1
    k, j2 = 0, len(q)
    while j2 > 0 and k < term:
        j2 -= 1
        k += q[j2] != "-"
    qmap: List[int] = []
    wi = start
    for a, b in zip(q, t):
        if a != "-":
            qmap.append(wi if b != "-" else -1)
        if b != "-":
            wi += 1
    return Glocal(start, end, _identity(q[:j], t[:j]), _identity(q[j2:], t[j2:]),
                  _identity(q, t), r.score, tuple(qmap[:edge]), tuple(reversed(qmap[-edge:])))


def oriented_ltrs(m: Member, g: Genomes, pad: int = 0) -> Tuple[str, str]:
    """(5' LTR, 3' LTR), biological 5'->3', each padded by `pad` bases of context."""
    left, _ = g.fetch(m.prefix, m.chrom, m.start - pad, m.l1 + pad)
    right, _ = g.fetch(m.prefix, m.chrom, m.r0 - pad, m.end + pad)
    return (left, right) if m.strand != "-" else (rc(right), rc(left))


def modal_length(lengths: Sequence[int]) -> float:
    """Mode of LTR length: densest 5%-wide log bin (smoothed by its neighbours), its median."""
    bins = Counter(int(math.log(x) / math.log(1.05)) for x in lengths if x > 0)
    top = max(bins, key=lambda b: bins[b] + 0.5 * (bins.get(b - 1, 0) + bins.get(b + 1, 0)))
    inbin = [x for x in lengths if x > 0 and abs(int(math.log(x) / math.log(1.05)) - top) <= 1]
    return float(np.median(inbin))


def select_references(members: Sequence[Member], strategy: str, seed_key: str,
                      exclude: FrozenSet[str] = frozenset(), n_young: int = 40,
                      n_random: int = 40) -> Tuple[List[Member], Optional[float]]:
    """Reference copies for a model, and the family's modal called LTR length.

    `modal`: both called LTRs within 15% of the modal length. `tsd`: copies whose
    called ends carry a Kmer2LTR TSD (independent evidence the call is right),
    falling back to `modal` below MIN_REFS. The youngest first, then a seeded
    random draw of the rest, so the choice does not depend on table order.
    """
    stranded = [m for m in members if m.stranded and m.key not in exclude]
    if len(stranded) < MIN_REFS:
        return [], None
    modal = modal_length([m.len_left for m in stranded] + [m.len_right for m in stranded])
    pool: List[Member] = []
    if strategy == "tsd":
        pool = [m for m in stranded if m.has_tsd]
    if len(pool) < MIN_REFS:
        pool = [m for m in stranded if abs(m.len_left - modal) <= 0.15 * modal
                and abs(m.len_right - modal) <= 0.15 * modal]
    if len(pool) < MIN_REFS:
        return [], modal
    age = lambda m: (math.inf if m.k2p is None else m.k2p, m.key)  # noqa: E731
    young = sorted(pool, key=age)[:n_young]
    taken = {m.key for m in young}
    rest = sorted((m for m in pool if m.key not in taken), key=lambda m: m.key)
    rng = random.Random(zlib.crc32(seed_key.encode()))
    return young + rng.sample(rest, min(n_random, len(rest))), modal


def mafft_align(seqs: Sequence[str], mafft: str = "mafft") -> List[str]:
    fd, path = tempfile.mkstemp(suffix=".fa")
    try:
        with os.fdopen(fd, "w") as fo:
            for i, s in enumerate(seqs):
                fo.write(f">{i}\n{s}\n")
        out = subprocess.run([mafft, "--auto", "--thread", "1", "--quiet", path],
                             capture_output=True, text=True, check=True).stdout
    finally:
        os.unlink(path)
    rows: Dict[int, List[str]] = {}
    cur = -1
    for line in out.splitlines():
        if line.startswith(">"):
            cur = int(line[1:])
            rows[cur] = []
        elif cur >= 0:
            rows[cur].append(line.strip().upper())
    return ["".join(rows[i]) for i in range(len(seqs))]


def consensus_from_alignment(aln: Sequence[str], sides: Sequence[int]) -> Optional[str]:
    """The LTR in an alignment of padded 5'-LTR (side 5) and 3'-LTR (side 3) rows.

    A column is LTR when both sides agree on a base carried by >= 40% of each.
    Beyond the outer end the 5' rows are flank; beyond the inner end the 3' rows
    are; either way the sides disagree. The span runs from the first to the last
    sustained block (>= 15 of 21 columns) so dips inside the LTR are ignored.
    """
    A = np.frombuffer("".join(aln).encode(), dtype="S1").reshape(len(aln), -1)
    side = np.asarray(sides)
    bases = np.array([b"A", b"C", b"G", b"T"])
    c5 = np.stack([(A[side == 5] == b).sum(0) for b in bases])
    c3 = np.stack([(A[side == 3] == b).sum(0) for b in bases])
    f5 = c5.max(0) / max(1, int((side == 5).sum()))
    f3 = c3.max(0) / max(1, int((side == 3).sum()))
    agree = (c5.argmax(0) == c3.argmax(0)) & (f5 >= 0.4) & (f3 >= 0.4)
    smooth = np.convolve(agree.astype(float), np.ones(21) / 21, mode="same")
    hit = np.flatnonzero(smooth >= 15 / 21)
    if len(hit) == 0:
        return None
    lo, hi, n = int(hit[0]), int(hit[-1]), len(agree)
    if hi - lo < 50:
        return None
    while lo > 0 and agree[lo - 1]:
        lo -= 1
    while hi < n - 1 and agree[hi + 1]:
        hi += 1
    while lo < hi and not agree[lo]:
        lo += 1
    while hi > lo and not agree[hi]:
        hi -= 1
    tot = c5 + c3
    occ = tot.sum(0) / len(aln)
    return "".join(bases[tot[:, j].argmax()].decode()
                   for j in range(lo, hi + 1) if occ[j] >= 0.5)


def polish_ends(seq: str, refs: Sequence[Member], g: Genomes, k: int = 8, thr: float = 0.7,
                rounds: int = 3) -> Tuple[str, Tuple[int, int]]:
    """Move the model's ends to where each reference's two LTRs stop agreeing.

    That is the definition of an LTR end. The model is placed in both LTRs of
    every reference and the base pairs at and just beyond each end are compared:
    inside the LTR they agree (1 - divergence); beyond it one side is flank and
    the other internal region. Extend while >= `thr` of the references agree,
    trim while fewer do. A model base aligned to a gap counts as disagreement,
    so bases no reference carries are trimmed. Returns (seq, (bases added at the
    start, bases added at the end)); negative = trimmed.
    """
    total_s = total_e = 0
    for _ in range(rounds):
        placed = []
        for m in refs:
            s5, s3 = oriented_ltrs(m, g, PAD)
            g5, g3 = glocal(seq, s5, edge=k), glocal(seq, s3, edge=k)
            if min(g5.id_all, g3.id_all) >= 0.6:
                placed.append((s5, s3, g5, g3))
        if len(placed) < MIN_REFS:
            break

        def vote(pairs) -> Tuple[float, str]:
            votes: Counter = Counter()
            hits = 0
            for pair in pairs:
                if pair is not None and pair[0] == pair[1] and pair[0] != "N":
                    hits += 1
                    votes[pair[0]] += 1
            return hits / len(placed), (votes.most_common(1)[0][0] if votes else "N")

        def inside(i: int, at_start: bool) -> Tuple[float, str]:
            """Agreement at model base i counted from that end (0 = the outermost)."""
            pairs = []
            for s5, s3, g5, g3 in placed:
                a = (g5.head if at_start else g5.tail)[i]
                b = (g3.head if at_start else g3.tail)[i]
                pairs.append((s5[a], s3[b]) if a >= 0 and b >= 0 else None)
            return vote(pairs)

        def outside(j: int, at_start: bool) -> Tuple[float, str]:
            """Agreement j bases beyond that end of the model (1 = adjacent)."""
            pairs = []
            for s5, s3, g5, g3 in placed:
                i5 = g5.start - j if at_start else g5.end + j
                i3 = g3.start - j if at_start else g3.end + j
                ok = 0 <= i5 < len(s5) and 0 <= i3 < len(s3)
                pairs.append((s5[i5], s3[i3]) if ok else None)
            return vote(pairs)

        def trim_to(at_start: bool) -> int:
            """Model bases to drop from this end: up to the first 3 agreeing in a row.

            A lone agreeing base is not trusted: the aligner parks an overhanging
            base on any matching base nearby, in both LTRs at once.
            """
            for t in range(k - 2):
                if all(inside(t + d, at_start)[0] >= thr for d in range(3)):
                    return t
            return 0

        ds = 0
        while ds < k and outside(ds + 1, True)[0] >= thr:
            ds += 1
        if ds == 0:
            ds = -trim_to(True)
        de = 0
        while de < k and outside(de + 1, False)[0] >= thr:
            de += 1
        if de == 0:
            de = -trim_to(False)
        if ds == 0 and de == 0:
            break
        pre = "".join(outside(j, True)[1] for j in range(ds, 0, -1)) if ds > 0 else ""
        post = "".join(outside(j, False)[1] for j in range(1, de + 1)) if de > 0 else ""
        seq = pre + seq[(-ds if ds < 0 else 0):(len(seq) + de if de < 0 else len(seq))] + post
        total_s += ds
        total_e += de
    return seq, (total_s, total_e)


def consensus_model(family: str, refs: Sequence[Member], g: Genomes, modal: Optional[float],
                    mafft: str = "mafft", model_id: Optional[str] = None) -> Optional[Model]:
    """One model from `refs`: MAFFT their padded 5'/3' LTRs, take the agreed span, polish."""
    if len(refs) < MIN_REFS or modal is None:
        return None
    seqs: List[str] = []
    sides: List[int] = []
    for m in refs:
        s5, s3 = oriented_ltrs(m, g, PAD)
        seqs += [s5, s3]
        sides += [5, 3]
    seq = consensus_from_alignment(mafft_align(seqs, mafft), sides)
    if seq is None:
        return None
    seq, _ = polish_ends(seq, refs, g)
    return Model(model_id or f"{family}:consensus", family, seq, modal, len(refs))


def ratio_ok(model: Model, max_ratio: float) -> bool:
    """A model much longer than the family's usual called LTR is not an LTR model."""
    return len(model.seq) <= max_ratio * model.modal_len


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p)."""
    if k <= 0:
        return 1.0
    return float(sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1)))


def tsd_enrichment(hits: int, n: int, p0: float, min_n: int, alpha: float = 0.01) -> str:
    """Family QC on TSDs at proposed ends: 'pass', 'fail', or 'untested' (n < min_n).

    `p0` is the displaced-flank null rate, floored at 0.5% so one lucky family
    cannot pass on a null of zero.
    """
    if n < min_n:
        return "untested"
    return "pass" if binom_sf(hits, n, max(p0, 0.005)) < alpha else "fail"
```

- [ ] **Step 5: Run the tests**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT tests/test_ltr_model.py -q`
Expected: 14 passed. If `test_polish_ends_restores_clipped_termini` reports `(ds, de) == (3, 2)` or similar off-by-one, inspect `agree()` offsets against the spike's `polish_ends` (`/data2/chris/poa_LTR/ltrquest_run/reboundary_spike/reboundary_pilot.py`) — the offsets there were validated on real data.

- [ ] **Step 6: Commit**

```bash
git add src/ltrquest/ltr_model.py tests/test_ltr_model.py tests/reboundary_fixtures.py
git commit -m "ltr_model: family LTR models from 5'/3' agreement, reference choice, QC maths

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: `ltr_place.py` — proposals from a model: core, terminal anchor, gates, obstacles, TSD probe

**Files:**
- Create: `src/ltrquest/ltr_place.py`
- Test: `tests/test_ltr_place.py`

**Interfaces:**
- Consumes: `ltr_model` names from Task 3; `kmer2ltr.Api` from Task 1 (for `tsd_at`).
- Produces:
  - `Params(min_identity=0.8, min_whole_identity=0.6, min_ext=5, anchor=True, anchor_len=30, max_indel=5000, n_models=3)`
  - `Core(left, l1, r0, right, id_left, id_right, whole, score)`, `Hit(q_start, q_end, w_start, w_end, identity, score)`
  - `Proposal(member, model_id, orient, left, right, left_src, right_src, id_left, id_right, whole, raw_left, raw_right, inner_left, inner_right, gate_ok)` with properties `ext_left`, `ext_right`, `ext5`, `ext3`, `raw_ext5`, `raw_ext3`
  - `local(query, window) -> Optional[Hit]`
  - `core(m, fwd, g) -> Optional[Core]`
  - `anchor_left(m, fwd, g, p) -> Optional[(pos, identity)]`, `anchor_right(...)`
  - `candidate_models(m, models, g, n) -> List[Model]`
  - `propose(m, models, g, p) -> Optional[Proposal]` (None = the model does not move any end outward)
  - `obstacle(m, prop, g) -> (left_label, right_label)` — forward frame; `"."` for an end that did not move
  - `tsd_at(api, g, m, left, right, shift=0) -> str` (TSD motif, `"."` none, `"NA"` unreadable)

- [ ] **Step 1: Write the failing tests**

Create `tests/test_ltr_place.py`:

```python
"""Placing a family model on one element: extend only, gates, anchor, obstacles."""

from __future__ import annotations

import dataclasses

import pytest

from ltrquest import ltr_model as lm
from ltrquest import ltr_place as lp

from reboundary_fixtures import FAMILY, build, members


@pytest.fixture(scope="module")
def fx(tmp_path_factory):
    return build(tmp_path_factory.mktemp("syn_place"))


@pytest.fixture(scope="module")
def g(fx):
    return lm.Genomes({fx.prefix: str(fx.genome)})


@pytest.fixture(scope="module")
def truth_model(fx):
    """The planted LTR itself: placement tests should not depend on MAFFT."""
    return lm.Model(f"{FAMILY}:truth", FAMILY, fx.ltr, 400.0, 13)


def member(fx, kind):
    e = fx.kind(kind)
    return [m for m in members(fx) if m.name == e.name][0], e


def test_local_reports_query_and_window_coordinates():
    q = "ACGTACGTTTGACCAGTTAG"
    w = "GGGGG" + q[2:] + "CCCC"
    h = lp.local(q, w)
    assert (h.q_start, h.q_end, h.w_start, h.w_end) == (2, len(q) - 1, 5, 5 + len(q) - 3)
    assert h.identity == 1.0


def test_intact_copy_is_not_a_candidate(fx, g, truth_model):
    m = [x for x in members(fx) if x.name == fx.elements[0].name][0]
    assert lp.propose(m, [truth_model], g, lp.Params()) is None


def test_deletion_near_the_left_end_is_recovered_by_the_core(fx, g, truth_model):
    m, e = member(fx, "del_left")
    p = lp.propose(m, [truth_model], g, lp.Params())
    assert p.gate_ok and (p.left, p.right) == (e.true_start, e.end)
    assert p.left_src == "core" and p.right_src == "." and p.ext5 == 50 and p.ext3 == 0


def test_mutation_patch_near_the_right_end_is_recovered(fx, g, truth_model):
    m, e = member(fx, "patch_right")
    p = lp.propose(m, [truth_model], g, lp.Params())
    assert p.gate_ok and (p.left, p.right) == (e.start, e.true_end) and p.ext3 == 70


def test_minus_strand_copy_extends_at_its_biological_five_prime_end(fx, g, truth_model):
    m, e = member(fx, "del_left_minus")
    p = lp.propose(m, [truth_model], g, lp.Params())
    assert p.orient == "-" and (p.left, p.right) == (e.start, e.true_end)
    assert p.ext5 == 50 and p.ext3 == 0


def test_large_insertion_is_crossed_by_the_terminal_anchor(fx, g, truth_model):
    m, e = member(fx, "ins_left")
    p = lp.propose(m, [truth_model], g, lp.Params())
    assert p.gate_ok and p.left == e.true_start and p.left_src == "anchor"


def test_without_the_anchor_the_insertion_fails_its_gate(fx, g, truth_model):
    m, e = member(fx, "ins_left")
    p = lp.propose(m, [truth_model], g, lp.Params(anchor=False))
    assert p is None or (not p.gate_ok and p.left == m.start)


def test_never_trims_an_over_extended_call(fx, g, truth_model):
    base = [x for x in members(fx) if x.family == FAMILY and x.len_left == 400][0]
    wide = dataclasses.replace(base, start=base.start - 30, end=base.end + 30)
    assert lp.propose(wide, [truth_model], g, lp.Params()) is None


def test_identity_gate_blocks_a_weak_end(fx, g, truth_model):
    m, e = member(fx, "del_left")
    p = lp.propose(m, [truth_model], g, lp.Params(min_identity=1.01))
    assert p is not None and not p.gate_ok and p.left == m.start and p.raw_left == e.true_start


def test_obstacle_names_the_deletion(fx, g, truth_model):
    m, _ = member(fx, "del_left")
    p = lp.propose(m, [truth_model], g, lp.Params())
    left, right = lp.obstacle(m, p, g)
    assert left.startswith("gap:") and int(left.split(":")[1]) >= 50 and right == "."


def test_tsd_probe_reads_the_planted_duplication_and_not_the_null(fx, g, k2l_api):
    m, e = member(fx, "del_left")
    found = lp.tsd_at(k2l_api, g, m, e.true_start, e.true_end)
    assert found not in (".", "NA") and len(found) >= 5
    assert lp.tsd_at(k2l_api, g, m, m.start, m.end) == "."
```

- [ ] **Step 2: Run to verify failure**

Run: `PYT tests/test_ltr_place.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'ltrquest.ltr_place'`.

- [ ] **Step 3: Implement `src/ltrquest/ltr_place.py`**

```python
"""Placing a family LTR model on one element: where would its LTRs really end?

The model is aligned end to end in a window around each called LTR (the core).
When the core's outer bases do not match -- a large insertion sits between the
true end and the rest of the LTR -- the model's outermost `anchor_len` bases are
searched for locally up to `max_indel` beyond the call (the anchor).

Only outward moves are ever proposed. An end that would move inward stays where
the caller put it: in the feasibility spike, trimming was wrong most of the time.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import parasail

from .ltr_model import (GAP_EXTEND, GAP_OPEN, Genomes, Member, Model, canonical_kmers, glocal,
                        jaccard, matrix, rc)


@dataclass(frozen=True)
class Params:
    min_identity: float = 0.8       # outer-TERM identity a moved end must reach
    min_whole_identity: float = 0.6
    min_ext: int = 5                # an end must move at least this far to count
    anchor: bool = True
    anchor_len: int = 30
    max_indel: int = 5000
    n_models: int = 3               # models tried per element when a family has several


@dataclass(frozen=True)
class Core:
    left: int
    l1: int
    r0: int
    right: int
    id_left: float
    id_right: float
    whole: float
    score: int


@dataclass(frozen=True)
class Hit:
    q_start: int
    q_end: int
    w_start: int
    w_end: int
    identity: float
    score: int


@dataclass(frozen=True)
class Proposal:
    member: Member
    model_id: str
    orient: str
    left: int
    right: int
    left_src: str       # core | anchor | '.' (did not move)
    right_src: str
    id_left: float
    id_right: float
    whole: float
    raw_left: int       # before the gates, for reporting
    raw_right: int
    inner_left: int     # the model's inner ends, for the obstacle report
    inner_right: int
    gate_ok: bool

    @property
    def ext_left(self) -> int:
        return self.member.start - self.left

    @property
    def ext_right(self) -> int:
        return self.right - self.member.end

    @property
    def ext5(self) -> int:
        return self.ext_left if self.orient == "+" else self.ext_right

    @property
    def ext3(self) -> int:
        return self.ext_right if self.orient == "+" else self.ext_left

    @property
    def raw_ext5(self) -> int:
        left, right = self.member.start - self.raw_left, self.raw_right - self.member.end
        return max(0, left if self.orient == "+" else right)

    @property
    def raw_ext3(self) -> int:
        left, right = self.member.start - self.raw_left, self.raw_right - self.member.end
        return max(0, right if self.orient == "+" else left)


def local(query: str, window: str) -> Optional[Hit]:
    """Best local alignment of `query` in `window` (0-based, inclusive coordinates)."""
    if not query or not window:
        return None
    r = parasail.sw_trace_scan_sat(query, window, GAP_OPEN, GAP_EXTEND, matrix())
    if r.score <= 0:
        return None
    q, t = r.traceback.query, r.traceback.ref
    q_len = sum(1 for c in q if c != "-")
    t_len = sum(1 for c in t if c != "-")
    matches = sum(1 for a, b in zip(q, t) if a == b and a != "-")
    return Hit(r.end_query - q_len + 1, r.end_query, r.end_ref - t_len + 1, r.end_ref,
               matches / max(1, q_len), r.score)


def core(m: Member, fwd: str, g: Genomes) -> Optional[Core]:
    """The model (forward frame) glocal-aligned around each called LTR."""
    lc = len(fwd)
    extra = max(100, lc - min(m.len_left, m.len_right) + 100)
    wl, wl0 = g.fetch(m.prefix, m.chrom, m.start - extra, m.l1 + extra)
    wr, wr0 = g.fetch(m.prefix, m.chrom, m.r0 - extra, m.end + extra)
    if len(wl) < lc // 2 or len(wr) < lc // 2:
        return None
    a, b = glocal(fwd, wl), glocal(fwd, wr)
    return Core(wl0 + a.start, wl0 + a.end, wr0 + b.start, wr0 + b.end, a.id_first, b.id_last,
                min(a.id_all, b.id_all), a.score + b.score)


def anchor_left(m: Member, fwd: str, g: Genomes, p: Params) -> Optional[Tuple[int, float]]:
    q = fwd[:p.anchor_len]
    w, w0 = g.fetch(m.prefix, m.chrom, m.start - p.max_indel - len(q), m.start + 200)
    h = local(q, w)
    if h is None or h.q_start > 2 or (h.q_end - h.q_start + 1) < 0.9 * len(q):
        return None
    return w0 + h.w_start - h.q_start, h.identity


def anchor_right(m: Member, fwd: str, g: Genomes, p: Params) -> Optional[Tuple[int, float]]:
    q = fwd[-p.anchor_len:]
    w, w0 = g.fetch(m.prefix, m.chrom, m.end - 200, m.end + p.max_indel + len(q))
    h = local(q, w)
    if h is None or h.q_end < len(q) - 3 or (h.q_end - h.q_start + 1) < 0.9 * len(q):
        return None
    return w0 + h.w_end + (len(q) - 1 - h.q_end), h.identity


def candidate_models(m: Member, models: Sequence[Model], g: Genomes, n: int) -> List[Model]:
    """With several models, the `n` whose k-mers best match the element's called left LTR."""
    if len(models) <= 1:
        return list(models)
    s, _ = g.fetch(m.prefix, m.chrom, m.start, m.l1)
    ks = canonical_kmers(s)
    return sorted(models, key=lambda mo: (-jaccard(ks, mo.kmers), mo.model_id))[:n]


def propose(m: Member, models: Sequence[Model], g: Genomes, p: Params) -> Optional[Proposal]:
    best = None
    for model in candidate_models(m, models, g, p.n_models):
        for o in ((m.strand,) if m.stranded else ("+", "-")):
            fwd = model.seq if o == "+" else rc(model.seq)
            c = core(m, fwd, g)
            if c is not None and (best is None or c.score > best[0].score):
                best = (c, model, o, fwd)
    if best is None:
        return None
    c, model, o, fwd = best

    left, left_id, left_src = c.left, c.id_left, "core"
    if p.anchor and c.id_left < p.min_identity:
        a = anchor_left(m, fwd, g, p)
        if a is not None and a[1] >= p.min_identity and a[0] <= c.left + 2:
            left, left_id, left_src = a[0], a[1], "anchor"
    right, right_id, right_src = c.right, c.id_right, "core"
    if p.anchor and c.id_right < p.min_identity:
        a = anchor_right(m, fwd, g, p)
        if a is not None and a[1] >= p.min_identity and a[0] >= c.right - 2:
            right, right_id, right_src = a[0], a[1], "anchor"

    moves_left = left <= m.start - p.min_ext
    moves_right = right >= m.end + p.min_ext
    if not (moves_left or moves_right):
        return None
    whole_ok = c.whole >= p.min_whole_identity
    gate_left = moves_left and left_id >= p.min_identity and whole_ok
    gate_right = moves_right and right_id >= p.min_identity and whole_ok
    return Proposal(
        member=m, model_id=model.model_id, orient=o,
        left=left if gate_left else m.start, right=right if gate_right else m.end,
        left_src=left_src if gate_left else ".", right_src=right_src if gate_right else ".",
        id_left=left_id, id_right=right_id, whole=c.whole, raw_left=left, raw_right=right,
        inner_left=c.l1, inner_right=c.r0, gate_ok=gate_left or gate_right)


def _gap_near(q: str, t: str, c: int) -> int:
    best, n = 0, len(q)
    for i in range(max(0, c - 15), min(n, c + 16)):
        for s in (q, t):
            if s[i] == "-":
                j = i
                while j < n and s[j] == "-":
                    j += 1
                k = i
                while k > 0 and s[k - 1] == "-":
                    k -= 1
                best = max(best, j - k)
    return best


def _mismatch(q: str, t: str, cols) -> float:
    m = u = 0
    for j in cols:
        if 0 <= j < len(q) and q[j] != "-" and t[j] != "-":
            u += 1
            m += q[j] != t[j]
            if u == 30:
                break
    return m / u if u else 0.0


def obstacle(m: Member, prop: Proposal, g: Genomes) -> Tuple[str, str]:
    """What sits between the element's two LTRs at each old outer end that moved.

    The new left LTR (outer end to the model's inner end) is aligned to the new
    right LTR; the old end is mapped onto that alignment. Reported, never used
    to decide anything: 'gap:<bp>' for an indel >= 5 bp within 15 columns,
    'mm:<rate>' for >= 20% mismatches in the 30 bp just outside, else 'none'.
    """
    a, a0 = g.fetch(m.prefix, m.chrom, prop.left, prop.inner_left)
    b, b0 = g.fetch(m.prefix, m.chrom, prop.inner_right, prop.right)
    if not a or not b:
        return ".", "."
    r = parasail.nw_trace_scan_sat(a, b, 10, 1, matrix())
    q, t = r.traceback.query, r.traceback.ref
    col_a, col_b, ia, ib = {}, {}, a0, b0
    for j, (x, y) in enumerate(zip(q, t)):
        if x != "-":
            col_a[ia] = j
            ia += 1
        if y != "-":
            col_b[ib] = j
            ib += 1

    def label(col: Optional[int], outward: int) -> str:
        if col is None:
            return "."
        gap = _gap_near(q, t, col)
        if gap >= 5:
            return f"gap:{gap}"
        cols = range(col - 1, -1, -1) if outward < 0 else range(col + 1, len(q))
        mm = _mismatch(q, t, cols)
        return f"mm:{mm:.2f}" if mm >= 0.2 else "none"

    left = label(col_a.get(m.start), -1) if prop.left < m.start else "."
    right = label(col_b.get(m.end), +1) if prop.right > m.end else "."
    return left, right


def tsd_at(api, g: Genomes, m: Member, left: int, right: int, shift: int = 0) -> str:
    """Kmer2LTR's TSD call for an element spanning left..right; `shift` displaces the right flank."""
    a, a0 = g.fetch(m.prefix, m.chrom, left - 10, left + 9)
    b, _ = g.fetch(m.prefix, m.chrom, right - 9 + shift, right + 10 + shift)
    if a0 != left - 10 or len(a) != 20 or len(b) != 20:
        return "NA"
    hit = api.find_tsd(a + b, 10, 30, api.TSD_K, api.TSD_SHIFTS)
    return hit[0] if hit else "."
```


- [ ] **Step 4: Run the tests**

Run: `PYT tests/test_ltr_place.py -q`
Expected: 11 passed.
If `test_deletion_near_the_left_end_is_recovered_by_the_core` fails with `left_src == "anchor"`, the core preferred aligning the model's first bases into flank over paying for the 60 bp gap: print `lp.core(m, fx.ltr, g)` and check `id_left`; the anchor is then doing its job and the assertion on `left_src` may be relaxed to `in ("core", "anchor")` — the position assertion must not be relaxed.

- [ ] **Step 5: Commit**

```bash
git add src/ltrquest/ltr_place.py tests/test_ltr_place.py
git commit -m "ltr_place: outward-only proposals from a family model, with terminal anchor

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---
### Task 5: `reboundary_io.py` — tables, records, nesting, conflicts, sidecar, atomic rewrite

**Files:**
- Create: `src/ltrquest/reboundary_io.py`
- Test: `tests/test_reboundary_io.py`

**Interfaces:**
- Consumes: `annotate.discover_depth_tables/element_key/header_names/read_table/target_mode`, `detect.revcomp`, `kmer2ltr.COLUMNS`, `ltr_model.Member`, `table.Columns/as_int/as_float`.
- Produces:
  - `SIDECAR_SUFFIX = "_reboundary.tsv"`, `SIDECAR_COLUMNS` (22 names, spec §5.6 order)
  - `CleanTable(path, fasta, depth, header, rows, cols)`; `load_clean_tables(indir, prefix) -> List[CleanTable]`
  - `members_from(prefix, tables) -> List[Member]`
  - `iter_fasta(path)`, `fetch_records(paths, names) -> Dict[name, seq]`
  - `forward(record, orientation)`, `stored(fwd, orientation)`, `sanitize(seq)`, `paint(record, orientation, rec_start, spans, letter)`
  - `rekey_nest(value, old2new)`, `hosts_of(nest_status) -> List[key]`
  - `SpanIndex(members).overlapping(prefix, chrom, lo, hi)`; `conflict(index, m, left, right) -> Optional[str]`; `mutual_conflicts(items: [(Member, left, right)]) -> Set[uid]`
  - `Accepted(member, new_name, fields, record, left, right)`
  - `Commit` with `.open(path)`, `.commit()`, `.abort()`; `atomic_write_text(path, text)`
  - `rewrite(tables, accepted: Dict[old_key, Accepted], letters, commit, wrap=60)`
  - `sidecar_text(rows) -> str`; `Rebound(old_name, new_name, ext5, ext3)`; `read_map(path) -> Dict[new_key, Rebound]`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reboundary_io.py`:

```python
"""Reading and rewriting a run's clean side: nothing raw, nothing half-written."""

from __future__ import annotations

import dataclasses
import os

import pytest

from ltrquest import reboundary_io as rio
from ltrquest.detect import revcomp as revcomp_record
from ltrquest.kmer2ltr import COLUMNS
from ltrquest.reconcile import IUPAC_DEPTH_SEQ

from reboundary_fixtures import build, members


@pytest.fixture
def fx(tmp_path):
    return build(tmp_path / "syn")


def test_load_clean_tables_finds_both_depths_and_their_fastas(fx):
    tables = rio.load_clean_tables(str(fx.indir), fx.prefix)
    assert [t.depth for t in tables] == [0, 1]
    assert all(os.path.isfile(t.fasta) for t in tables)


def test_load_refuses_a_table_that_was_never_annotated(fx):
    path = fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv"
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    cut = head.index("family")
    path.write_text("\n".join("\t".join(l.split("\t")[:cut] + l.split("\t")[cut + 1:])
                              for l in lines) + "\n")
    with pytest.raises(ValueError, match="family"):
        rio.load_clean_tables(str(fx.indir), fx.prefix)


def test_members_from_reads_called_coordinates(fx):
    got = {m.name: m for m in rio.members_from(fx.prefix, rio.load_clean_tables(str(fx.indir),
                                                                                 fx.prefix))}
    e = fx.kind("del_left")
    m = got[e.name]
    assert (m.start, m.end, m.l1, m.r0, m.strand, m.depth) == (e.start, e.end, e.l1, e.r0, "+", 0)
    assert got[fx.kind("host").name].depth == 1 and len(got) == len(fx.elements)


def test_fetch_records_returns_only_the_named(fx):
    names = {fx.kind("del_left").name}
    got = rio.fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")], names)
    assert set(got) == names


def test_forward_and_stored_leave_depth_letters_alone():
    rec = "ACGTNRDYACGT"
    assert rio.forward(rec, "-") == revcomp_record(rec)
    assert rio.stored(rio.forward(rec, "-"), "-") == rec
    assert rio.forward(rec, "+") == rec


def test_sanitize_turns_every_ambiguity_into_n():
    assert rio.sanitize("acgtNRDY") == "ACGTNNNN"


def test_paint_mirrors_on_reverse_stored_records():
    rec = "A" * 20
    assert rio.paint(rec, "+", 101, [(103, 105)], "R") == "AA" + "RRR" + "A" * 15
    assert rio.paint(rec, "-", 101, [(103, 105)], "R") == "A" * 15 + "RRR" + "AA"


def test_rekey_nest_and_hosts():
    value = "nest-outer:c:1-9;nest-inner:c:50-900"
    assert rio.rekey_nest(value, {"c:1-9": "c:1-12"}) == "nest-outer:c:1-12;nest-inner:c:50-900"
    assert rio.hosts_of(value) == ["c:50-900"]
    assert rio.rekey_nest(".", {"c:1-9": "c:1-12"}) == "."


def test_conflict_rules(fx):
    ms = {m.name: m for m in members(fx)}
    index = rio.SpanIndex(ms.values())
    inner = ms[fx.kind("nested_del_left").name]
    host = ms[fx.kind("host").name]
    assert rio.conflict(index, inner, inner.start - 50, inner.end) is None
    assert rio.conflict(index, inner, host.start - 10, inner.end) == "host_exceeded"
    c = ms[fx.kind("conflict_left").name]
    assert rio.conflict(index, c, fx.kind("conflict_left").true_start, c.end) == "overlaps_element"
    x = fx.kind("neighbour")
    assert rio.conflict(index, c, x.start - 5, c.end) == "engulfs_element"


def test_mutual_conflicts_keep_the_first_claim(fx):
    a = members(fx)[0]
    b = dataclasses.replace(a, name="chrS:900000-900500#x", start=a.end + 100, end=a.end + 600,
                            l1=a.end + 200, r0=a.end + 500)
    lost = rio.mutual_conflicts([(a, a.start, a.end + 60), (b, b.start - 60, b.end)])
    assert lost == {b.uid}


def test_sidecar_round_trips_through_read_map(tmp_path):
    row = {c: "." for c in rio.SIDECAR_COLUMNS}
    row.update(old_seq_id="c:150-900#LTR/Gypsy/X", new_seq_id="c:100-900#LTR/Gypsy/X",
               decision="extended", ext5="50", ext3="0")
    rejected = dict(row, new_seq_id=".", decision="rejected", reason="gate_identity",
                    old_seq_id="c:2000-2500#LTR/Gypsy/X")
    path = tmp_path / "p_reboundary.tsv"
    path.write_text(rio.sidecar_text([row, rejected]))
    got = rio.read_map(str(path))
    assert set(got) == {"c:100-900"}
    assert got["c:100-900"] == rio.Rebound("c:150-900#LTR/Gypsy/X", "c:100-900#LTR/Gypsy/X", 50, 0)


def test_commit_is_all_or_nothing(tmp_path):
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    a.write_text("old-a")
    c = rio.Commit()
    with c.open(str(a)) as fh:
        fh.write("new-a")
    with c.open(str(b)) as fh:
        fh.write("new-b")
    assert a.read_text() == "old-a" and not b.exists()
    c.abort()
    assert a.read_text() == "old-a" and not b.exists() and not list(tmp_path.glob("*.rb.tmp"))
    c2 = rio.Commit()
    with c2.open(str(a)) as fh:
        fh.write("new-a")
    c2.commit()
    assert a.read_text() == "new-a"


def test_rewrite_renames_rekeys_and_paints_the_host(fx):
    tables = rio.load_clean_tables(str(fx.indir), fx.prefix)
    ms = {m.name: m for m in rio.members_from(fx.prefix, tables)}
    e = fx.kind("nested_del_left")
    m = ms[e.name]
    new_name = f"chrS:{e.true_start}-{e.end}#LTR/Gypsy/Synth"
    old_row = [r for r in tables[0].rows if r[0] == e.name][0]
    fields = [new_name] + old_row[1:len(COLUMNS)]
    genome = "".join(fx.genome.read_text().splitlines()[1:])
    record = genome[e.true_start - 1:e.end]
    acc = rio.Accepted(m, new_name, fields, record, e.true_start, e.end)
    raw_before = (fx.indir / f"{fx.prefix}_depth0_ltr.tsv").read_bytes()
    c = rio.Commit()
    rio.rewrite(tables, {m.key: acc}, IUPAC_DEPTH_SEQ, c)
    c.commit()

    t0 = (fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv").read_text()
    assert new_name in t0 and e.name not in t0
    t1 = (fx.indir / f"{fx.prefix}_depth1_clean_ltr.tsv").read_text()
    assert f"nest-outer:chrS:{e.true_start}-{e.end}" in t1
    fa0 = rio.fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")], {new_name})
    assert fa0[new_name] == record
    host = fx.kind("host")
    fa1 = rio.fetch_records([str(fx.indir / f"{fx.prefix}_depth1_clean_ltr.fa")], {host.name})
    painted = fa1[host.name][e.true_start - host.start:e.end - host.start + 1]
    assert set(painted) == {"N"}
    assert (fx.indir / f"{fx.prefix}_depth0_ltr.tsv").read_bytes() == raw_before
```

- [ ] **Step 2: Run to verify failure**

Run: `PYT tests/test_reboundary_io.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'ltrquest.reboundary_io'`.

- [ ] **Step 3: Implement `src/ltrquest/reboundary_io.py`**

```python
"""Reading and rewriting a run's clean side for ltrquest.reboundary.

Only the `_clean_` depth tables and FASTAs are ever rewritten. The raw
`_depth<N>_ltr.*` files are detection outputs that ltrquest.record hands to a
later run, so they keep the calls as detected. Every rewritten file is staged
beside its target and renamed only once every file has been written, so a
failure leaves the run exactly as it was.
"""

from __future__ import annotations

import bisect
import os
import re
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, NamedTuple, Optional, Sequence, Set, Tuple

from .annotate import discover_depth_tables, element_key, header_names, read_table, target_mode
from .detect import revcomp as revcomp_record
from .kmer2ltr import COLUMNS
from .ltr_model import Member
from .table import Columns, as_float, as_int

SIDECAR_SUFFIX = "_reboundary.tsv"
SIDECAR_COLUMNS = [
    "old_seq_id", "new_seq_id", "family", "method", "model_id", "decision", "reason",
    "ext5", "ext3", "end_source5", "end_source3", "id_outer5", "id_outer3", "credit_bits",
    "k2l_status", "tsd_called", "tsd_new", "tsd_null", "k2p_called", "k2p_new",
    "obstacle5", "obstacle3",
]
NEST_RE = re.compile(r"^(nest-outer|nest-inner):(.+)$")


@dataclass
class CleanTable:
    path: str
    fasta: str
    depth: int
    header: List[str]
    rows: List[List[str]]
    cols: Columns


def load_clean_tables(indir: str, prefix: str) -> List[CleanTable]:
    out: List[CleanTable] = []
    for t in discover_depth_tables(prefix, indir, variants=("clean",)):
        header, rows = read_table(t.path)
        if header is None:
            raise ValueError(f"{t.path}: no header line")
        cols = Columns.of(header_names(header))
        if cols.names[:len(COLUMNS)] != COLUMNS:
            raise ValueError(f"{t.path}: its first {len(COLUMNS)} columns are not Kmer2LTR's, "
                             f"and re-boundarying replaces them by position")
        missing = [c for c in ("strand", "family", "nest_status", "orientation") if c not in cols]
        if missing:
            raise ValueError(f"{t.path}: no {', '.join(missing)} column; run ltrquest.annotate "
                             f"first (re-boundarying models families)")
        out.append(CleanTable(t.path, t.path[:-len(".tsv")] + ".fa", t.depth, header, rows, cols))
    return out


def members_from(prefix: str, tables: Sequence[CleanTable]) -> List[Member]:
    out: List[Member] = []
    for t in tables:
        for row in t.rows:
            key = element_key(row[0]) if row else None
            if key is None:
                continue
            chrom, span = key.rsplit(":", 1)
            start, end = (int(x) for x in span.split("-"))
            l5e, l3s = as_int(t.cols.get(row, "ltr5_end")), as_int(t.cols.get(row, "ltr3_start"))
            if l5e is None or l3s is None:
                continue
            out.append(Member(
                prefix=prefix, name=row[0], chrom=chrom, start=start, end=end,
                l1=start + l5e - 1, r0=start + l3s - 1,
                strand=t.cols.get(row, "strand"), orientation=t.cols.get(row, "orientation", "+"),
                family=t.cols.get(row, "family"), depth=t.depth,
                k2p=as_float(t.cols.get(row, "k2p")), tsd=t.cols.get(row, "tsd"),
                nest_status=t.cols.get(row, "nest_status")))
    return out


def iter_fasta(path: str) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].rstrip("\r\n"), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def fetch_records(paths: Iterable[str], names: Set[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for path in paths:
        if names and os.path.isfile(path):
            for name, seq in iter_fasta(path):
                if name in names:
                    out[name] = seq
    return out


def forward(record: str, orientation: str) -> str:
    """A stored record turned to the genome-forward frame (depth letters untouched)."""
    return record if orientation != "-" else revcomp_record(record)


def stored(fwd: str, orientation: str) -> str:
    return fwd if orientation != "-" else revcomp_record(fwd)


def sanitize(seq: str) -> str:
    """Every non-ACGT base as N, as Kmer2LTR's own FASTA reader does."""
    return re.sub(r"[^ACGT]", "N", seq.upper())


def paint(record: str, orientation: str, rec_start: int, spans: Sequence[Tuple[int, int]],
          letter: str) -> str:
    """Overwrite genomic spans (1-based, inclusive) inside a stored record with a depth letter."""
    chars, n = list(record), len(record)
    for a, b in spans:
        lo, hi = max(0, a - rec_start), min(n, b - rec_start + 1)
        if hi <= lo:
            continue
        if orientation == "-":
            lo, hi = n - hi, n - lo
        for i in range(lo, hi):
            chars[i] = letter
    return "".join(chars)


def rekey_nest(value: str, old2new: Dict[str, str]) -> str:
    if value in ("", ".") or not old2new:
        return value
    out = []
    for tok in value.split(";"):
        m = NEST_RE.match(tok)
        if m and m.group(2) in old2new:
            tok = f"{m.group(1)}:{old2new[m.group(2)]}"
        out.append(tok)
    return ";".join(out)


def hosts_of(nest_status: str) -> List[str]:
    """Keys of the elements this one is nested inside."""
    return [m.group(2) for m in (NEST_RE.match(t) for t in nest_status.split(";"))
            if m and m.group(1) == "nest-inner"]


class SpanIndex:
    """Every element's span per (prefix, chrom), for overlap queries."""

    def __init__(self, members: Iterable[Member]):
        by: Dict[Tuple[str, str], List[Tuple[int, int, str]]] = defaultdict(list)
        for m in members:
            by[(m.prefix, m.chrom)].append((m.start, m.end, m.key))
        self._by = {k: sorted(v) for k, v in by.items()}
        self._starts = {k: [s for s, _, _ in v] for k, v in self._by.items()}
        self._longest = {k: max(e - s + 1 for s, e, _ in v) for k, v in self._by.items()}

    def overlapping(self, prefix: str, chrom: str, lo: int, hi: int) -> List[Tuple[int, int, str]]:
        k = (prefix, chrom)
        spans = self._by.get(k, [])
        out = []
        for s, e, key in spans[bisect.bisect_left(self._starts.get(k, []), lo - self._longest.get(k, 0)):]:
            if s > hi:
                break
            if e >= lo:
                out.append((s, e, key))
        return out


def _added(m: Member, left: int, right: int) -> List[Tuple[int, int]]:
    return (([(left, m.start - 1)] if left < m.start else [])
            + ([(m.end + 1, right)] if right > m.end else []))


def conflict(index: SpanIndex, m: Member, left: int, right: int) -> Optional[str]:
    """Why extending `m` to left..right would clash with another element, or None."""
    added = _added(m, left, right)
    for s, e, key in index.overlapping(m.prefix, m.chrom, left, right):
        if key == m.key:
            continue
        if s <= m.start and e >= m.end:           # a host of the call as it stands
            if not (s <= left and e >= right):
                return "host_exceeded"
            continue
        if s >= m.start and e <= m.end:           # nested inside the call already
            continue
        for a, b in added:
            if s >= a and e <= b:
                return "engulfs_element"
            if s <= b and e >= a:
                return "overlaps_element"
    return None


def mutual_conflicts(items: Sequence[Tuple[Member, int, int]]) -> Set[str]:
    """uids whose added bases overlap an earlier candidate's added bases (genome, chrom, start order)."""
    taken: Dict[Tuple[str, str], List[Tuple[int, int]]] = defaultdict(list)
    lost: Set[str] = set()
    for m, left, right in sorted(items, key=lambda t: (t[0].prefix, t[0].chrom, t[0].start)):
        added = _added(m, left, right)
        k = (m.prefix, m.chrom)
        if any(a <= d and b >= c for a, b in added for c, d in taken[k]):
            lost.add(m.uid)
            continue
        taken[k].extend(added)
    return lost


@dataclass(frozen=True)
class Accepted:
    member: Member
    new_name: str
    fields: List[str]       # the 29 Kmer2LTR columns, rebased
    record: str             # stored-orientation FASTA record
    left: int
    right: int


class Commit:
    """Stage each rewritten file beside its target; rename them all only at the end."""

    def __init__(self):
        self._staged: List[Tuple[str, str]] = []

    def open(self, path: str):
        directory = os.path.dirname(os.path.abspath(path))
        fd, tmp = tempfile.mkstemp(prefix=os.path.basename(path) + ".", suffix=".rb.tmp",
                                   dir=directory)
        self._staged.append((tmp, path))
        return os.fdopen(fd, "w")

    def commit(self) -> None:
        for tmp, path in self._staged:
            os.chmod(tmp, target_mode(path))
            os.replace(tmp, path)
        self._staged = []

    def abort(self) -> None:
        for tmp, _ in self._staged:
            try:
                os.unlink(tmp)
            except OSError:
                pass
        self._staged = []


def atomic_write_text(path: str, text: str) -> None:
    c = Commit()
    with c.open(path) as fh:
        fh.write(text)
    c.commit()


def rewrite(tables: Sequence[CleanTable], accepted: Dict[str, Accepted], letters: Sequence[str],
            commit: Commit, wrap: int = 60) -> None:
    """Stage one genome's clean tables and FASTAs with `accepted` (old key -> result) applied."""
    old2new = {k: element_key(a.new_name) for k, a in accepted.items()}
    jobs: Dict[str, List[Tuple[List[Tuple[int, int]], str]]] = defaultdict(list)
    for a in accepted.values():
        m = a.member
        letter = letters[m.depth] if m.depth < len(letters) else "X"
        for host in hosts_of(m.nest_status):
            jobs[host].append((_added(m, a.left, a.right), letter))
    where: Dict[str, Tuple[int, str]] = {}
    for t in tables:
        for row in t.rows:
            k = element_key(row[0]) if row else None
            if k is not None:
                start = accepted[k].left if k in accepted else int(k.rsplit(":", 1)[1].split("-")[0])
                where[k] = (start, t.cols.get(row, "orientation", "+"))
    for t in tables:
        i_nest = t.cols.require("nest_status")
        rows = []
        for row in t.rows:
            k = element_key(row[0]) if row else None
            if k in accepted:
                row = list(accepted[k].fields) + row[len(COLUMNS):]
            else:
                row = list(row)
            row[i_nest] = rekey_nest(row[i_nest], old2new)
            rows.append(row)
        if os.path.isfile(t.fasta):
            with commit.open(t.fasta) as out:
                for name, seq in iter_fasta(t.fasta):
                    k = element_key(name)
                    if k in accepted:
                        name, seq = accepted[k].new_name, accepted[k].record
                    for spans, letter in jobs.get(k, ()):
                        start, orientation = where[k]
                        seq = paint(seq, orientation, start, spans, letter)
                    out.write(f">{name}\n")
                    for i in range(0, len(seq), wrap):
                        out.write(seq[i:i + wrap] + "\n")
        with commit.open(t.path) as out:
            out.write("\t".join(t.header) + "\n")
            for row in rows:
                out.write("\t".join(row) + "\n")


def sidecar_text(rows: Sequence[Dict[str, str]]) -> str:
    lines = ["#" + "\t".join(SIDECAR_COLUMNS)]
    lines += ["\t".join(str(r.get(c, ".")) for c in SIDECAR_COLUMNS) for r in rows]
    return "\n".join(lines) + "\n"


class Rebound(NamedTuple):
    old_name: str
    new_name: str
    ext5: int
    ext3: int


def read_map(path: str) -> Dict[str, Rebound]:
    """New element key -> Rebound, for the extended rows of a sidecar."""
    header, rows = read_table(path)
    cols = Columns.of(header_names(header))
    out: Dict[str, Rebound] = {}
    for row in rows:
        if cols.get(row, "decision") != "extended":
            continue
        new = cols.get(row, "new_seq_id")
        key = element_key(new)
        if key is not None:
            out[key] = Rebound(cols.get(row, "old_seq_id"), new,
                               as_int(cols.get(row, "ext5"), 0), as_int(cols.get(row, "ext3"), 0))
    return out
```

- [ ] **Step 4: Run the tests**

Run: `PYT tests/test_reboundary_io.py -q`
Expected: 13 passed.

- [ ] **Step 5: Commit**

```bash
git add src/ltrquest/reboundary_io.py tests/test_reboundary_io.py
git commit -m "reboundary_io: clean-side loading, record splicing, conflicts, sidecar, atomic rewrite

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: `reboundary.py` — the pooled driver with Kmer2LTR arbitration, and the CLI

**Files:**
- Create: `src/ltrquest/reboundary.py`
- Test: `tests/test_reboundary.py`

**Interfaces:**
- Consumes: Tasks 1, 3, 4, 5.
- Produces:
  - `Settings(method, references, min_copies, credit, max_ratio, qc_min_n, untested, subfamily_jaccard, mutation_rate, mafft, threads, place: Params)` with `family_key()`
  - `METHODS` tuple; `_G`, `_API`, `_init(paths, tools_dir)` (worker set-up, also used by the benchmark)
  - `build_models(family, members, g, s, exclude=frozenset()) -> (models, modal, status)`
  - `FamilyResult`; `family_job((family, members, s)) -> FamilyResult`
  - `qc(results, s) -> (status_by_family, p0)`; `family_accepts(status, s) -> bool`
  - `credit_for(p, s) -> float`
  - `Verdict(uid, status, accepted, k2l_status, tsd, k2p, credit)`; `arbitrate((p, record, credit, mu, min_ext)) -> Verdict`; `rebase(fields, f5, name)`; `length_ok(l5, l3, aln, left, right)`
  - `sidecar_row(p, fr, s) -> Dict[str, str]`
  - `run(indir, prefixes, genomes, s, tools_dir, dump=None, cache=None, verbose=False) -> Dict[str, int]` (counts: `extended` and every reject reason)
  - `resolve_mutation_rate(indir, prefixes, given) -> float`; `build_parser()`, `settings_from(args)`, `main(argv=None) -> int`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_reboundary.py`:

```python
"""The driver: family QC, Kmer2LTR arbitration, and what a real rewrite leaves behind."""

from __future__ import annotations

import json
import os

import pytest

from ltrquest import reboundary as rb
from ltrquest.detect import revcomp as revcomp_record
from ltrquest.kmer2ltr import COLUMNS
from ltrquest.ltr_model import Member, Model
from ltrquest.ltr_place import Proposal
from ltrquest.reboundary_io import fetch_records

from reboundary_fixtures import build

TOOLS = os.environ.get("LTRQUEST_TOOLS_DIR", ".")


def fake_member(**kw):
    base = dict(prefix="p", name="c:100-900#LTR/Gypsy/X", chrom="c", start=100, end=900,
                l1=199, r0=801, strand="+", orientation="+", family="f", depth=0, k2p=0.01,
                tsd=".")
    base.update(kw)
    return Member(**base)


def fake_proposal(m, left, right, gate=True, id_left=0.9, id_right=0.9):
    return Proposal(m, "f:consensus", "+", left, right, "core", "core", id_left, id_right, 0.9,
                    left, right, m.l1, m.r0, gate)


def test_rebase_restates_coordinates_against_the_cut_record():
    f = ["NA"] * len(COLUMNS)
    i = {c: k for k, c in enumerate(COLUMNS)}
    f[i["ltr5_start"]], f[i["ltr5_end"]] = "11", "210"
    f[i["ltr3_start"]], f[i["ltr3_end"]] = "811", "1010"
    f[i["flank5_len"]], f[i["flank3_len"]] = "10", "5"
    out = rb.rebase(f, 10, "c:110-1109#x")
    assert out[i["seq_id"]] == "c:110-1109#x"
    assert [out[i[c]] for c in ("ltr5_start", "ltr5_end", "ltr3_start", "ltr3_end", "seq_len")] \
        == ["1", "200", "801", "1000", "1000"]
    assert out[i["flank5_len"]] == out[i["flank3_len"]] == "0"


def test_length_gate_matches_detection():
    assert rb.length_ok(400, 390, 405, 1, 3000)
    assert not rb.length_ok(400, 150, 405, 1, 3000)      # ratio < 0.65
    assert not rb.length_ok(90, 90, 95, 1, 3000)         # LTR < 100
    assert not rb.length_ok(400, 400, 400, 1, 200)       # element < 300


def test_credit_for_constant_and_model():
    m = fake_member()
    p = fake_proposal(m, 50, 900, id_left=0.9)
    assert rb.credit_for(p, rb.Settings(credit="200")) == 200.0
    s = rb.Settings(credit="model")
    assert rb.credit_for(p, s) == pytest.approx(2 * s.place.anchor_len * 0.9)


def test_qc_passes_families_whose_proposals_gain_tsds():
    fam = rb.FamilyResult("f", 12, [Model("f:consensus", "f", "A" * 400, 400.0, 12)], 400.0, "ok")
    for k in range(6):
        m = fake_member(name=f"c:{1000 * k + 100}-{1000 * k + 900}#x", start=1000 * k + 100,
                        end=1000 * k + 900, l1=1000 * k + 199, r0=1000 * k + 801)
        fam.proposals.append(fake_proposal(m, m.start - 50, m.end))
        fam.tsd_new[m.uid], fam.tsd_null[m.uid] = "ACGTA", "."
    status, p0 = rb.qc([fam], rb.Settings())
    assert status["f"] == "pass" and p0 == 0.0
    assert rb.family_accepts("pass", rb.Settings())
    assert rb.family_accepts("untested", rb.Settings(untested="accept"))
    assert not rb.family_accepts("untested", rb.Settings(untested="skip"))
    assert not rb.family_accepts("fail", rb.Settings())


def test_mutation_rate_comes_from_the_run_record(tmp_path):
    (tmp_path / "p.detect.json").write_text(json.dumps({"settings": {"mutation_rate": "1.3e-8"}}))
    assert rb.resolve_mutation_rate(str(tmp_path), ["p"], None) == pytest.approx(1.3e-8)
    assert rb.resolve_mutation_rate(str(tmp_path), ["p"], 5e-9) == 5e-9
    assert rb.resolve_mutation_rate(str(tmp_path / "none"), ["p"], None) == 3e-8


def test_a_non_numeric_credit_is_refused():
    with pytest.raises(SystemExit):
        rb.main(["--prefix", "p", "--genome", "g.fa", "--credit", "lots"])


def sidecar(fx):
    lines = (fx.indir / f"{fx.prefix}_reboundary.tsv").read_text().splitlines()
    head = lines[0].lstrip("#").split("\t")
    rows = [dict(zip(head, l.split("\t"))) for l in lines[1:]]
    return {r["old_seq_id"].split("#")[0]: r for r in rows}


@pytest.fixture(scope="module")
def ran(tmp_path_factory, k2l_api, mafft):
    fx = build(tmp_path_factory.mktemp("syn_run"))
    raw = {p.name: p.read_bytes() for p in fx.indir.glob(f"{fx.prefix}_depth*_ltr.*")
           if "_clean_" not in p.name}
    counts = rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)],
                    rb.Settings(threads=2, mafft=mafft), TOOLS)
    return fx, counts, raw


@pytest.mark.parametrize("kind", ["del_left", "patch_right", "del_left_minus", "nested_del_left"])
def test_truncated_calls_are_extended_to_their_true_ends(ran, kind):
    fx, _, _ = ran
    e = fx.kind(kind)
    row = sidecar(fx)[e.key]
    assert row["decision"] == "extended", row
    assert row["new_seq_id"].split("#")[0] == f"chrS:{e.true_start}-{e.true_end}"
    assert row["tsd_new"] not in (".", "NA")


def test_a_clash_with_a_neighbour_is_rejected(ran):
    fx, _, _ = ran
    assert sidecar(fx)[fx.kind("conflict_left").key]["reason"] == "overlaps_element"


def test_an_insertion_beyond_ltrquests_pair_rules_is_rejected_with_a_reason(ran):
    fx, _, _ = ran
    row = sidecar(fx)[fx.kind("ins_left").key]
    assert row["decision"] == "rejected" and row["end_source5"] == "anchor"
    assert row["reason"] in ("kmer2ltr_reverted", "length_filter", "kmer2ltr_not_pass")


def test_intact_copies_are_not_candidates(ran):
    fx, _, _ = ran
    side = sidecar(fx)
    assert not [e for e in fx.elements if e.kind in ("normal", "decayed") and e.key in side]


def test_rewritten_rows_come_from_kmer2ltr_and_keep_the_annotation(ran):
    fx, _, _ = ran
    e = fx.kind("del_left")
    new = f"chrS:{e.true_start}-{e.true_end}#LTR/Gypsy/Synth"
    lines = (fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv").read_text().splitlines()
    head = lines[0].lstrip("#").split("\t")
    row = dict(zip(head, [l for l in lines if l.startswith(new + "\t")][0].split("\t")))
    assert row["ltr5_start"] == "1" and row["seq_len"] == str(e.true_end - e.true_start + 1)
    assert row["flank5_len"] == row["flank3_len"] == "0" and row["status"] == "pass"
    assert row["orientation"] == "+" and row["strand"] == "+" and row["family"] == e.family
    assert row["tsd"] not in (".", "NA")


def test_a_minus_record_is_stored_reverse_complemented(ran):
    fx, _, _ = ran
    e = fx.kind("del_left_minus")
    new = f"chrS:{e.true_start}-{e.true_end}#LTR/Gypsy/Synth"
    genome = "".join(fx.genome.read_text().splitlines()[1:])
    rec = fetch_records([str(fx.indir / f"{fx.prefix}_depth0_clean_ltr.fa")], {new})[new]
    assert rec == revcomp_record(genome[e.true_start - 1:e.true_end])


def test_the_host_masks_the_recovered_bases_and_points_at_the_new_name(ran):
    fx, _, _ = ran
    e, host = fx.kind("nested_del_left"), fx.kind("host")
    rec = fetch_records([str(fx.indir / f"{fx.prefix}_depth1_clean_ltr.fa")], {host.name})[host.name]
    assert set(rec[e.true_start - host.start:e.start - host.start]) == {"N"}
    t1 = (fx.indir / f"{fx.prefix}_depth1_clean_ltr.tsv").read_text()
    assert f"nest-outer:chrS:{e.true_start}-{e.true_end}" in t1


def test_raw_tables_are_untouched(ran):
    fx, _, raw = ran
    for name, data in raw.items():
        assert (fx.indir / name).read_bytes() == data, name


def test_genomes_from_different_family_namespaces_are_refused(tmp_path, k2l_api, mafft):
    fx = build(tmp_path / "syn_ns")
    path = fx.indir / f"{fx.prefix}_depth0_clean_ltr.tsv"
    lines = path.read_text().splitlines()
    lines[1] = lines[1].replace("\tsyn_fam00001\t", "\tother_fam00001\t")
    path.write_text("\n".join(lines) + "\n")
    with pytest.raises(SystemExit, match="namespaces"):
        rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)], rb.Settings(threads=1, mafft=mafft),
               TOOLS)


def test_a_second_run_extends_nothing_and_the_cli_works(tmp_path, k2l_api, mafft):
    fx = build(tmp_path / "syn_twice")
    args = ["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
            "-t", "2", "--mafft", mafft, "--tools-dir", TOOLS]
    assert rb.main(args) == 0
    first = sidecar(fx)
    assert sum(r["decision"] == "extended" for r in first.values()) == 4
    counts = rb.run(str(fx.indir), [fx.prefix], [str(fx.genome)],
                    rb.Settings(threads=2, mafft=mafft), TOOLS)
    assert counts.get("extended", 0) == 0
```

- [ ] **Step 2: Run to verify failure**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT tests/test_reboundary.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'ltrquest.reboundary'`.

- [ ] **Step 3: Implement `src/ltrquest/reboundary.py`**

```python
"""Family-guided re-boundarying of truncated LTR-RT calls (opt-in stage).

LTRharvest and LTR_FINDER extend an LTR pair outward from a seed and stop at the
first obstacle between the element's two LTRs -- an indel, a patch of mutations.
Kmer2LTR can only trim. So an element whose LTR carries an obstacle near its end
is called short. This stage builds each family's LTR model from its full-length
copies (pooled over every genome), places it on every member, and proposes
outward-only ends. Candidates pass family QC and conflict checks, then Kmer2LTR
re-scores each widened record with the family's evidence as `tsd_credit`: it
settles the final ends and every Kmer2LTR column. Only the `_clean_` tables and
FASTAs are rewritten; `<prefix>_reboundary.tsv` records every candidate and why
it was or was not extended, and doubles as the old->new key map for the GFF3.

Usage:
  ltrquest-reboundary --indir RUN --prefix P1 [P2 ...] --genome G1 [G2 ...] [options]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import shutil
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from multiprocessing import Pool
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

from . import kmer2ltr as k2l
from .kmer2ltr import COLUMNS
from .ltr_model import (Genomes, Member, Model, consensus_model, ratio_ok, select_references,
                        tsd_enrichment)
from .ltr_place import Params, Proposal, obstacle, propose, tsd_at
from .reboundary_io import (SIDECAR_SUFFIX, Accepted, Commit, SpanIndex, conflict,
                            fetch_records, forward, load_clean_tables, members_from,
                            mutual_conflicts, rewrite, sanitize, sidecar_text, stored)
from .reconcile import IUPAC_DEPTH_SEQ

METHODS = ("consensus",)
_I = {c: i for i, c in enumerate(COLUMNS)}


def log(msg: str) -> None:
    print(f"[reboundary] {msg}", flush=True)


def warn(msg: str) -> None:
    print(f"[reboundary] WARNING: {msg}", file=sys.stderr, flush=True)


@dataclass(frozen=True)
class Settings:
    method: str = "consensus"
    references: str = "modal"
    min_copies: int = 10
    credit: str = "200"
    max_ratio: float = 1.15
    qc_min_n: int = 5
    untested: str = "accept"
    subfamily_jaccard: float = 0.5
    mutation_rate: float = 3e-8
    mafft: str = "mafft"
    threads: int = 8
    place: Params = field(default_factory=Params)

    def family_key(self) -> str:
        """Everything the family phase depends on, so a cached phase is reused only when valid."""
        # min_copies is left out on purpose: it only chooses which families run, and a
        # cached phase covering more families is reused for fewer (see _family_phase).
        keep = ("method", "references", "subfamily_jaccard", "max_ratio", "place")
        d = {k: v for k, v in asdict(self).items() if k in keep}
        return hashlib.sha1(json.dumps(d, sort_keys=True).encode()).hexdigest()


_G: Optional[Genomes] = None
_API = None


def _init(paths: Dict[str, str], tools_dir: str) -> None:
    """Worker set-up: genome handles and the Kmer2LTR API, once per process."""
    global _G, _API
    _G = Genomes(paths)
    _API = k2l.api(tools_dir)


def build_models(family: str, members: Sequence[Member], g: Genomes, s: Settings,
                 exclude: FrozenSet[str] = frozenset()) -> Tuple[List[Model], Optional[float], str]:
    refs, modal = select_references(members, s.references, family, exclude)
    if not refs:
        return [], modal, "too_few_references"
    if s.method == "consensus":
        model = consensus_model(family, refs, g, modal, s.mafft)
        return ([model] if model else []), modal, ("ok" if model else "no_ltr_span")
    raise ValueError(f"unknown method {s.method!r}")


@dataclass
class FamilyResult:
    family: str
    n_members: int
    models: List[Model]
    modal: Optional[float]
    status: str
    proposals: List[Proposal] = field(default_factory=list)
    tsd_new: Dict[str, str] = field(default_factory=dict)
    tsd_null: Dict[str, str] = field(default_factory=dict)
    obstacles: Dict[str, Tuple[str, str]] = field(default_factory=dict)


def family_job(args) -> FamilyResult:
    family, members, s = args
    models, modal, status = build_models(family, members, _G, s)
    fr = FamilyResult(family, len(members), models, modal, status)
    if not models:
        return fr
    kept = [mo for mo in models if ratio_ok(mo, s.max_ratio)]
    if not kept:
        fr.status = "ratio_failed"   # proposals still computed, so the sidecar says what failed
    for m in members:
        p = propose(m, kept or models, _G, s.place)
        if p is None:
            continue
        fr.proposals.append(p)
        if p.gate_ok:
            fr.tsd_new[m.uid] = tsd_at(_API, _G, m, p.left, p.right)
            fr.tsd_null[m.uid] = tsd_at(_API, _G, m, p.left, p.right, shift=1000)
            fr.obstacles[m.uid] = obstacle(m, p, _G)
    return fr


def _found(tsd: str) -> bool:
    return tsd not in ("", ".", "NA")


def qc(results: Sequence[FamilyResult], s: Settings) -> Tuple[Dict[str, str], float]:
    """Per-family QC status and the pooled displaced-flank null rate it was tested against.

    Evidence: TSDs at the proposed ends of candidates whose called ends had none
    (TSDs are a credibility signal the finders did not select on; TG..CA is not
    used because they did).
    """
    null_hits = null_n = 0
    for fr in results:
        for p in fr.proposals:
            u = p.member.uid
            if p.gate_ok and not p.member.has_tsd and u in fr.tsd_null:
                null_n += 1
                null_hits += _found(fr.tsd_null[u])
    p0 = null_hits / null_n if null_n else 0.0
    status: Dict[str, str] = {}
    for fr in results:
        if fr.status != "ok":
            status[fr.family] = fr.status
            continue
        ev = [fr.tsd_new[p.member.uid] for p in fr.proposals
              if p.gate_ok and not p.member.has_tsd and p.member.uid in fr.tsd_new]
        status[fr.family] = tsd_enrichment(sum(map(_found, ev)), len(ev), p0, s.qc_min_n)
    return status, p0


def family_accepts(status: str, s: Settings) -> bool:
    return status == "pass" or (status == "untested" and s.untested == "accept")


def credit_for(p: Proposal, s: Settings) -> float:
    """Bits of family evidence handed to Kmer2LTR that the proposed termini are the ends."""
    if s.credit == "model":
        ids = (([p.id_left] if p.left < p.member.start else [])
               + ([p.id_right] if p.right > p.member.end else []))
        return 2.0 * s.place.anchor_len * min(ids) if ids else 0.0
    return float(s.credit)


@dataclass(frozen=True)
class Verdict:
    uid: str
    status: str
    accepted: Optional[Accepted]
    k2l_status: str
    tsd: str
    k2p: str
    credit: float


def window(m: Member, left: int, right: int):
    """Kmer2LTR's four reference cuts around a record (genome._cuts), from pyfaidx."""
    pad, probe = _API.PAD, _API.PROBE
    return _API.Window(_G.fetch(m.prefix, m.chrom, left - pad, left - 1)[0],
                       _G.fetch(m.prefix, m.chrom, left, left + probe - 1)[0],
                       _G.fetch(m.prefix, m.chrom, right - probe + 1, right)[0],
                       _G.fetch(m.prefix, m.chrom, right + 1, right + pad)[0])


def rebase(fields: List[str], f5: int, name: str) -> List[str]:
    """Kmer2LTR's row restated against the record cut to its own bounds (detect.rebase_to_trimmed)."""
    out = list(fields)
    out[_I["seq_id"]] = name
    for c in ("ltr5_end", "ltr3_start", "ltr3_end"):
        out[_I[c]] = str(int(out[_I[c]]) - f5)
    out[_I["ltr5_start"]] = "1"
    out[_I["seq_len"]] = out[_I["ltr3_end"]]
    out[_I["flank5_len"]] = out[_I["flank3_len"]] = "0"
    return out


def length_ok(l5, l3, aln, left: int, right: int) -> bool:
    """Detection's Step 8a gate (detect.filter_kmer2ltr_in_place) on the settled pair."""
    if l5 is None or l3 is None or aln is None:
        return False
    hi = max(l5, l3, aln)
    return (hi > 0 and min(l5, l3) >= 100 and aln >= 90
            and min(l5, l3, aln) / hi >= 0.65 and right - left >= 300)


def arbitrate(args) -> Verdict:
    """Kmer2LTR re-scores the widened record; it settles the final ends and every column."""
    p, record, credit, mu, min_ext = args
    m = p.member
    left_seg = _G.fetch(m.prefix, m.chrom, p.left, m.start - 1)[0] if p.left < m.start else ""
    right_seg = _G.fetch(m.prefix, m.chrom, m.end + 1, p.right)[0] if p.right > m.end else ""
    fwd = left_seg + forward(record, m.orientation) + right_seg
    if len(fwd) != p.right - p.left + 1:
        return Verdict(m.uid, "record_mismatch", None, "NA", "NA", "NA", credit)
    suffix = "#" + m.name.split("#", 1)[1] if "#" in m.name else ""
    seq = sanitize(fwd)
    ctx = _API.orient(seq, window(m, p.left, p.right))
    res = _API.classify(f"{m.chrom}:{p.left}-{p.right}{suffix}", seq, period_rule="outermost",
                        mutation_rate=mu, tsd_credit=credit)
    if ctx is not None:
        res = _API.annotate(res, seq, ctx, _API.Options())
    if res.status != "pass":
        return Verdict(m.uid, "kmer2ltr_not_pass", None, res.status, "NA", "NA", credit)
    row = _API.format_row(res).split("\t")
    f5, f3 = res.flank5_len, res.flank3_len
    fl, fr = p.left + f5, p.right - f3
    tsd, k2p = row[_I["tsd"]], row[_I["k2p"]]
    if fl > m.start or fr < m.end or (m.start - fl < min_ext and fr - m.end < min_ext):
        return Verdict(m.uid, "kmer2ltr_reverted", None, res.status, tsd, k2p, credit)
    if not length_ok(res.ltr5_len, res.ltr3_len, res.aln_len, fl, fr):
        return Verdict(m.uid, "length_filter", None, res.status, tsd, k2p, credit)
    fields = rebase(row, f5, f"{m.chrom}:{fl}-{fr}{suffix}")
    fields[_I["orientation"]] = m.orientation   # a storage fact: the record stays stored as it was
    acc = Accepted(m, fields[0], fields, stored(fwd[f5:len(fwd) - f3], m.orientation), fl, fr)
    return Verdict(m.uid, "extended", acc, res.status, tsd, k2p, credit)


def sidecar_row(p: Proposal, fr: FamilyResult, s: Settings) -> Dict[str, str]:
    m, u = p.member, p.member.uid
    ob_l, ob_r = fr.obstacles.get(u, (".", "."))
    plus = p.orient == "+"
    return {
        "old_seq_id": m.name, "new_seq_id": ".", "family": m.family, "method": s.method,
        "model_id": p.model_id, "decision": "rejected", "reason": ".",
        "ext5": str(p.ext5 if p.gate_ok else p.raw_ext5),
        "ext3": str(p.ext3 if p.gate_ok else p.raw_ext3),
        "end_source5": p.left_src if plus else p.right_src,
        "end_source3": p.right_src if plus else p.left_src,
        "id_outer5": f"{(p.id_left if plus else p.id_right):.3f}",
        "id_outer3": f"{(p.id_right if plus else p.id_left):.3f}",
        "credit_bits": ".", "k2l_status": ".",
        "tsd_called": m.tsd or ".", "tsd_new": fr.tsd_new.get(u, "."),
        "tsd_null": fr.tsd_null.get(u, "."),
        "k2p_called": "NA" if m.k2p is None else f"{m.k2p:g}", "k2p_new": ".",
        "obstacle5": ob_l if plus else ob_r, "obstacle3": ob_r if plus else ob_l,
    }


def _row_order(row: Dict[str, str]):
    chrom, span = row["old_seq_id"].split("#", 1)[0].rsplit(":", 1)
    return chrom, int(span.split("-")[0])


def _family_phase(pool, eligible, s: Settings, cache: Optional[str], verbose: bool):
    key = s.family_key()
    if cache and os.path.isfile(cache):
        with open(cache, "rb") as fh:
            saved_key, results = pickle.load(fh)
        wanted = {f for f, _ in eligible}
        if saved_key == key and wanted <= {r.family for r in results}:
            log(f"family phase reused from {cache}")
            return [r for r in results if r.family in wanted]
    results = []
    for fr in pool.imap_unordered(family_job, [(f, ms, s) for f, ms in eligible]):
        results.append(fr)
        if verbose:
            lens = ",".join(str(len(mo.seq)) for mo in fr.models) or "-"
            log(f"{fr.family}: {fr.n_members} copies, model {lens} bp, modal {fr.modal}, "
                f"{fr.status}, {len(fr.proposals)} proposal(s)")
    results.sort(key=lambda fr: fr.family)
    if cache:
        with open(cache, "wb") as fh:
            pickle.dump((key, results), fh)
    return results


def _dump(path: str, results, status, rows) -> None:
    with open(path, "w") as fh:
        fh.write("uid\tfamily\tfamily_status\tgate_ok\tcalled_tsd\ttsd_new\ttsd_null\t"
                 "decision\treason\n")
        for fr in results:
            for p in fr.proposals:
                u = p.member.uid
                r = rows[u]
                fh.write("\t".join([u.replace("\t", "|"), fr.family, status[fr.family],
                                    str(int(p.gate_ok)), str(int(p.member.has_tsd)),
                                    fr.tsd_new.get(u, "."), fr.tsd_null.get(u, "."),
                                    r["decision"], r["reason"]]) + "\n")


def run(indir: str, prefixes: Sequence[str], genomes: Sequence[str], s: Settings,
        tools_dir: str, dump: Optional[str] = None, cache: Optional[str] = None,
        verbose: bool = False) -> Dict[str, int]:
    t0 = time.time()
    if len(prefixes) != len(genomes):
        raise SystemExit("reboundary: --prefix and --genome must pair up one to one")
    if s.method not in METHODS:
        raise SystemExit(f"reboundary: unknown method {s.method!r}")
    if shutil.which(s.mafft) is None:
        raise SystemExit(f"reboundary: mafft not found ({s.mafft}); it is in environment.yml")
    for gpath in genomes:
        if not os.path.isfile(gpath):
            raise SystemExit(f"reboundary: genome not found: {gpath}")
    k2l.api(tools_dir)     # fail here, and clone at most once, before any worker starts
    tables = {p: load_clean_tables(indir, p) for p in prefixes}
    for p, t in tables.items():
        if not t:
            raise SystemExit(f"reboundary: no {p}_depth<N>_clean_ltr.tsv in {indir}")
    members = [m for p in prefixes for m in members_from(p, tables[p])]
    families: Dict[str, List[Member]] = defaultdict(list)
    for m in members:
        if m.family not in ("", ".", "NA"):
            families[m.family].append(m)
    spaces = {f.rsplit("_fam", 1)[0] for f in families if "_fam" in f}
    if len(spaces) > 1:
        raise SystemExit(f"reboundary: the prefixes carry {len(spaces)} family namespaces "
                         f"({', '.join(sorted(spaces))}); pool only genomes whose families "
                         f"were clustered together in one run")
    eligible = sorted(((f, ms) for f, ms in families.items() if len(ms) >= s.min_copies),
                      key=lambda kv: (-len(kv[1]), kv[0]))
    log(f"start: {len(members)} elements in {len(prefixes)} genome(s), {len(families)} "
        f"families, {len(eligible)} with >= {s.min_copies} copies; method {s.method}")

    rows: Dict[str, Dict[str, str]] = {}
    by_uid: Dict[str, Proposal] = {}
    with Pool(s.threads, initializer=_init, initargs=(dict(zip(prefixes, genomes)), tools_dir)) \
            as pool:
        results = _family_phase(pool, eligible, s, cache, verbose)
        status, p0 = qc(results, s)
        log(f"family phase: {sum(len(fr.proposals) for fr in results)} candidates; "
            f"QC {dict(Counter(status.values()))}; displaced-flank null {p0:.3f}")
        survivors: List[Proposal] = []
        for fr in results:
            accept = family_accepts(status[fr.family], s)
            for p in fr.proposals:
                u = p.member.uid
                by_uid[u] = p
                rows[u] = sidecar_row(p, fr, s)
                if not p.gate_ok:
                    rows[u]["reason"] = "gate_identity"
                elif not accept:
                    rows[u]["reason"] = ("family_untested" if status[fr.family] == "untested"
                                         else "family_qc_failed")
                else:
                    survivors.append(p)
        index = SpanIndex(members)
        clear: List[Proposal] = []
        for p in survivors:
            why = conflict(index, p.member, p.left, p.right)
            if why:
                rows[p.member.uid]["reason"] = why
            else:
                clear.append(p)
        lost = mutual_conflicts([(p.member, p.left, p.right) for p in clear])
        for u in lost:
            rows[u]["reason"] = "overlaps_element"
        final = [p for p in clear if p.member.uid not in lost]
        jobs = []
        for prefix in prefixes:
            mine = [p for p in final if p.member.prefix == prefix]
            recs = fetch_records([t.fasta for t in tables[prefix]], {p.member.name for p in mine})
            for p in mine:
                rec = recs.get(p.member.name)
                if rec is None:
                    rows[p.member.uid]["reason"] = "record_missing"
                else:
                    jobs.append((p, rec, credit_for(p, s), s.mutation_rate, s.place.min_ext))
        log(f"arbitration: {len(jobs)} candidate(s) to Kmer2LTR")
        verdicts = pool.map(arbitrate, jobs, chunksize=8) if jobs else []

    accepted: Dict[str, Dict[str, Accepted]] = defaultdict(dict)
    for v in verdicts:
        row = rows[v.uid]
        row["credit_bits"], row["k2l_status"] = f"{v.credit:g}", v.k2l_status
        if v.accepted is None:
            row["reason"] = v.status
            continue
        a = v.accepted
        ext_l, ext_r = a.member.start - a.left, a.right - a.member.end
        plus = by_uid[v.uid].orient == "+"
        row.update(decision="extended", reason=".", new_seq_id=a.new_name, tsd_new=v.tsd,
                   k2p_new=v.k2p, ext5=str(ext_l if plus else ext_r),
                   ext3=str(ext_r if plus else ext_l))
        accepted[a.member.prefix][a.member.key] = a

    commit = Commit()
    try:
        for prefix in prefixes:
            if accepted[prefix]:
                rewrite(tables[prefix], accepted[prefix], IUPAC_DEPTH_SEQ, commit)
            mine = sorted((r for u, r in rows.items() if u.split("\t", 1)[0] == prefix),
                          key=_row_order)
            with commit.open(os.path.join(indir, prefix + SIDECAR_SUFFIX)) as fh:
                fh.write(sidecar_text(mine))
        commit.commit()
    except BaseException:
        commit.abort()
        raise
    if dump:
        _dump(dump, results, status, rows)
    counts = Counter(r["decision"] if r["decision"] == "extended" else r["reason"]
                     for r in rows.values())
    log(f"done in {time.time() - t0:.0f} s: {len(rows)} candidates, "
        f"{counts.get('extended', 0)} extended; "
        + ", ".join(f"{k} {n}" for k, n in sorted(counts.items()) if k != "extended"))
    return dict(counts)


def resolve_mutation_rate(indir: str, prefixes: Sequence[str], given: Optional[float]) -> float:
    if given is not None:
        return given
    for p in prefixes:
        try:
            with open(os.path.join(indir, p + ".detect.json")) as fh:
                return float(json.load(fh)["settings"]["mutation_rate"])
        except (OSError, KeyError, ValueError, TypeError):
            continue
    warn("no --mutation-rate and no <prefix>.detect.json to read one from; using 3e-8")
    return 3e-8


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="ltrquest-reboundary",
        description="Extend truncated LTR-RT calls to the ends their family's LTR model "
                    "supports. Rewrites the _clean_ depth tables and FASTAs in place and "
                    "writes <prefix>_reboundary.tsv.")
    ap.add_argument("--indir", default=".", help="directory holding the run (default: .)")
    ap.add_argument("--prefix", nargs="+", required=True,
                    help="genome prefix(es); all are pooled, since families span genomes")
    ap.add_argument("--genome", nargs="+",
                    help="original (unmasked) genome FASTA per prefix, in the same order")
    ap.add_argument("--method", choices=METHODS, default="consensus",
                    help="family model (default: consensus)")
    ap.add_argument("--references", choices=("modal", "tsd"), default="modal",
                    help="copies models are built from: modal-length (default) or TSD-bearing")
    ap.add_argument("--min-copies", type=int, default=10,
                    help="families with fewer copies are left alone (default 10)")
    ap.add_argument("--min-identity", type=float, default=0.8,
                    help="identity over the outer 30 bp a moved end needs (default 0.8)")
    ap.add_argument("--anchor-len", type=int, default=30,
                    help="bp of the model's end searched for past a large indel (default 30)")
    ap.add_argument("--no-anchor", action="store_true", help="do not search past large indels")
    ap.add_argument("--max-indel", type=int, default=5000,
                    help="how far past the call the anchor looks, bp (default 5000)")
    ap.add_argument("--credit", default="200",
                    help="bits of family evidence given to Kmer2LTR: a number, or 'model' "
                         "(default 200)")
    ap.add_argument("--max-ratio", type=float, default=1.15,
                    help="family QC: model length / modal called LTR length ceiling (default 1.15)")
    ap.add_argument("--qc-min-n", type=int, default=5,
                    help="family QC: candidates needed to test TSD enrichment (default 5)")
    ap.add_argument("--untested", choices=("accept", "skip"), default="accept",
                    help="families with too few candidates to test (default: accept)")
    ap.add_argument("--mutation-rate", type=float, default=None,
                    help="per site per year, for k2p_time (default: the run's own, else 3e-8)")
    ap.add_argument("--tools-dir", default=os.environ.get("LTRQUEST_TOOLS_DIR", "ltrquest_tools"),
                    help="Kmer2LTR checkout location, cloned there if absent "
                         "(default: $LTRQUEST_TOOLS_DIR or ./ltrquest_tools)")
    ap.add_argument("--mafft", default="mafft", help="mafft executable (default: mafft on PATH)")
    ap.add_argument("-t", "--threads", type=int, default=8, help="worker processes (default 8)")
    ap.add_argument("--dump-proposals", default=None,
                    help="also write every proposal to this TSV (benchmarks)")
    ap.add_argument("--cache", default=None,
                    help="pickle the family phase here; reuse it when the settings allow")
    ap.add_argument("-v", "--verbose", action="store_true", help="one line per family")
    return ap


def settings_from(args) -> Settings:
    if args.credit != "model":
        try:
            float(args.credit)
        except ValueError:
            raise SystemExit(f"reboundary: --credit must be a number or 'model' "
                             f"(got {args.credit!r})")
    return Settings(
        method=args.method, references=args.references, min_copies=args.min_copies,
        credit=args.credit, max_ratio=args.max_ratio, qc_min_n=args.qc_min_n,
        untested=args.untested,
        mutation_rate=resolve_mutation_rate(args.indir, args.prefix, args.mutation_rate),
        mafft=args.mafft, threads=args.threads,
        place=Params(min_identity=args.min_identity, anchor=not args.no_anchor,
                     anchor_len=args.anchor_len, max_indel=args.max_indel))


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    s = settings_from(args)
    if not args.genome:
        raise SystemExit("reboundary: --genome is required")
    run(args.indir, args.prefix, args.genome, s, args.tools_dir, args.dump_proposals,
        args.cache, args.verbose)
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run the tests**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT tests/test_reboundary.py -q`
Expected: 19 passed (6 unit, 13 integration).
Debugging aids if an integration test fails:
- read the sidecar: `column -t` on `<tmp>/syn_run*/syn_LTRs_reboundary.tsv` (the path is in the pytest failure output);
- `kmer2ltr_reverted` on `del_left` means credit 200 was not enough for this synthetic pair: re-run that one case in a Python shell with `rb.arbitrate` and `credit` 1000 to confirm it is the credit, then report — do **not** raise the default here; Task 10 tunes credit on real data;
- `family_untested`/`family_qc_failed` for the fixture family means the TSD probe missed the planted TSDs: check `tsd_at` on `fx.kind("del_left").true_start/true_end`.

Run the whole non-slow suite: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT -m "not slow" -q 2>&1 | tail -1`
Expected: no failures.

- [ ] **Step 5: Commit**

```bash
git add src/ltrquest/reboundary.py tests/test_reboundary.py
git commit -m "reboundary: pooled driver with family QC, conflicts and Kmer2LTR arbitration

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: GFF3 follows renamed elements (`--reboundary-map`)

**Files:**
- Modify: `src/ltrquest/gff3.py` (`strand_provenance`, `build_element_blocks`, `convert`, `main`)
- Test: `tests/test_gff3.py`

**Interfaces:**
- Consumes: `reboundary_io.Rebound`, `reboundary_io.read_map`, `SIDECAR_COLUMNS`.
- Produces: `gff3.convert(..., reboundary_map: Optional[str] = None)`; CLI flag `--reboundary-map PATH`; attributes `boundary_source=family_model` and `boundary_shift=<ext5>,<ext3>` on re-bounded elements only.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_gff3.py`:

```python
def _renamed_run(tmp_path, with_map: bool):
    from ltrquest.kmer2ltr import COLUMNS
    from ltrquest.reboundary_io import SIDECAR_COLUMNS
    tail = ["strand", "family", "domains", "nest_status"]
    old, new = "chr1:150-3000#LTR/Gypsy/Tekay", "chr1:100-3000#LTR/Gypsy/Tekay"
    vals = {c: "NA" for c in COLUMNS}
    vals.update(seq_id=new, seq_len="2901", status="pass", ltr5_start="1", ltr5_end="450",
                ltr3_start="2551", ltr3_end="2901", orientation="+", tsd="ACGTA",
                tsd_offset="0,0")
    row = [vals[c] for c in COLUMNS] + ["+", "merged_fam00001", ".", "."]
    (tmp_path / "p_depth0_clean_ltr.tsv").write_text(
        "#" + "\t".join(COLUMNS + tail) + "\n" + "\t".join(row) + "\n")
    cluster = tmp_path / "merged_all_ltr.consensus_id0.75_cluster.tsv"
    cluster.write_text(f"{old}\t{old}\n")
    kw = dict(consensus_cluster=str(cluster), family_prefix="merged")
    if with_map:
        side = {c: "." for c in SIDECAR_COLUMNS}
        side.update(old_seq_id=old, new_seq_id=new, decision="extended", ext5="50", ext3="0")
        path = tmp_path / "p_reboundary.tsv"
        path.write_text("#" + "\t".join(SIDECAR_COLUMNS) + "\n"
                        + "\t".join(side[c] for c in SIDECAR_COLUMNS) + "\n")
        kw["reboundary_map"] = str(path)
    assert gff3_module.convert("p", str(tmp_path), **kw) == 0
    text = (tmp_path / "p_all_depth_LTR_cleaned.gff3").read_text()
    return [l for l in text.splitlines() if f"\t{gff3_module.LTR_TYPE}\t" in l][0]


def test_reboundary_map_keeps_family_attributes_and_marks_the_shift(tmp_path):
    line = _renamed_run(tmp_path, with_map=True)
    assert "\t100\t3000\t" in line
    assert "family_size=1" in line
    assert "boundary_source=family_model" in line and "boundary_shift=50,0" in line


def test_without_the_map_a_renamed_element_loses_its_cluster_attributes(tmp_path):
    line = _renamed_run(tmp_path, with_map=False)
    assert "family_size" not in line and "boundary_source" not in line
```

`tests/test_gff3.py` already imports the module as `gff3_module`; the helper uses that name.

- [ ] **Step 2: Run to verify failure**

Run: `PYT tests/test_gff3.py -q -k reboundary_map`
Expected: FAIL — `TypeError: convert() got an unexpected keyword argument 'reboundary_map'`.

- [ ] **Step 3: Implement**

In `src/ltrquest/gff3.py`:

1. Imports: add `from .reboundary_io import Rebound, read_map`.

2. `strand_provenance` — add a trailing parameter and translate keys around the cascade:

```python
def strand_provenance(prefix: str, indir: str, tables, families,
                      verbose: bool = False,
                      recovered_strands: Optional[str] = None,
                      alias: Optional[Dict[str, str]] = None
                      ) -> Dict[str, Tuple[str, str]]:
```

and, inside, immediately after `elements = collect_elements(load_unannotated(t.path) for t in tables)`:

```python
    # Re-bounded elements are named by their new span, but every tier below is
    # keyed on the span the element was detected with.
    if alias:
        elements = {alias.get(k, k): v for k, v in elements.items()}
```

and replace the final `return {...}` with:

```python
    out = {key: (strand[key], source[key]) for key in elements
           if key in strand and key in source}
    if alias:
        back = {old: new for new, old in alias.items()}
        out = {back.get(k, k): v for k, v in out.items()}
    return out
```

3. `build_element_blocks` — add a trailing parameter `rebound: Optional[Dict[str, Rebound]] = None`; at the top of the body:

```python
    rebound = rebound or {}
    alias_key = {k: element_key(r.old_name) for k, r in rebound.items()}
    alias_name = {r.new_name: r.old_name for r in rebound.values()}
```

replace `family = lookup_family(name, key, family_by_name, family_by_coord)` with:

```python
            family = lookup_family(alias_name.get(name, name), alias_key.get(key, key),
                                   family_by_name, family_by_coord)
            rb = rebound.get(key)
```

and in the `render_attributes([...])` list insert, right after the `("tsd_offset", ...)` pair:

```python
                ("boundary_source", "family_model" if rb else ""),
                ("boundary_shift", list_value(f"{rb.ext5},{rb.ext3}") if rb else ""),
```

4. `convert` — add a trailing parameter `reboundary_map: Optional[str] = None` and replace the two calls:

```python
    rebound = read_map(reboundary_map) if reboundary_map else {}
    alias = {k: element_key(r.old_name) for k, r in rebound.items()}
    provenance = strand_provenance(prefix, indir, tables, families, verbose,
                                   recovered_strands, alias)
    element_blocks, skipped = build_element_blocks(prefix, tables, ranker,
                                                   provenance, families, verbose, rebound)
```

5. `main` — add the flag next to `--recovered-strands`:

```python
    parser.add_argument("--reboundary-map", default=None,
                        help="<prefix>_reboundary.tsv from ltrquest-reboundary: lets family "
                             "and strand lookups find elements it renamed, and marks them "
                             "boundary_source=family_model. Passed by path, never globbed.")
```

and pass `reboundary_map=args.reboundary_map` in the `convert(...)` call.

`element_key` is already imported in gff3.py from annotate (it is used in `build_element_blocks`); if not, add it to the `from .annotate import (...)` list.

- [ ] **Step 4: Run the tests**

Run: `PYT tests/test_gff3.py -q`
Expected: all pass, including the two new tests.
Run: `PYT -m "not slow" -q 2>&1 | tail -1` — no failures.

- [ ] **Step 5: Commit**

```bash
git add src/ltrquest/gff3.py tests/test_gff3.py
git commit -m "gff3: follow re-bounded elements through --reboundary-map; mark them

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: Post-hoc mode — backup once, re-run from the originals, regenerate GFF3 and plots, restore

**Files:**
- Modify: `src/ltrquest/reboundary.py` (new functions + `build_parser`/`main`)
- Test: `tests/test_reboundary.py`

**Interfaces:**
- Consumes: Task 6 `run`, Task 7 `gff3.convert(..., reboundary_map=...)`.
- Produces: `BACKUP_SUFFIX = "_pre_reboundary"`; `prepare_posthoc(indir, prefix) -> str`; `restore(indir, prefix) -> None`; `regenerate(indir, prefix, genome, plots) -> None`; CLI flags `--posthoc`, `--no-plots`, `--restore`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_reboundary.py`:

```python
def test_posthoc_backs_up_once_reruns_from_originals_and_restores(tmp_path, k2l_api, mafft):
    fx = build(tmp_path / "syn_posthoc")
    originals = {p.name: p.read_bytes() for p in fx.indir.glob(f"{fx.prefix}_depth*_clean_ltr.*")}
    (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").write_text("##gff-version 3\n#old\n")
    (fx.indir / f"{fx.prefix}_plots").mkdir()
    (fx.indir / f"{fx.prefix}_plots" / "x.pdf").write_text("old plot")
    base = ["--indir", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
            "-t", "2", "--mafft", mafft, "--tools-dir", TOOLS, "--posthoc", "--no-plots"]

    assert rb.main(base) == 0
    backup = fx.indir / f"{fx.prefix}{rb.BACKUP_SUFFIX}"
    assert {p.name for p in backup.iterdir()} >= set(originals) | {
        f"{fx.prefix}_all_depth_LTR_cleaned.gff3", f"{fx.prefix}_plots"}
    gff = (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").read_text()
    assert "boundary_source=family_model" in gff and "#old" not in gff
    assert not (fx.indir / f"{fx.prefix}_plots").exists()

    assert rb.main(base + ["--no-anchor"]) == 0      # starts again from the backup
    e = fx.kind("del_left")
    assert sidecar(fx)[e.key]["decision"] == "extended"

    assert rb.main(["--indir", str(fx.indir), "--prefix", fx.prefix, "--restore"]) == 0
    assert not backup.exists() and not (fx.indir / f"{fx.prefix}_reboundary.tsv").exists()
    for name, data in originals.items():
        assert (fx.indir / name).read_bytes() == data, name
    assert "#old" in (fx.indir / f"{fx.prefix}_all_depth_LTR_cleaned.gff3").read_text()
    assert (fx.indir / f"{fx.prefix}_plots" / "x.pdf").read_text() == "old plot"


def test_restore_without_a_backup_says_so(tmp_path):
    with pytest.raises(SystemExit, match="nothing to restore"):
        rb.main(["--indir", str(tmp_path), "--prefix", "p", "--restore"])


def test_an_interrupted_backup_stops_the_run(tmp_path):
    (tmp_path / f"p{rb.BACKUP_SUFFIX}.partial").mkdir()
    with pytest.raises(SystemExit, match="interrupted"):
        rb.prepare_posthoc(str(tmp_path), "p")
```

- [ ] **Step 2: Run to verify failure**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT tests/test_reboundary.py -q -k "posthoc or restore or interrupted"`
Expected: FAIL — `AttributeError: module 'ltrquest.reboundary' has no attribute 'BACKUP_SUFFIX'` (and argparse rejecting `--posthoc`).

- [ ] **Step 3: Implement**

In `src/ltrquest/reboundary.py` add `import glob` and `import subprocess` to the imports, and add after `resolve_mutation_rate`:

```python
BACKUP_SUFFIX = "_pre_reboundary"


def _clean_files(directory: str, prefix: str) -> List[str]:
    return sorted(glob.glob(os.path.join(directory, f"{prefix}_depth*_clean_ltr.tsv"))
                  + glob.glob(os.path.join(directory, f"{prefix}_depth*_clean_ltr.fa")))


def _derived(indir: str, prefix: str) -> List[str]:
    names = (f"{prefix}_all_depth_LTR_cleaned.gff3", f"{prefix}_all_depth_protein_LTR_cleaned.gff3",
             f"{prefix}_plots", prefix + SIDECAR_SUFFIX)
    return [os.path.join(indir, n) for n in names if os.path.lexists(os.path.join(indir, n))]


def _remove(path: str) -> None:
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    elif os.path.lexists(path):
        os.remove(path)


def prepare_posthoc(indir: str, prefix: str) -> str:
    """Move the originals into <prefix>_pre_reboundary/ once; put fresh copies back to work on.

    If the backup already exists this run starts again from it, so post-hoc
    runs with different settings never stack on each other.
    """
    bdir = os.path.join(indir, prefix + BACKUP_SUFFIX)
    if os.path.isdir(bdir + ".partial"):
        raise SystemExit(f"reboundary: {bdir}.partial exists: an earlier backup was interrupted. "
                         f"Move its files back into {indir}, delete it, then re-run.")
    if not os.path.isdir(bdir):
        tmp = bdir + ".partial"
        os.makedirs(tmp)
        for path in _clean_files(indir, prefix) + _derived(indir, prefix):
            shutil.move(path, os.path.join(tmp, os.path.basename(path)))
        os.rename(tmp, bdir)
        log(f"{prefix}: originals moved to {os.path.basename(bdir)}/")
    else:
        for path in _clean_files(indir, prefix) + _derived(indir, prefix):
            _remove(path)
        log(f"{prefix}: starting again from {os.path.basename(bdir)}/")
    originals = _clean_files(bdir, prefix)
    if not originals:
        raise SystemExit(f"reboundary: {bdir} holds no {prefix}_depth<N>_clean_ltr tables")
    for path in originals:
        shutil.copy2(path, os.path.join(indir, os.path.basename(path)))
    return bdir


def restore(indir: str, prefix: str) -> None:
    bdir = os.path.join(indir, prefix + BACKUP_SUFFIX)
    if not os.path.isdir(bdir):
        raise SystemExit(f"reboundary: nothing to restore: no {bdir}")
    for path in _clean_files(indir, prefix) + _derived(indir, prefix):
        _remove(path)
    for name in os.listdir(bdir):
        shutil.move(os.path.join(bdir, name), os.path.join(indir, name))
    os.rmdir(bdir)
    log(f"{prefix}: originals restored from {os.path.basename(bdir)}/")


def regenerate(indir: str, prefix: str, genome: str, plots: bool) -> None:
    """Rewrite the GFF3s (through the key map) and, unless told not to, the plots."""
    from . import gff3
    cons = sorted(glob.glob(os.path.join(indir, "*_all_ltr.consensus_id*_cluster.tsv")))
    if len(cons) > 1:
        raise SystemExit(f"reboundary: {len(cons)} consensus cluster tables in {indir}; "
                         f"expected at most one")
    fam = {}
    if cons:
        fam = dict(consensus_cluster=cons[0],
                   family_prefix=os.path.basename(cons[0]).split("_all_ltr.consensus_id")[0])
    rec = os.path.join(indir, prefix + "_strand_recovery.tsv")
    side = os.path.join(indir, prefix + SIDECAR_SUFFIX)
    if gff3.convert(prefix, indir, genome, recovered_strands=rec if os.path.isfile(rec) else None,
                    reboundary_map=side if os.path.isfile(side) else None, **fam) != 0:
        raise SystemExit(f"reboundary: GFF3 regeneration failed for {prefix}")
    if plots:
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts", "plots.sh")
        env = dict(os.environ, LTRQUEST_PYTHON=sys.executable)
        done = subprocess.run(["bash", script, "--prefix", prefix, "--genome", genome,
                               "--indir", indir], env=env)
        if done.returncode != 0:
            warn(f"plotting reported failures for {prefix}; the tables and GFF3 are unaffected")
```

In `build_parser()`, before `return ap`, add:

```python
    ap.add_argument("--posthoc", action="store_true",
                    help="update a finished run in place: originals go to "
                         "<prefix>_pre_reboundary/ (once; later runs start from it), then the "
                         "GFF3s and plots are rewritten")
    ap.add_argument("--no-plots", action="store_true", help="with --posthoc: skip the plots")
    ap.add_argument("--restore", action="store_true",
                    help="put <prefix>_pre_reboundary/ back in place and stop")
```

Replace `main` with:

```python
def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.restore:
        for p in args.prefix:
            restore(args.indir, p)
        return 0
    s = settings_from(args)
    if not args.genome:
        raise SystemExit("reboundary: --genome is required")
    if len(args.genome) != len(args.prefix):
        raise SystemExit("reboundary: --prefix and --genome must pair up one to one")
    if args.posthoc:
        for p in args.prefix:
            prepare_posthoc(args.indir, p)
    run(args.indir, args.prefix, args.genome, s, args.tools_dir, args.dump_proposals,
        args.cache, args.verbose)
    if args.posthoc:
        for p, gpath in zip(args.prefix, args.genome):
            regenerate(args.indir, p, gpath, plots=not args.no_plots)
    return 0
```

- [ ] **Step 4: Run the tests**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT tests/test_reboundary.py -q`
Expected: 22 passed.
Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT -m "not slow" -q 2>&1 | tail -1` — no failures.

- [ ] **Step 5: Commit**

```bash
git add src/ltrquest/reboundary.py tests/test_reboundary.py
git commit -m "reboundary: post-hoc mode with one-time backup, regeneration and restore

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---
### Task 9: Benchmark harness (B1 planted obstacles, B2 real run, B3 false changes) and the first results

**Files:**
- Create: `benchmarks/reboundary/bench.py`
- Create: `benchmarks/reboundary/README.md`
- Test: `tests/test_bench_smoke.py`
- Outputs (not in the repo): `/data2/chris/poa_LTR/ltrquest_run/reboundary_bench/`

**Interfaces:**
- Consumes: `reboundary._init/build_models/credit_for/arbitrate/run/Settings`, `ltr_model.Genomes/Member/ratio_ok/tsd_enrichment`, `ltr_place.Params/propose`, `reboundary_io.load_clean_tables/members_from`.
- Produces: `bench.py b1|b2|collect`; each run writes `<out>/summary.json` (keys used later: B1 `exact_or_1, exact, over, wrong, false_change, n_truth, table`; B2 `extended, tsd_gain, tsd_null, cv_tsd, cv_null, cv_n, b3_real_changed, b3_real_tsd_kept, ext_median, ext_q90, dk2p_median, dk2p_q90, wall_s, peak_rss_gb`) and a TSV (`b1.tsv` or `proposals.tsv` + the sidecars under `<out>/run/`).

- [ ] **Step 1: Write the harness**

Create `benchmarks/reboundary/bench.py`:

```python
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
    return dict(exact_or_1=frac(obst["exact"] + obst["within1"], n_obst),
                exact=frac(obst["exact"], n_obst), over=frac(obst["over"], n_obst),
                wrong=frac(obst["wrong"], n_obst), unchanged=frac(obst["unchanged"], n_obst),
                false_change=frac(ctrl["false_change"], sum(ctrl.values())), table=table)


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
        rows.append(dict(kind=kind, size=size, delta=delta, end=end, copy=copy, family=t.family,
                         model=status,
                         proposal="none" if prop is None else ("ok" if prop.gate_ok else "gated"),
                         verdict=v.status if v else ".", truth_left=tr[0], truth_right=tr[1],
                         called_left=m.start, called_right=m.end, final_left=left,
                         final_right=right, outcome=outcome(m, tr, left, right, kind)))
    write_tsv(os.path.join(args.out, "b1.tsv"), rows)
    summary = summarize_b1(rows)
    summary.update(settings=asdict(s), n_truth=len(truth), n_contigs=len(synth),
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
    return dict(cv_n=n, cv_tsd=frac(hits, n), cv_null=frac(null, n), cv_p0=round(p0, 4))


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
    summary = dict(settings=asdict(s), candidates=len(side), extended=len(ext), counts=counts,
                   n_fresh=len(fresh),
                   tsd_gain=frac(sum(found(r["tsd_new"]) for r in fresh), len(fresh)),
                   tsd_null=frac(sum(found(r["tsd_null"]) for r in fresh), len(fresh)),
                   b3_real_changed=frac(len(had), tsd_calls),
                   b3_real_tsd_kept=frac(sum(found(r["tsd_new"]) for r in had), len(had)),
                   ext_median=q(exts, 0.5), ext_q90=q(exts, 0.9),
                   dk2p_median=q(dk, 0.5), dk2p_q90=q(dk, 0.9),
                   wall_s=round(wall), peak_rss_gb=round(rss_kb / 1e6, 2))
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
            out[k] = dict(candidates=len(rs), extended=len(e),
                          tsd_gain=frac(sum(found(r["tsd_new"]) for r in f), len(f)))
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
```

`rb.place_params(s)` does not exist yet: add it to `src/ltrquest/reboundary.py` now (Task 10 extends it), right after `build_models`:

```python
def place_params(s: Settings) -> Params:
    """Placement parameters for a run (the method may adjust them; see Task 10)."""
    return s.place
```

and use it in `family_job`: replace `p = propose(m, kept or models, _G, s.place)` with `p = propose(m, kept or models, _G, place_params(s))`.

- [ ] **Step 2: Write the smoke test**

Create `tests/test_bench_smoke.py`:

```python
"""The benchmark harness runs end to end on the synthetic run."""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import pytest

from reboundary_fixtures import build

BENCH = Path(__file__).resolve().parents[1] / "benchmarks" / "reboundary" / "bench.py"


def load():
    spec = importlib.util.spec_from_file_location("bench", BENCH)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["bench"] = mod          # multiprocessing pickles bench._b1_family by name
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.skipif(not BENCH.is_file(), reason="benchmarks are not in this checkout")
def test_b1_b2_and_collect_run_on_the_synthetic_run(tmp_path, k2l_api, mafft, monkeypatch,
                                                    capsys):
    bench = load()
    fx = build(tmp_path / "syn")
    common = ["--run", str(fx.indir), "--prefix", fx.prefix, "--genome", str(fx.genome),
              "--tools-dir", os.environ.get("LTRQUEST_TOOLS_DIR", "."), "--threads", "2",
              "--mafft", mafft]
    monkeypatch.setattr(sys, "argv", ["bench.py", "b1", *common, "--out",
                                      str(tmp_path / "t" / "b1_x"), "--n", "3"])
    bench.main()
    s1 = json.loads((tmp_path / "t" / "b1_x" / "summary.json").read_text())
    assert s1["n_truth"] == 3 and s1["exact_or_1"] is not None and s1["false_change"] == 0.0
    monkeypatch.setattr(sys, "argv", ["bench.py", "b2", *common, "--out",
                                      str(tmp_path / "t" / "b2_x")])
    bench.main()
    s2 = json.loads((tmp_path / "t" / "b2_x" / "summary.json").read_text())
    assert s2["extended"] == 4 and s2["tsd_gain"] == 1.0
    monkeypatch.setattr(sys, "argv", ["bench.py", "collect", "--out", str(tmp_path / "t")])
    bench.main()
    out = capsys.readouterr().out
    assert "| x |" in out and "B1 planted obstacles" in out and "B2 real run" in out
```

- [ ] **Step 3: Run the smoke test**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT tests/test_bench_smoke.py tests/test_reboundary.py -q`
Expected: all pass. If `s1["false_change"]` is not 0.0, open `<tmp>/t/b1_x/b1.tsv` and look at the `control` rows: a changed control means an intact, correctly-called element was moved — a real bug in placement or arbitration, fix it before continuing.

- [ ] **Step 4: Write the benchmark README**

Create `benchmarks/reboundary/README.md`:

````markdown
# Re-boundarying benchmarks

What `ltrquest-reboundary` is measured on, how to reproduce it, and the numbers
behind its defaults. Design: `docs/superpowers/specs/2026-09-22-reboundary-design.md` §9.

## The three benchmarks

- **B1 — planted obstacles (ground truth).** Elements whose called ends carry a
  Kmer2LTR TSD (independent evidence the call is right) are copied onto synthetic
  contigs with 3 kb of their real flanks. One obstacle is planted in one LTR,
  `delta` bp from an outer end: a 30 bp patch at 30% divergence, a deletion
  (10/50/200 bp), an insertion (50/300/1,000/5,000 bp of random sequence, or 300/1,000 bp
  cut from another family's element — a TE-like insertion), or nothing. The call is then
  cut at the obstacle exactly as LTRharvest / LTR_FINDER stop. Each truth element is
  left out of its own family's model (leave-one-out). Scored: recovered outer end
  exact / within 1 / within 5 bp, over-extended, wrong, left unchanged.
- **B2 — a real run.** The run's clean tables are copied and re-bounded. Measured:
  elements extended; among those with no TSD when called, the TSD rate at the new
  ends against a null with the right flank displaced by 1 kb; the same with family
  QC chosen on one half of each family's candidates and scored on the other (`cv`),
  so TSD-based QC cannot inflate its own validation; K2P shifts; time and memory.
- **B3 — false changes.** B1: untouched truth elements (no obstacle, no truncation)
  that get altered. B2: calls that already had a TSD and were altered, and how many
  of those keep a TSD.

## Decision rule (spec §9)

Maximise B1 exact-or-±1, subject to B3 false change ≤ 0.5% and B2 cross-validated
TSD gain ≥ 74.2% (the spike's figure). Ties go to the simpler, then the faster setting.

## Reproduce

```bash
cd /data2/chris/poa_LTR/ltrquest_run          # any finished LTRquest run directory
export PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH
export PYTHONPATH=/data2/chris/poa_LTR/LTRquest/src
B=/data2/chris/poa_LTR/LTRquest/benchmarks/reboundary/bench.py
P="annua.nuclear_LTRs chaixii.nuclear_LTRs infirma.nuclear_LTRs supina.nuclear_LTRs"
G="annua.nuclear.fa chaixii.nuclear.fa infirma.nuclear.fa supina.nuclear.fa"
T=/data2/chris/poa_LTR/ltrquest_tools          # Kmer2LTR checkout at aa25f46
python $B b1 --run . --prefix $P --genome $G --tools-dir $T --out reboundary_bench/b1_NAME --n 1000
python $B b2 --run . --prefix $P --genome $G --tools-dir $T --out reboundary_bench/b2_NAME \
    --cache reboundary_bench/cache_NAME.pkl
python $B collect --out reboundary_bench
```

## Results

Data: 4 *Poa* genomes (annua, chaixii, infirma, supina), LTRquest run of 2026-09-19
(212,360 elements). Kmer2LTR aa25f46. Tables are pasted from `bench.py collect`.
````

- [ ] **Step 5: Commit the harness**

```bash
git add benchmarks/reboundary/bench.py benchmarks/reboundary/README.md tests/test_bench_smoke.py src/ltrquest/reboundary.py
git commit -m "benchmarks: reboundary harness (planted obstacles, real run, false changes)

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 6: Launch the first B1 and B2 on the Poa run (consensus defaults), detached**

```bash
cd /data2/chris/poa_LTR/ltrquest_run && mkdir -p reboundary_bench
cat > reboundary_bench/run_first.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cd /data2/chris/poa_LTR/ltrquest_run
export PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH
export PYTHONPATH=/data2/chris/poa_LTR/LTRquest/src
B=/data2/chris/poa_LTR/LTRquest/benchmarks/reboundary/bench.py
P=(annua.nuclear_LTRs chaixii.nuclear_LTRs infirma.nuclear_LTRs supina.nuclear_LTRs)
G=(annua.nuclear.fa chaixii.nuclear.fa infirma.nuclear.fa supina.nuclear.fa)
T=/data2/chris/poa_LTR/ltrquest_tools
/usr/bin/time -v python $B b1 --run . --prefix "${P[@]}" --genome "${G[@]}" --tools-dir $T \
  --out reboundary_bench/first/b1_consensus --n 1000 --threads 32
/usr/bin/time -v python $B b2 --run . --prefix "${P[@]}" --genome "${G[@]}" --tools-dir $T \
  --out reboundary_bench/first/b2_consensus --cache reboundary_bench/first/cache_consensus.pkl \
  --threads 32
python $B collect --out reboundary_bench/first
echo FIRST_DONE
EOF
chmod +x reboundary_bench/run_first.sh
setsid nohup reboundary_bench/run_first.sh > reboundary_bench/run_first.log 2>&1 < /dev/null &
sleep 3; ps -o pid,ppid,sid,etime,cmd -u "$USER" | grep -E 'run_first|bench.py' | grep -v grep
```
Expected: the script and a `python .../bench.py b1` process listed; the script's SID equals its own PID (its own session). Wait for `FIRST_DONE` in `reboundary_bench/run_first.log` (use a Monitor with an `until grep -q -e FIRST_DONE -e Traceback -e Error ...` loop; do not poll with sleep). Expected duration: B1 10–30 min, B2 30–90 min.

- [ ] **Step 7: Record the first results**

Append the `collect` output from the log to the `## Results` section of `benchmarks/reboundary/README.md` under a heading `### First pass: consensus model, spike-derived defaults (2026-MM-DD)` (today's date), followed by 2–4 sentences stating: B1 exact-or-±1, which obstacle kinds fail (read the per-kind table printed by `b1`), B3 false change, B2 cv TSD gain vs null, wall time and peak RSS against the ≤1 h / ≤16 GB target.

Write `/data2/chris/poa_LTR/ltrquest_run/reboundary_bench/memo.md` (terse key:value + the exact `run_first.sh` commands, per Global Constraints), then:

```bash
cd /data2/chris/poa_LTR/LTRquest
git add benchmarks/reboundary/README.md
git commit -m "benchmarks: first reboundary results on the Poa run

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

**Checkpoint — report to the user** before Task 10: the first-pass table, and anything that falls short of the decision rule (e.g. B3 > 0.5%, cv TSD gain < 74%). If B3 false change is above 0.5% stop and debug (systematic-debugging) before tuning anything.

---

### Task 10: `subfamily` and `nearest` methods, `references=tsd`, and the tuning sweeps that set the defaults

**Files:**
- Modify: `src/ltrquest/ltr_model.py` (add `subfamily_models`)
- Modify: `src/ltrquest/ltr_place.py` (add `combine` to `Params`, median combining in `propose`, `nearest_templates`)
- Modify: `src/ltrquest/reboundary.py` (`METHODS`, `build_models`, `place_params`, CLI `--subfamily-jaccard`)
- Create: `benchmarks/reboundary/tune.sh`
- Modify: `benchmarks/reboundary/README.md`
- Test: `tests/test_ltr_model.py`, `tests/test_ltr_place.py`, `tests/test_reboundary.py`

**Interfaces:**
- Consumes: Tasks 3–9.
- Produces: `ltr_model.subfamily_models(family, refs, g, modal, jaccard_min, mafft="mafft") -> List[Model]`; `ltr_place.Params.combine: str = "best"` (`"best"` | `"median"`); `ltr_place.nearest_templates(consensus, refs, g, tol=2) -> List[Model]`; `reboundary.METHODS = ("consensus", "subfamily", "nearest")`; `reboundary.place_params(s)` returns `replace(s.place, combine="median")` for `nearest`, else `s.place`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_ltr_model.py`:

```python
def test_subfamily_models_on_a_one_type_family_give_the_family_ltr(fx, g, fam, mafft):
    refs, modal = lm.select_references(fam, "modal", FAMILY, n_young=150, n_random=150)
    models = lm.subfamily_models(FAMILY, refs, g, modal, 0.5, mafft)
    assert models and all(abs(len(mo.seq) - 400) <= 2 for mo in models)
    assert models[0].seq[:20] == fx.ltr[:20]
```

Append to `tests/test_ltr_place.py`:

```python
def test_nearest_templates_are_references_whose_ends_the_consensus_confirms(fx, g, truth_model):
    refs, _ = lm.select_references([m for m in members(fx) if m.family == FAMILY], "modal", FAMILY)
    tpl = lp.nearest_templates(truth_model, refs, g)
    assert len(tpl) >= 10 and all(t.model_id.startswith(f"{FAMILY}:tpl:") for t in tpl)
    assert all(len(t.seq) == 400 for t in tpl)          # the truncated ref (ins_left) is not one


def test_median_combining_agrees_with_best_on_a_clean_case(fx, g, truth_model):
    refs, _ = lm.select_references([m for m in members(fx) if m.family == FAMILY], "modal", FAMILY)
    tpl = lp.nearest_templates(truth_model, refs, g)
    m, e = member(fx, "del_left")
    p = lp.propose(m, tpl, g, lp.Params(combine="median"))
    assert p.gate_ok and p.left == e.true_start
```

Append to `tests/test_reboundary.py`:

```python
@pytest.mark.parametrize("method", ["consensus", "subfamily", "nearest"])
def test_every_method_recovers_the_deletion(tmp_path, mafft, method):
    from ltrquest.ltr_model import Genomes
    from ltrquest.ltr_place import propose
    from reboundary_fixtures import FAMILY, members
    fx = build(tmp_path / "syn_methods")
    g = Genomes({fx.prefix: str(fx.genome)})
    fam = [m for m in members(fx) if m.family == FAMILY]
    s = rb.Settings(method=method, mafft=mafft)
    models, _, status = rb.build_models(FAMILY, fam, g, s)
    assert status == "ok" and models
    e = fx.kind("del_left")
    m = [x for x in fam if x.name == e.name][0]
    p = propose(m, models, g, rb.place_params(s))
    assert p.gate_ok and p.left == e.true_start


def test_nearest_uses_median_combining():
    assert rb.place_params(rb.Settings(method="nearest")).combine == "median"
    assert rb.place_params(rb.Settings()).combine == "best"
```

- [ ] **Step 2: Run to verify failure**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT tests/test_ltr_model.py tests/test_ltr_place.py tests/test_reboundary.py -q -k "subfamily or nearest or median or every_method"`
Expected: FAIL — missing `subfamily_models`, `nearest_templates`, `Params.combine`.

- [ ] **Step 3: Implement**

Append to `src/ltrquest/ltr_model.py`:

```python
def subfamily_models(family: str, refs: Sequence[Member], g: Genomes, modal: Optional[float],
                     jaccard_min: float, mafft: str = "mafft") -> List[Model]:
    """One consensus per cluster of similar reference LTRs (>= MIN_REFS copies), plus the rest.

    Clusters grow greedily around centroids taken youngest-first (select_references
    returns the youngest first), on canonical k-mer Jaccard of each reference's 5'
    LTR. With no cluster big enough, this is the family consensus.
    """
    profiles = [canonical_kmers(oriented_ltrs(m, g)[0]) for m in refs]
    centroids: List[int] = []
    clusters: List[List[int]] = []
    for i, prof in enumerate(profiles):
        for c, members in zip(centroids, clusters):
            if jaccard(prof, profiles[c]) >= jaccard_min:
                members.append(i)
                break
        else:
            centroids.append(i)
            clusters.append([i])
    models: List[Model] = []
    big = [cl for cl in clusters if len(cl) >= MIN_REFS]
    for n, cl in enumerate(big, start=1):
        mo = consensus_model(family, [refs[i] for i in cl[:80]], g, modal, mafft,
                             f"{family}:sub{n}")
        if mo is not None:
            models.append(mo)
    rest = [refs[i] for cl in clusters if len(cl) < MIN_REFS for i in cl]
    if big and len(rest) >= MIN_REFS:
        mo = consensus_model(family, rest[:80], g, modal, mafft, f"{family}:rest")
        if mo is not None:
            models.append(mo)
    if not models:
        mo = consensus_model(family, list(refs)[:80], g, modal, mafft)
        if mo is not None:
            models.append(mo)
    return models
```

In `src/ltrquest/ltr_place.py`:

1. Add `import dataclasses` and `import statistics`, and extend the model import to `from .ltr_model import (GAP_EXTEND, GAP_OPEN, Genomes, Member, Model, canonical_kmers, glocal, jaccard, matrix, oriented_ltrs, rc)`.
2. Add a field to `Params` after `n_models`: `combine: str = "best"      # best | median of the top n_models placements`.
3. In `propose`, replace the model/orientation loop and the `c, model, o, fwd = best` line with:

```python
    cands = []
    for model in candidate_models(m, models, g, p.n_models):
        top = None
        for o in ((m.strand,) if m.stranded else ("+", "-")):
            fwd = model.seq if o == "+" else rc(model.seq)
            c = core(m, fwd, g)
            if c is not None and (top is None or c.score > top[0].score):
                top = (c, model, o, fwd)
        if top is not None:
            cands.append(top)
    if not cands:
        return None
    cands.sort(key=lambda t: -t[0].score)
    c, model, o, fwd = cands[0]
    if p.combine == "median" and len(cands) >= 3:
        c = dataclasses.replace(c, left=int(statistics.median(x[0].left for x in cands)),
                                right=int(statistics.median(x[0].right for x in cands)))
```

4. Append:

```python
def nearest_templates(consensus: Model, refs: Sequence[Member], g: Genomes,
                      tol: int = 2) -> List[Model]:
    """References whose called outer ends the family consensus confirms, as models of their own.

    A template's termini are then checked by its family rather than taken on trust.
    Falls back to the consensus itself when no reference qualifies.
    """
    out: List[Model] = []
    for m in refs:
        fwd = consensus.seq if m.strand != "-" else rc(consensus.seq)
        c = core(m, fwd, g)
        if c is None or abs(c.left - m.start) > tol or abs(c.right - m.end) > tol:
            continue
        five, _ = oriented_ltrs(m, g)
        out.append(Model(f"{consensus.family}:tpl:{m.key}", consensus.family, five,
                         consensus.modal_len, 1))
    return out or [consensus]
```

In `src/ltrquest/reboundary.py`:

1. `METHODS = ("consensus", "subfamily", "nearest")`.
2. Imports: add `subfamily_models` to the `.ltr_model` import, `nearest_templates` to the `.ltr_place` import, and `replace` to the `dataclasses` import.
3. Replace `build_models` and `place_params` with:

```python
def build_models(family: str, members: Sequence[Member], g: Genomes, s: Settings,
                 exclude: FrozenSet[str] = frozenset()) -> Tuple[List[Model], Optional[float], str]:
    many = s.method == "subfamily"
    refs, modal = select_references(members, s.references, family, exclude,
                                    n_young=150 if many else 40, n_random=150 if many else 40)
    if not refs:
        return [], modal, "too_few_references"
    if s.method == "subfamily":
        models = subfamily_models(family, refs, g, modal, s.subfamily_jaccard, s.mafft)
        return models, modal, ("ok" if models else "no_ltr_span")
    model = consensus_model(family, refs, g, modal, s.mafft)
    if model is None:
        return [], modal, "no_ltr_span"
    if s.method == "nearest":
        return nearest_templates(model, refs, g), modal, "ok"
    return [model], modal, "ok"


def place_params(s: Settings) -> Params:
    """`nearest` takes the median of its top templates' placements (spec 5.1)."""
    return replace(s.place, combine="median") if s.method == "nearest" else s.place
```

4. In `build_parser`, add after `--references`:

```python
    ap.add_argument("--subfamily-jaccard", type=float, default=0.5,
                    help="with --method subfamily: LTR 11-mer Jaccard to join a cluster "
                         "(default 0.5)")
```

and in `settings_from` pass `subfamily_jaccard=args.subfamily_jaccard`.

- [ ] **Step 4: Run the tests**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT -m "not slow" -q 2>&1 | tail -1`
Expected: no failures.

- [ ] **Step 5: Commit the methods**

```bash
git add src/ltrquest/ltr_model.py src/ltrquest/ltr_place.py src/ltrquest/reboundary.py tests/
git commit -m "reboundary: subfamily and nearest-template models; median combining

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 6: Write the tuning script**

Create `benchmarks/reboundary/tune.sh` (executable):

```bash
#!/usr/bin/env bash
# Tuning sweeps for ltrquest-reboundary (design spec section 9, "Tuning order").
# Run from the LTRquest run directory. Each stage writes reboundary_bench/tune<STAGE>/.
#   tune.sh 1                         method x references            (B1 n=500 + B2)
#   METHOD=.. REFS=.. tune.sh 2       min_identity x anchor_len x max_indel  (B1 n=500)
#   KNOBS="--set ..." tune.sh 3       credit                         (B2 cached + B1)
#   KNOBS="--set ..." tune.sh 4a      max_ratio x min_copies          (B2 cached)
#   KNOBS="--set ..." tune.sh 4b      qc_min_n x untested             (B2 cached)
# Finished cells are skipped, so a stage can be re-run after an interruption.
set -euo pipefail
STAGE="$1"
REPO=/data2/chris/poa_LTR/LTRquest
export PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH
export PYTHONPATH=$REPO/src
P=(annua.nuclear_LTRs chaixii.nuclear_LTRs infirma.nuclear_LTRs supina.nuclear_LTRs)
G=(annua.nuclear.fa chaixii.nuclear.fa infirma.nuclear.fa supina.nuclear.fa)
OUT=reboundary_bench/tune${STAGE}
mkdir -p "$OUT"
bench() {
  python "$REPO/benchmarks/reboundary/bench.py" "$@" --run . --prefix "${P[@]}" \
    --genome "${G[@]}" --tools-dir /data2/chris/poa_LTR/ltrquest_tools --threads "${THREADS:-32}"
}
done_or() { [[ -s "$1/summary.json" ]]; }
case "$STAGE" in
  1)
    for m in ${METHODS:-consensus subfamily nearest}; do
      for r in modal tsd; do
        n="${m}_${r}"
        done_or "$OUT/b1_$n" || bench b1 --out "$OUT/b1_$n" --n 500 --set method=$m --set references=$r
        done_or "$OUT/b2_$n" || bench b2 --out "$OUT/b2_$n" --cache "$OUT/cache_$n.pkl" \
          --set method=$m --set references=$r
      done
    done ;;
  2)
    : "${METHOD:?set METHOD to the stage-1 winner}" "${REFS:?set REFS to the stage-1 winner}"
    for id in 0.7 0.8; do for al in 30 50 80; do for mi in 1000 5000 20000; do
      n="id${id}_al${al}_mi${mi}"
      done_or "$OUT/b1_$n" || bench b1 --out "$OUT/b1_$n" --n 500 --set method=$METHOD \
        --set references=$REFS --set min_identity=$id --set anchor_len=$al --set max_indel=$mi
    done; done; done ;;
  3)
    : "${KNOBS:?set KNOBS to the stage-2 winner as --set pairs}"
    for c in 0 50 200 1000 model; do
      done_or "$OUT/b2_credit_$c" || bench b2 --out "$OUT/b2_credit_$c" --cache "$OUT/cache.pkl" \
        $KNOBS --set credit=$c
      done_or "$OUT/b1_credit_$c" || bench b1 --out "$OUT/b1_credit_$c" --n 500 $KNOBS --set credit=$c
    done ;;
  4a)
    : "${KNOBS:?set KNOBS to the stage-3 winner as --set pairs}"
    # min_copies 5 first: its cached family phase covers the 10 and 20 runs as well.
    for mr in 1.10 1.15 1.25; do for mc in 5 10 20; do
      n="ratio${mr}_min${mc}"
      done_or "$OUT/b2_$n" || bench b2 --out "$OUT/b2_$n" --cache "$OUT/cache_ratio${mr}.pkl" \
        $KNOBS --set max_ratio=$mr --set min_copies=$mc
    done; done ;;
  4b)
    : "${KNOBS:?set KNOBS to the stage-4a winner as --set pairs}"
    for qn in 3 5 10; do for u in accept skip; do
      n="n${qn}_${u}"
      done_or "$OUT/b2_$n" || bench b2 --out "$OUT/b2_$n" --cache "$OUT/cache.pkl" \
        $KNOBS --set qc_min_n=$qn --set untested=$u
    done; done ;;
  *) echo "usage: tune.sh 1|2|3|4a|4b" >&2; exit 2 ;;
esac
python "$REPO/benchmarks/reboundary/bench.py" collect --out "$OUT"
echo "STAGE_${STAGE}_DONE"
```

```bash
chmod +x benchmarks/reboundary/tune.sh
git add benchmarks/reboundary/tune.sh
git commit -m "benchmarks: tuning sweeps for reboundary

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 7: Stage 1 — method × references**

Three methods run in parallel, 32 threads each (the server has 256 cores; check `uptime` first and use `THREADS=16` if the load average is above ~150):

```bash
cd /data2/chris/poa_LTR/ltrquest_run
for m in consensus subfamily nearest; do
  METHODS=$m setsid nohup /data2/chris/poa_LTR/LTRquest/benchmarks/reboundary/tune.sh 1 \
    > reboundary_bench/tune1_$m.log 2>&1 < /dev/null &
done
sleep 3; ps -o pid,ppid,sid,etime,cmd -u "$USER" | grep tune.sh | grep -v grep
```
Wait for `STAGE_1_DONE` in all three logs (Monitor with an until-loop that also matches `Traceback|Error`). Then:
`python /data2/chris/poa_LTR/LTRquest/benchmarks/reboundary/bench.py collect --out reboundary_bench/tune1`

Choose the winner by the decision rule: among settings with B1 `B3 false` ≤ 0.005 and B2 `cv TSD` ≥ 0.742, the highest B1 `exact/±1`; within 0.01 of the best, prefer consensus < subfamily < nearest and modal < tsd. If no setting reaches cv TSD ≥ 0.742, take the highest cv TSD among those with B3 ≤ 0.005 and say so in the README.

- [ ] **Step 8: Stage 2 — placement knobs (B1 only)**

```bash
METHOD=<winner method> REFS=<winner refs> setsid nohup \
  /data2/chris/poa_LTR/LTRquest/benchmarks/reboundary/tune.sh 2 > reboundary_bench/tune2.log 2>&1 < /dev/null &
```
Winner: highest B1 `exact/±1` with `B3 false` ≤ 0.005; ties → min_identity 0.8, anchor_len 30, max_indel 5000 (the current defaults).

- [ ] **Step 9: Stage 3 — credit**

```bash
KNOBS="--set method=<m> --set references=<r> --set min_identity=<id> --set anchor_len=<al> --set max_indel=<mi>" \
  setsid nohup /data2/chris/poa_LTR/LTRquest/benchmarks/reboundary/tune.sh 3 > reboundary_bench/tune3.log 2>&1 < /dev/null &
```
Winner: highest B1 `exact/±1` subject to B1 `over + wrong` ≤ 0.02, `B3 false` ≤ 0.005 and B2 `cv TSD` ≥ 0.742.

- [ ] **Step 10: Stage 4 — family QC and family size**

```bash
KNOBS="<stage-3 KNOBS> --set credit=<c>" setsid nohup \
  /data2/chris/poa_LTR/LTRquest/benchmarks/reboundary/tune.sh 4a > reboundary_bench/tune4a.log 2>&1 < /dev/null &
# after STAGE_4a_DONE, with the 4a winner added to KNOBS:
KNOBS="<stage-3 KNOBS> --set credit=<c> --set max_ratio=<r> --set min_copies=<n>" setsid nohup \
  /data2/chris/poa_LTR/LTRquest/benchmarks/reboundary/tune.sh 4b > reboundary_bench/tune4b.log 2>&1 < /dev/null &
```
Winner of each: most B2 `extended` subject to `cv TSD` ≥ 0.742 (else the highest cv TSD).

- [ ] **Step 11: Set the defaults and record the evidence**

Change the defaults to the winners in exactly these places (keep them identical):
- `src/ltrquest/reboundary.py` `Settings` fields (`method`, `references`, `min_copies`, `credit`, `max_ratio`, `qc_min_n`, `untested`, `subfamily_jaccard`);
- `src/ltrquest/ltr_place.py` `Params` fields (`min_identity`, `anchor_len`, `max_indel`);
- the matching `default=` and help-text defaults in `reboundary.build_parser()`.

Add to `benchmarks/reboundary/README.md` under `## Results`: one `### Stage N` subsection per stage with its `collect` table and one sentence naming the winner and why; then a `### Chosen defaults` table (knob | value | stage | deciding numbers). Update `reboundary_bench/memo.md` with the tune commands.

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT -m "not slow" -q 2>&1 | tail -1` — no failures (a default change can move a fixture test; if one fails, check whether the new default is legitimately stricter on the synthetic case before touching the test, and report it).

```bash
git add src/ltrquest/reboundary.py src/ltrquest/ltr_place.py benchmarks/reboundary/README.md
git commit -m "reboundary: defaults chosen by the benchmark sweeps

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

**Checkpoint — report to the user**: the chosen defaults with their numbers, and anything the benchmark showed the module cannot do (e.g. kb-scale insertions rejected by the 0.65 length-ratio gate). Do not relax LTRquest's own length gate without the user's decision.

---

### Task 11: Driver integration (`ltrquest --reboundary`) and packaging

**Files:**
- Modify: `src/ltrquest/scripts/ltrquest.sh`
- Modify: `pyproject.toml`
- Modify: `tests/test_packaging.py`
- Create: `tests/test_driver_reboundary.py`

**Interfaces:**
- Consumes: CLI `ltrquest-reboundary` / `python -m ltrquest.reboundary` (Tasks 6, 8), `gff3 --reboundary-map` (Task 7).
- Produces: driver flag `--reboundary`; shell function `run_reboundary_stage`; shell function `reorient_recovered`; console script `ltrquest-reboundary`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_driver_reboundary.py`:

```python
"""The driver's --reboundary wiring: flag, early check, stage order, file hygiene."""

from __future__ import annotations

import subprocess


def run(*args, **kw):
    return subprocess.run(args, capture_output=True, text=True, **kw)


def test_help_documents_the_flag(driver):
    result = run("bash", str(driver), "--help")
    assert "--reboundary" in result.stdout + result.stderr


def test_the_flag_parses(driver):
    assert run("bash", str(driver), "--reboundary", "--help").returncode == 0


def test_a_missing_mafft_fails_before_any_work(driver, toy_genome):
    result = run("bash", str(driver), "--genome", str(toy_genome), "--reboundary",
                 env={"PATH": "/usr/bin:/bin", "LTRQUEST_PYTHON": "/nonexistent/python"})
    assert result.returncode != 0 and "mafft" in result.stderr


def test_it_is_not_a_detection_setting(driver):
    text = driver.read_text()
    block = text[text.index("detection_settings() {"):text.index("plan_reuse() {")]
    assert "REBOUNDARY" not in block


def test_carry_forward_and_promotion_drop_its_files(driver):
    text = driver.read_text()
    assert "*_reboundary.*|*_pre_reboundary) continue;;" in text
    assert 'rm -f "${dst}/${p}_reboundary.tsv"' in text
    assert 'rm -rf "${dst}/${p}_pre_reboundary"' in text


def test_stage_order(driver):
    text = driver.read_text()
    stage = text[text.index("run_annotation_stage() {"):text.index("run_reboundary_stage() {")]
    off_gff3 = stage.index('"${GFF3[@]}" --prefix "$p" --indir . \\\n      --genome "${p}.input_genome.fa" "${fam_opts[@]}" "${rec_opts[@]}"\n')
    on_skip = stage.index('if [[ "$REBOUNDARY" == true ]]; then')
    assert on_skip < off_gff3                            # with the flag, loop 1 stops before gff3
    assert stage.index("run_reboundary_stage\n") > off_gff3   # then pooled, then gff3 with the map
    assert '--reboundary-map "${p}_reboundary.tsv"' in stage
```

In `tests/test_packaging.py` add `"reboundary"` to `STAGES` and `"ltrquest-reboundary"` to `CONSOLE_SCRIPTS`.

- [ ] **Step 2: Run to verify failure**

Run: `PYT tests/test_driver_reboundary.py -q`
Expected: FAIL (`--reboundary` unknown; strings absent).

- [ ] **Step 3: Implement the driver changes**

In `src/ltrquest/scripts/ltrquest.sh`:

(a) Defaults — after `STRAND_RECOVERY_PPT=false` add:

```bash

# Family-guided re-boundarying of truncated calls (opt-in post-processing)
REBOUNDARY=false
```

(b) `usage()` — after the two `--strand-recovery-ppt` help lines add:

```
  --reboundary          Extend LTR-RT calls that stop short of their true ends
                        (an indel or a mutation-dense patch near an LTR end stops
                        LTRharvest/LTR_FINDER early). Each family's LTR model,
                        built from its full-length copies across all genomes,
                        proposes the ends; Kmer2LTR re-scores every changed
                        element. Extends only, never trims. Rewrites the _clean_
                        tables and FASTAs and writes <out_prefix>_reboundary.tsv.
                        Off by default. Needs mafft. For a finished run use
                        ltrquest-reboundary --posthoc instead.
```

(c) Argument parsing — after `--strand-recovery-ppt) STRAND_RECOVERY_PPT=true; shift;;` add:

```bash
    --reboundary) REBOUNDARY=true; shift;;
```

(d) Validation — after the strand-recovery validation `fi` (the one closing `elif [[ "$STRAND_RECOVERY_PPT" == true ]]; then die ...`) add:

```bash

# --reboundary builds family LTR models with mafft; fail now, not after detection.
if [[ "$REBOUNDARY" == true ]]; then
  command -v mafft >/dev/null 2>&1 \
    || die "--reboundary needs mafft on PATH to build family LTR models; it is in environment.yml"
fi
```

(e) Stage commands — after `RECORD=(    "$PY" -m ltrquest.record    )` add:

```bash
REBOUND=(   "$PY" -m ltrquest.reboundary )
```

(f) `promote_fp_outputs` — after `rm -f "${dst}/${p}_strand_recovery.tsv"` add:

```bash
    # Same for re-boundarying: a sidecar or post-hoc backup left by an earlier
    # run describes that run's elements, not these.
    rm -f "${dst}/${p}_reboundary.tsv"
    rm -rf "${dst}/${p}_pre_reboundary"
```

(g) `carry_forward_genome` — after the `*_strand_recovery.*) continue;;` case line add:

```bash
      # Re-boundarying's sidecar and a post-hoc backup belong to one attempt too.
      *_reboundary.*|*_pre_reboundary) continue;;
```

(h) `run_annotation_stage` — replace everything from the line `    echo ""` that precedes `echo "Writing LTR-RT GFF3 for ${p}..."` down to the function's closing `}` with:

```bash
    if [[ "$REBOUNDARY" == true ]]; then
      # Re-orient now, not after the GFF3: re-boundarying renames the elements
      # the apply phase looks up by name, and the GFF3 reads neither the FASTAs
      # nor the orientation column, so it can come after.
      if (( ${#rec_opts[@]} > 0 )); then
        reorient_recovered "$p"
      fi
      continue
    fi

    echo ""
    echo "============================================================"
    echo "Writing LTR-RT GFF3 for ${p}..."
    set -x
    "${GFF3[@]}" --prefix "$p" --indir . \
      --genome "${p}.input_genome.fa" "${fam_opts[@]}" "${rec_opts[@]}"
    set +x

    # Phase 2: bring the depth FASTAs into line with the strand column the
    # annotator just wrote, so a recovered minus element is stored in coding
    # sense exactly like a TEsorter2-called one.
    if (( ${#rec_opts[@]} > 0 )); then
      reorient_recovered "$p"
    fi
  done

  [[ "$REBOUNDARY" == true ]] || return 0
  run_reboundary_stage
  local -a rb_opts=()
  for i in "${!OUT_PREFIXES[@]}"; do
    p="${OUT_PREFIXES[$i]}"
    depth_tables_for "$p" tsv annot_tsvs
    (( ${#annot_tsvs[@]} > 0 )) || continue
    rec_opts=()
    if [[ -n "$STRAND_RECOVERY" && -s "${p}_strand_recovery.tsv" ]]; then
      rec_opts=( --recovered-strands "${p}_strand_recovery.tsv" )
    fi
    rb_opts=()
    if [[ -s "${p}_reboundary.tsv" ]]; then
      rb_opts=( --reboundary-map "${p}_reboundary.tsv" )
    fi
    echo ""
    echo "============================================================"
    echo "Writing LTR-RT GFF3 for ${p}..."
    set -x
    "${GFF3[@]}" --prefix "$p" --indir . \
      --genome "${p}.input_genome.fa" "${fam_opts[@]}" "${rec_opts[@]}" "${rb_opts[@]}"
    set +x
  done
}

reorient_recovered() {
  local p="$1"
  echo ""
  echo "============================================================"
  echo "Re-orienting ${p} depth FASTAs to the recovered strand..."
  set -x
  "${RECOVER[@]}" --phase apply --prefix "$p" --indir . || {
    set +x
    echo "WARNING: could not re-orient the depth FASTAs for ${p}; the tables" >&2
    echo "and GFF3s still carry the recovered strand." >&2
  }
  set +x
}

# Pooled over every genome: a family's LTR model is built from all its copies,
# whichever genome they sit in. After the annotator (it needs the family and
# strand columns) and before the GFF3 (which must learn the names it changes).
run_reboundary_stage() {
  local i p
  local -a prefixes=() genomes=() cleaned=()
  for i in "${!OUT_PREFIXES[@]}"; do
    p="${OUT_PREFIXES[$i]}"
    shopt -s nullglob
    cleaned=( "${p}"_depth*_clean_ltr.tsv )
    shopt -u nullglob
    (( ${#cleaned[@]} > 0 )) || continue
    prefixes+=( "$p" )
    # The ORIGINAL genome, as for strand recovery: from FP attempt 2 onward the
    # staged input is the hard-masked FASTA.
    genomes+=( "${abs_genomes[$i]:-${p}.input_genome.fa}" )
  done
  if (( ${#prefixes[@]} == 0 )); then
    echo "WARNING: no _clean_ depth tables; skipping --reboundary." >&2
    return 0
  fi
  resolve_merged_tools_dir
  ensure_kmer2ltr_dir
  echo ""
  echo "============================================================"
  echo "Re-boundarying truncated LTR-RT calls (${#prefixes[@]} genome(s), pooled)..."
  set -x
  if "${REBOUND[@]}" --indir . --prefix "${prefixes[@]}" --genome "${genomes[@]}" \
       --threads "$THREADS" --mutation-rate "$MUTATION_RATE" --tools-dir "$TOOLS_DIR"; then
    set +x
  else
    set +x
    echo "WARNING: re-boundarying failed; the calls are kept as detected." >&2
    for p in "${prefixes[@]}"; do rm -f "${p}_reboundary.tsv"; done
  fi
}
```

The removed apply block (the old `if (( ${#rec_opts[@]} > 0 )); then echo ... "${RECOVER[@]}" --phase apply ... fi` inside the loop) is now `reorient_recovered`, which prints the same lines and runs the same command, so with `--reboundary` off the stage behaves exactly as before.

(i) `pyproject.toml` — under `[project.scripts]`, after the `ltrquest-recover-strand` line add:

```toml
ltrquest-reboundary = "ltrquest.reboundary:main"
```

- [ ] **Step 4: Run the tests**

Run: `PYT tests/test_driver_reboundary.py tests/test_packaging.py -q`
Expected: `test_driver_reboundary.py` all pass; `bash -n` of the driver passes (`test_scripts_parse`). The two parametrized `TestConsoleScripts` cases for `ltrquest-reboundary` fail until the package is re-installed, because they look for the installed console script — this is expected here. **Ask the user** to run `mamba run -n ltrquest pip install --no-deps /data2/chris/poa_LTR/LTRquest` (their install, per their rules), then re-run `PYT tests/test_packaging.py -q` and confirm all pass. Do not install it yourself.

Run: `PYT -m "not slow" -q 2>&1 | tail -1` — no other failures.

- [ ] **Step 5: Commit**

```bash
git add src/ltrquest/scripts/ltrquest.sh pyproject.toml tests/test_packaging.py tests/test_driver_reboundary.py
git commit -m "ltrquest.sh: --reboundary stage between annotation and GFF3; packaging

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 12: Nextflow parity

**Files:**
- Create: `modules/local/ltrquest/reboundary/main.nf`, `modules/local/ltrquest/reboundary/environment.yml`
- Modify: `modules/local/ltrquest/gff3/main.nf`, `workflows/ltrquest.nf`, `nextflow.config`, `nextflow_schema.json`, `conf/modules.config`

**Interfaces:**
- Consumes: `ltrquest-reboundary` CLI, `ltrquest-gff3 --reboundary-map`.
- Produces: `params.reboundary` (default false); process `LTRQUEST_REBOUNDARY`; `LTRQUEST_GFF3` input gains `path(reboundary_map)`.

- [ ] **Step 1: Write the module**

`modules/local/ltrquest/reboundary/environment.yml`:

```yaml
name: ltrquest_reboundary
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - mafft=7.526
  - pyfaidx
  - numpy
  # `parasail` on Bioconda is the C library; `import parasail` needs this one.
  - parasail-python
  - pywfa
```

`modules/local/ltrquest/reboundary/main.nf`:

```groovy
// Family-guided re-boundarying. One pooled task: a family's LTR model is built
// from its copies in every genome, so the samples cannot be split here the way
// the per-genome stages are. It rewrites the _clean_ tables and FASTAs and
// writes one <prefix>_reboundary.tsv per sample, which LTRQUEST_GFF3 reads to
// follow the elements it renamed.

process LTRQUEST_REBOUNDARY {
    tag "pooled"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/cwb14/ltrquest:1.0.1'

    input:
    val(ids)
    path(tables, stageAs: 'in/*')
    path(fastas, stageAs: 'in/*')
    path(genomes, stageAs: 'genomes/g??/*')

    output:
    path("*_depth*_clean_ltr.tsv"), emit: tsv
    path("*_depth*_clean_ltr.fa") , emit: fasta
    path("*_reboundary.tsv")      , emit: sidecar
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_list = (genomes instanceof List ? genomes : [genomes]).join(' ')
    // Rewritten in place, so copied out of the staging directory first: editing
    // through Nextflow's input symlinks would corrupt the upstream task's outputs.
    """
    cp in/*_clean_ltr.tsv in/*_clean_ltr.fa .

    ltrquest-reboundary \\
        --indir . \\
        --prefix ${ids.join(' ')} \\
        --genome ${genome_list} \\
        --threads ${task.cpus} \\
        --mutation-rate ${params.mutation_rate} \\
        --tools-dir \${LTRQUEST_TOOLS_DIR:-/opt/ltrquest/tools} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: \$(python -c 'import ltrquest; print(ltrquest.__version__)')
        mafft: \$(mafft --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    cp in/*_clean_ltr.tsv in/*_clean_ltr.fa .
    for id in ${ids.join(' ')}; do
      printf '#old_seq_id\\tnew_seq_id\\tdecision\\n' > \${id}_reboundary.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ltrquest: 1.0.1
    END_VERSIONS
    """
}
```

- [ ] **Step 2: Extend the GFF3 module**

In `modules/local/ltrquest/gff3/main.nf`:
- input tuple becomes `tuple val(meta), path(tables), path(workdirs), path(genome), path(consensus_cluster), path(recovered_strands), path(reboundary_map)`;
- in `script:` add `def rbmap = reboundary_map ? "--reboundary-map ${reboundary_map}" : ''` next to `def recovered = ...`, and add the line `        ${rbmap} \\` after `        ${recovered} \\` in the command.

- [ ] **Step 3: Wire the workflow**

In `workflows/ltrquest.nf`:
1. Add `include { LTRQUEST_REBOUNDARY   } from '../modules/local/ltrquest/reboundary/main'` after the GFF3 include.
2. Immediately after the `if (params.strand_recovery) { ... } else { ... }` block that sets `ch_tables` / `ch_fastas`, insert:

```groovy
    // Re-boundarying, when asked for: one pooled task over every sample, between
    // the annotated tables and the GFF3. Off, an empty map flows in its place.
    ch_meta_by_id = ch_bundle.map { meta, _t, _f, _w, _g, _c -> [ meta.id, meta ] }
    if (params.reboundary) {
        ch_rb_in = ch_tables
            .join(ch_fastas)
            .join(ch_bundle.map { meta, _t, _f, _w, genome, _c -> [ meta, genome ] })
            .toSortedList { a, b -> a[0].id <=> b[0].id }
        LTRQUEST_REBOUNDARY(
            ch_rb_in.map { rows -> rows.collect { it[0].id } },
            ch_rb_in.map { rows -> rows.collect { it[1] }.flatten() },
            ch_rb_in.map { rows -> rows.collect { it[2] }.flatten() },
            ch_rb_in.map { rows -> rows.collect { it[3] } }
        )
        ch_versions = ch_versions.mix(LTRQUEST_REBOUNDARY.out.versions)
        ch_tables = LTRQUEST_REBOUNDARY.out.tsv.flatten()
            .map { f -> [ f.name.replaceFirst(/_depth\d+_clean_ltr\.tsv$/, ''), f ] }
            .groupTuple()
            .join(ch_meta_by_id)
            .map { _id, files, meta -> [ meta, files ] }
        ch_fastas = LTRQUEST_REBOUNDARY.out.fasta.flatten()
            .map { f -> [ f.name.replaceFirst(/_depth\d+_clean_ltr\.fa$/, ''), f ] }
            .groupTuple()
            .join(ch_meta_by_id)
            .map { _id, files, meta -> [ meta, files ] }
        ch_rbmap = LTRQUEST_REBOUNDARY.out.sidecar.flatten()
            .map { f -> [ f.name.replaceFirst(/_reboundary\.tsv$/, ''), f ] }
            .join(ch_meta_by_id)
            .map { _id, f, meta -> [ meta, f ] }
    } else {
        ch_rbmap = ch_bundle.map { meta, _t, _f, _w, _g, _c -> [ meta, [] ] }
    }
```

3. In the `LTRQUEST_GFF3(...)` call, change `}.join(ch_recovery)` to `}.join(ch_recovery).join(ch_rbmap)`.
4. At the top of the workflow's `main:` (next to the strand-recovery parameter checks) nothing is needed: `reboundary` is a boolean.

- [ ] **Step 4: Parameters, schema, publishing**

`nextflow.config` — after `strand_recovery_ppt        = false` add:

```groovy

    // Family-guided re-boundarying of truncated calls (opt-in)
    reboundary                 = false
```

`nextflow_schema.json` — after the `strand_recovery_ppt` property object add (mind the comma):

```json
                "reboundary": {
                    "type": "boolean",
                    "default": false,
                    "description": "Extend LTR-RT calls that stop short of their true ends, using each family's LTR model.",
                    "help_text": "Off unless set. LTRharvest and LTR_FINDER stop extending an LTR pair at the first indel or mutation-dense patch near an LTR end; each family's LTR model, built from its full-length copies across all samples, proposes the true ends and Kmer2LTR re-scores every changed element. Extends only, never trims. Writes <sample>_reboundary.tsv and marks changed elements boundary_source=family_model in the GFF3.",
                    "fa_icon": "fas fa-arrows-alt-h"
                },
```

`conf/modules.config`:
- `LTRQUEST_ANNOTATE`: `enabled: !params.strand_recovery` → `enabled: !params.strand_recovery && !params.reboundary`.
- `LTRQUEST_RECOVERSTRAND_APPLY`: add `enabled: !params.reboundary,` as the first publishDir key.
- Add:

```groovy
    withName: 'LTRQUEST_REBOUNDARY' {
        publishDir = [
            // The re-bounded clean tables replace the annotator's in the results,
            // plus each sample's sidecar.
            path: { "${params.outdir}/${file.name.replaceFirst(/(_depth\d+_clean_ltr\.(tsv|fa)|_reboundary\.tsv)$/, '')}" },
            mode: params.publish_dir_mode,
            pattern: '{*_depth*_clean_ltr.tsv,*_reboundary.tsv}',
            saveAs: { fn -> fn }
        ]
    }
```

If the closure over `file` is rejected by the installed Nextflow version, use `path: { "${params.outdir}/reboundary" }` instead and note it in the docs (Task 13): the pooled task has no single `meta.id`.

- [ ] **Step 5: Stub-run the DAG, both ways**

```bash
cd /data2/chris/poa_LTR/LTRquest
NF=$(command -v nextflow || ls /home/chris/bin/mambaforge/envs/nextflow/bin/nextflow 2>/dev/null || true)
echo "${NF:-no nextflow found}"
```
If no nextflow is found: stop this step and tell the user Nextflow is needed to verify the wiring (they install it); continue with Task 13 meanwhile and come back.
Otherwise:

```bash
$NF run . -profile test -stub-run --outdir /tmp/claude-rb-nf-off
$NF run . -profile test -stub-run --outdir /tmp/claude-rb-nf-on --reboundary
```
Expected: both complete; the second lists `LTRQUEST_REBOUNDARY` once and `LTRQUEST_GFF3` once per sample. Remove the two outdirs and `work/` afterwards (`rm -rf /tmp/claude-rb-nf-* work .nextflow*` — only these paths).

- [ ] **Step 6: Commit**

```bash
git add modules/local/ltrquest/reboundary modules/local/ltrquest/gff3/main.nf workflows/ltrquest.nf nextflow.config nextflow_schema.json conf/modules.config
git commit -m "nextflow: pooled LTRQUEST_REBOUNDARY between annotation and GFF3

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 13: Documentation

**Files:**
- Modify: `docs/outputs.md`, `README.md`, `CHANGELOG.md`

- [ ] **Step 1: `docs/outputs.md`**

Insert a new section before `## 8. Plots (<prefix>_plots/)` and renumber `## 8. Plots` → `## 9. Plots` and `## 9. Benchmarks` → `## 10. Benchmarks`; then `grep -n "section 8\|section 9\|§8\|§9" docs/outputs.md README.md` and update any reference to the renumbered sections.

````markdown
## 8. Re-boundarying (`--reboundary`)

LTRharvest and LTR_FINDER extend an LTR pair outward from a seed and stop at the
first obstacle between the element's two LTRs — an indel or a mutation-dense patch
near an LTR end. Kmer2LTR can only trim. So such an element is called short: in a
family alignment it starts late and ends early, and it has no TSD at its called ends.

`--reboundary` (or `ltrquest-reboundary --posthoc` on a finished run) builds each
family's LTR model from its full-length copies, pooled over every genome, and places
it on every member. The model's ends come only from where the family's 5′ and 3′ LTR
copies stop agreeing — never from a TSD or TG..CA. Ends move outward only. A family
whose model is implausibly long, or whose proposed ends show no TSD enrichment over a
displaced-flank null, is left alone. Each widened element is re-scored by Kmer2LTR,
which settles the final ends and every Kmer2LTR column (LTR coordinates, divergence,
K2P, time, TSD, motif, CIGAR).

What changes: the `_clean_` depth tables and FASTAs, the GFF3s and the plots. The raw
`<prefix>_depth<N>_ltr.{tsv,fa}` tables keep the calls as detected. Family labels do
not change; the pooled `<run>_all_ltr.*` clustering files describe the calls before
re-boundarying.

### 8.1 `<prefix>_reboundary.tsv`

One row per element whose model placement moved an end outward, extended or not.

| Column | Meaning |
|---|---|
| `old_seq_id`, `new_seq_id` | the call before and after (`.` when not extended) |
| `family`, `method`, `model_id` | the family and the model that placed it |
| `decision`, `reason` | `extended`, or `rejected` with one of: `gate_identity`, `family_qc_failed`, `family_untested`, `host_exceeded`, `engulfs_element`, `overlaps_element`, `record_missing`, `kmer2ltr_not_pass`, `kmer2ltr_reverted`, `length_filter` |
| `ext5`, `ext3` | bp added at the element's biological 5′ / 3′ end |
| `end_source5`, `end_source3` | `core` (model aligned end to end), `anchor` (model's outer bases found past a large indel), `.` (did not move) |
| `id_outer5`, `id_outer3` | identity of the model's outer 30 bp at each end |
| `credit_bits`, `k2l_status` | evidence handed to Kmer2LTR and its status |
| `tsd_called`, `tsd_new`, `tsd_null` | TSD at the called ends, at the new ends, and with the right flank displaced 1 kb (the null) |
| `k2p_called`, `k2p_new` | LTR-pair K2P before and after |
| `obstacle5`, `obstacle3` | what sits between the two LTRs at the old end: `gap:<bp>`, `mm:<rate>`, `none`, `.` |

### 8.2 GFF3

Re-bounded elements carry `boundary_source=family_model` and
`boundary_shift=<ext5>,<ext3>`. The GFF3 writer reads the sidecar (`--reboundary-map`)
so family and strand-provenance attributes still find renamed elements.

### 8.3 Finished runs

```bash
ltrquest-reboundary --posthoc --indir RUN --prefix P1 P2 ... --genome G1 G2 ... --threads 32
ltrquest-reboundary --restore --indir RUN --prefix P1 P2 ...     # undo
```

`--posthoc` moves the originals (clean tables and FASTAs, GFF3s, plots) into
`<prefix>_pre_reboundary/` once, re-bounds, then rewrites the GFF3s and plots. A
second `--posthoc` run starts again from that backup, so trying other settings never
stacks on an earlier result. Benchmarks and the evidence behind every default:
`benchmarks/reboundary/README.md`.
````

- [ ] **Step 2: `README.md`**

Insert before `## Nextflow`:

````markdown
## Re-boundarying truncated calls

Some LTR-RT calls stop short of the element's real ends: an indel or a patch of
mutations near one LTR end makes LTRharvest and LTR_FINDER stop extending there. In a
family alignment those copies start late and end early. `--reboundary` fixes them
using each family's own full-length copies as the guide, then lets Kmer2LTR re-score
every changed element. It only ever extends a call, never trims one.

```bash
ltrquest --genome A.fa B.fa --proteins prot.fa --threads 32 --reboundary
```

Already have a finished run? Update it in place (originals are kept in
`<prefix>_pre_reboundary/`; `--restore` puts them back):

```bash
ltrquest-reboundary --posthoc --indir my_run --prefix A_LTRs B_LTRs --genome A.fa B.fa --threads 32
```

Every candidate, and why it was or was not extended, is in `<prefix>_reboundary.tsv`.
Details: [docs/outputs.md §8](docs/outputs.md).
````

- [ ] **Step 3: `CHANGELOG.md`**

Under `## [Unreleased]` → `### Added`, add as the first bullet:

```markdown
- `--reboundary` / `ltrquest-reboundary`: extend LTR-RT calls that stop short of
  their true ends. Each family's LTR model (termini from 5′/3′-LTR agreement across
  its full-length copies, pooled over genomes) proposes outward-only ends; family QC,
  conflict rules and Kmer2LTR re-scoring decide. Rewrites the `_clean_` tables and
  FASTAs, marks changed elements `boundary_source=family_model` in the GFF3, and
  writes `<prefix>_reboundary.tsv`. `--posthoc` updates a finished run in place with
  a backup. Defaults set by `benchmarks/reboundary/`.
```

- [ ] **Step 4: Commit**

```bash
git add docs/outputs.md README.md CHANGELOG.md
git commit -m "docs: re-boundarying

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 14: Final verification and review

- [ ] **Step 1: Full suite**

Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT -m "not slow" -q 2>&1 | tail -3`
Expected: all pass (the two console-script cases pass only after the user's reinstall, Task 11).
Run: `PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH PYT -m slow -q 2>&1 | tail -3` — report the result as is (the slow end-to-end tests need the full toolchain; a skip or environment failure there is reported, not "fixed").

- [ ] **Step 2: Spec coverage walk**

Open the spec and, section by section (§2 goals, §4 files, §5.1–5.6, §6, §7, §8, §9, §10), write one line each in the final report: where it is implemented (file:function) or why it differs (e.g. the `anchor_len` candidates, the module split into four files).

- [ ] **Step 3: Code review**

Use superpowers:requesting-code-review on the branch diff `main...reboundary`. Fix confirmed findings (with tests), commit.

---

### Task 15: Apply to the Poa run (post-hoc) — only with the user's go-ahead

This modifies the user's real results directory (with a backup). **Ask first**, showing what will move into `<prefix>_pre_reboundary/` (clean tables/FASTAs, both GFF3s, the plots directories) and that analyses reading `*_plots/*_combined_ltr.tsv` (e.g. `ltr_biology/`) will see re-bounded elements afterwards.

- [ ] **Step 1: Launch (after a yes), detached**

```bash
cd /data2/chris/poa_LTR/ltrquest_run
cat > reboundary_posthoc.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cd /data2/chris/poa_LTR/ltrquest_run
export PATH=/home/chris/bin/mambaforge/envs/ltrquest/bin:$PATH
export PYTHONPATH=/data2/chris/poa_LTR/LTRquest/src
/usr/bin/time -v python -m ltrquest.reboundary --posthoc --indir . \
  --prefix annua.nuclear_LTRs chaixii.nuclear_LTRs infirma.nuclear_LTRs supina.nuclear_LTRs \
  --genome annua.nuclear.fa chaixii.nuclear.fa infirma.nuclear.fa supina.nuclear.fa \
  --threads 32 --tools-dir /data2/chris/poa_LTR/ltrquest_tools
echo POSTHOC_DONE
EOF
chmod +x reboundary_posthoc.sh
setsid nohup ./reboundary_posthoc.sh > reboundary_posthoc.log 2>&1 < /dev/null &
sleep 3; ps -o pid,ppid,sid,etime,cmd -u "$USER" | grep -E 'reboundary_posthoc|ltrquest.reboundary' | grep -v grep
```
Wait for `POSTHOC_DONE` (Monitor, matching `POSTHOC_DONE|Traceback|Error`).

- [ ] **Step 2: Check the result**

```bash
cd /data2/chris/poa_LTR/ltrquest_run
for p in annua chaixii infirma supina; do
  s=${p}.nuclear_LTRs_reboundary.tsv
  printf '%s\t' "$p"; awk -F'\t' 'NR>1{d[$6"/"$7]++} END{for(k in d) printf "%s=%d ", k, d[k]; print ""}' "$s"
done
grep -c 'boundary_source=family_model' *.nuclear_LTRs_all_depth_LTR_cleaned.gff3
ls -d *_pre_reboundary
```
Expected: extended counts in the thousands (spike: ~11% of copies in QC-pass families); GFF3 counts equal the extended counts; four backup directories.

- [ ] **Step 3: Memo**

Write `/data2/chris/poa_LTR/ltrquest_run/memo_reboundary.md`: date | purpose | env (ltrquest env + `PYTHONPATH` to the branch, commit hash from `git -C /data2/chris/poa_LTR/LTRquest rev-parse --short HEAD`) | inputs | the exact `reboundary_posthoc.sh` | per-genome extended counts and top reject reasons | wall time and peak RSS from `/usr/bin/time -v` | how to undo (`--restore`) | gotchas | `claude --resume` id.

- [ ] **Step 4: Report** the per-genome counts, TSD gain at new ends, and where the backup and memo are.
