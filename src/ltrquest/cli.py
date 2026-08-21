"""Console entry point for the ``ltrquest`` pipeline driver.

The driver itself is ``scripts/ltrquest.sh``: the round loop is a long-running
sequence of external tools, and bash is the right language for that. This module
exists to make it a first-class installed command, and to hand it two things it
cannot work out for itself:

* ``LTRQUEST_PYTHON`` - the interpreter that owns this package. Without it the
  driver would fall back to ``python3`` on ``PATH``, which in a conda/venv layout
  is not necessarily the one ``ltrquest`` was installed into.
* the absolute path of the packaged script, wherever pip put it.

``os.execv`` replaces this process rather than spawning a child, so signals,
exit codes and job control behave exactly as if bash had been invoked directly.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parent / "scripts"


def _script(name: str) -> str:
    path = _SCRIPTS / name
    if not path.is_file():
        sys.exit(
            f"ltrquest: packaged script is missing: {path}\n"
            "The installation is incomplete; try reinstalling the package."
        )
    return str(path)


def _exec(name: str, argv: list[str]) -> "None":
    bash = shutil.which("bash")
    if bash is None:
        sys.exit("ltrquest: bash is required but was not found on PATH.")

    env = dict(os.environ)
    env.setdefault("LTRQUEST_PYTHON", sys.executable)
    env["LTRQUEST_ARGV0"] = Path(sys.argv[0]).name or name
    os.execve(bash, [bash, _script(name), *argv], env)


def main() -> "None":
    """Run the full pipeline: ``ltrquest --genome ... [options]``."""
    _exec("ltrquest.sh", sys.argv[1:])


def plots() -> "None":
    """Run the plotting stage alone: ``ltrquest-plots --prefix ... --genome ...``."""
    _exec("plots.sh", sys.argv[1:])


if __name__ == "__main__":  # pragma: no cover - exercised via the console script
    main()
