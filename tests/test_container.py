"""The container entry point, the launcher, and the wiring that reaches them.

Two bugs made the 1.0.1 image awkward to use, and both are cheap to guard:

* `./ltrquest.sif --help` answered `FATAL: "--help": executable file not found
  in $PATH`. Apptainer builds a SIF's runscript from the image's ENTRYPOINT and
  CMD, and an image with only a CMD gets a runscript that *replaces* the command
  with the user's arguments -- so the first flag was run as a program.
* Apptainer mounts your home and working directories and nothing else. A working
  directory reached through a symlink out of home (`/home/you/data ->
  /scratch/you`) is not followed from inside, so Apptainer falls back to $HOME
  and every relative path in the command silently points somewhere else.

Nothing here needs a container runtime: the entry point and the launcher are
plain shell, driven against stubs.
"""

from __future__ import annotations

import os
import stat
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
ENTRYPOINT = REPO_ROOT / "bin" / "entrypoint.sh"
LAUNCHER = REPO_ROOT / "bin" / "ltrquest-container"

SHELLS = ["/bin/sh", "/bin/bash"]


def _stub(path: Path, body: str) -> Path:
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


@pytest.fixture
def echo_bin(tmp_path: Path) -> Path:
    """A directory of stub commands that print how they were called."""
    binned = tmp_path / "bin"
    binned.mkdir()
    for name in ("ltrquest", "ltrquest-gff3", "ltrquest-plots"):
        _stub(binned / name, f'#!/bin/sh\necho "{name}: $*"\n')
    return binned


def run_entrypoint(shell: str, echo_bin: Path, *args: str) -> subprocess.CompletedProcess:
    env = dict(os.environ, PATH=f"{echo_bin}:{os.environ['PATH']}")
    return subprocess.run(
        [shell, str(ENTRYPOINT), *args],
        capture_output=True, text=True, env=env,
    )


class TestScriptsAreWellFormed:
    @pytest.mark.parametrize("script", [ENTRYPOINT, LAUNCHER], ids=lambda p: p.name)
    @pytest.mark.parametrize("shell", SHELLS)
    def test_parses(self, script: Path, shell: str):
        # Both are POSIX sh on purpose: the entry point runs in whatever the
        # image's /bin/sh is, and the launcher runs on whatever the cluster has.
        result = subprocess.run([shell, "-n", str(script)], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr

    @pytest.mark.parametrize("script", [ENTRYPOINT, LAUNCHER], ids=lambda p: p.name)
    def test_is_executable(self, script: Path):
        assert script.stat().st_mode & 0o111, f"{script.name} must ship executable"

    @pytest.mark.parametrize("script", [ENTRYPOINT, LAUNCHER], ids=lambda p: p.name)
    def test_starts_with_a_posix_shebang(self, script: Path):
        assert script.read_text().startswith("#!/bin/sh")


class TestEntrypointDispatch:
    """A first argument that names a command is a command; anything else is an
    ltrquest argument. That is strictly more permissive than the empty
    ENTRYPOINT it replaced, so no documented invocation stops working."""

    @pytest.mark.parametrize("shell", SHELLS)
    def test_bare_flags_go_to_ltrquest(self, shell, echo_bin):
        result = run_entrypoint(shell, echo_bin, "--genome", "x.fa", "--threads", "20")
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == "ltrquest: --genome x.fa --threads 20"

    @pytest.mark.parametrize("shell", SHELLS)
    def test_no_arguments_still_reaches_ltrquest(self, shell, echo_bin):
        # Apptainer supplies the image's CMD (`--help`) in this case; here the
        # point is only that an empty argument list is not an error.
        result = run_entrypoint(shell, echo_bin)
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == "ltrquest:"

    @pytest.mark.parametrize("shell", SHELLS)
    def test_an_explicit_command_still_works(self, shell, echo_bin):
        # `docker run IMG ltrquest --help`, as the README has always documented.
        result = run_entrypoint(shell, echo_bin, "ltrquest", "--help")
        assert result.stdout.strip() == "ltrquest: --help"

    @pytest.mark.parametrize("shell", SHELLS)
    def test_the_other_stages_are_reachable(self, shell, echo_bin):
        result = run_entrypoint(shell, echo_bin, "ltrquest-gff3", "--help")
        assert result.stdout.strip() == "ltrquest-gff3: --help"

    @pytest.mark.parametrize("shell", SHELLS)
    def test_an_arbitrary_command_runs(self, shell, echo_bin):
        # This is how Nextflow drives the image: `docker run IMG /bin/bash -c ...`
        result = run_entrypoint(shell, echo_bin, "/bin/sh", "-c", "echo nextflow-task")
        assert result.stdout.strip() == "nextflow-task"

    @pytest.mark.parametrize("shell", SHELLS)
    def test_arguments_keep_their_spaces(self, shell, echo_bin):
        result = run_entrypoint(shell, echo_bin, "--out_prefix", "two words")
        assert result.stdout.strip() == "ltrquest: --out_prefix two words"

    @pytest.mark.parametrize("shell", SHELLS)
    def test_an_unresolvable_first_argument_is_an_ltrquest_argument(self, shell, echo_bin):
        result = run_entrypoint(shell, echo_bin, "not-a-command", "--threads", "2")
        assert result.stdout.strip() == "ltrquest: not-a-command --threads 2"


class TestLauncherMounts:
    """The launcher exists to compute `--bind`, so that is what is asserted."""

    @pytest.fixture
    def fake_apptainer(self, tmp_path: Path) -> Path:
        """Stands in for apptainer and reports the mounts and command it got."""
        binned = tmp_path / "stubbin"
        binned.mkdir()
        _stub(
            binned / "apptainer",
            '#!/bin/sh\n'
            'echo "BIND=$APPTAINER_BIND"\n'
            'echo "ARGS=$*"\n',
        )
        return binned

    def launch(self, fake_apptainer: Path, cwd: Path, *args: str, image=None):
        env = dict(os.environ, PATH=f"{fake_apptainer}:{os.environ['PATH']}")
        env.pop("APPTAINER_BIND", None)
        env.pop("SINGULARITY_BIND", None)
        if image is not None:
            env["LTRQUEST_SIF"] = str(image)
        return subprocess.run(
            [str(LAUNCHER), *args],
            capture_output=True, text=True, cwd=str(cwd), env=env,
        )

    @pytest.fixture
    def image(self, tmp_path: Path) -> Path:
        path = tmp_path / "ltrquest.sif"
        path.write_text("not a real image")
        return path

    def test_the_real_working_directory_is_mounted(self, fake_apptainer, tmp_path, image):
        work = tmp_path / "real" / "work"
        work.mkdir(parents=True)
        link = tmp_path / "link"
        link.symlink_to(tmp_path / "real")

        # Enter through the symlink, exactly as a shell in /home/you/data does.
        result = self.launch(fake_apptainer, link / "work", "--genome", "g.fa", image=image)
        assert result.returncode == 0, result.stderr

        bind = next(l for l in result.stdout.splitlines() if l.startswith("BIND="))
        assert str(work.resolve()) in bind, (
            "the launcher must mount where the working directory really is, "
            f"not the symlinked path: {bind}"
        )

    def test_directories_of_existing_inputs_are_mounted(self, fake_apptainer, tmp_path, image):
        work = tmp_path / "work"
        work.mkdir()
        elsewhere = tmp_path / "scratch" / "genomes"
        elsewhere.mkdir(parents=True)
        genome = elsewhere / "hg38.fa"
        genome.write_text(">chr1\nACGT\n")

        result = self.launch(fake_apptainer, work, "--genome", str(genome), image=image)
        assert result.returncode == 0, result.stderr
        bind = next(l for l in result.stdout.splitlines() if l.startswith("BIND="))
        assert str(elsewhere.resolve()) in bind, bind

    def test_flag_values_that_are_not_paths_are_left_alone(self, fake_apptainer, tmp_path, image):
        work = tmp_path / "work"
        work.mkdir()
        result = self.launch(fake_apptainer, work, "--sdust-args", "-w 64 -t 15", image=image)
        assert result.returncode == 0, result.stderr
        bind = next(l for l in result.stdout.splitlines() if l.startswith("BIND="))
        assert bind == f"BIND={work.resolve()}", bind

    def test_ltrquest_is_the_default_command(self, fake_apptainer, tmp_path, image):
        work = tmp_path / "work"
        work.mkdir()
        result = self.launch(fake_apptainer, work, "--genome", "g.fa", image=image)
        args = next(l for l in result.stdout.splitlines() if l.startswith("ARGS="))
        assert args.endswith("ltrquest --genome g.fa"), args

    def test_a_named_stage_is_passed_through(self, fake_apptainer, tmp_path, image):
        work = tmp_path / "work"
        work.mkdir()
        result = self.launch(fake_apptainer, work, "ltrquest-gff3", "--help", image=image)
        args = next(l for l in result.stdout.splitlines() if l.startswith("ARGS="))
        assert args.endswith("ltrquest-gff3 --help"), args

    def test_a_missing_image_is_explained_not_crashed(self, fake_apptainer, tmp_path):
        work = tmp_path / "work"
        work.mkdir()
        result = self.launch(fake_apptainer, work, "--genome", "g.fa")
        assert result.returncode == 1
        assert "apptainer pull" in result.stderr
        assert "LTRQUEST_SIF" in result.stderr


class TestImageWiring:
    """The Dockerfile is what turns bin/entrypoint.sh into the SIF's runscript."""

    @pytest.fixture(scope="class")
    def dockerfile(self) -> str:
        return (REPO_ROOT / "Dockerfile").read_text()

    def test_an_entrypoint_is_declared(self, dockerfile):
        # An image with only a CMD is what broke `./ltrquest.sif --genome x.fa`.
        assert 'ENTRYPOINT ["/opt/ltrquest/bin/entrypoint"]' in dockerfile
        assert "ENTRYPOINT []" not in dockerfile

    def test_the_entrypoint_is_installed_where_it_is_declared(self, dockerfile):
        assert "/opt/ltrquest/src/bin/entrypoint.sh" in dockerfile
        assert "/opt/ltrquest/bin/entrypoint" in dockerfile

    def test_the_launcher_ships_inside_the_image(self, dockerfile):
        # The README tells people to `cat` it out of the image, so it has to be
        # there, at that path.
        assert "/opt/ltrquest/bin/ltrquest-container" in dockerfile

    def test_bare_flags_are_the_default_command(self, dockerfile):
        assert 'CMD ["--help"]' in dockerfile


class TestEnvironmentCarriesItsBuildDependencies:
    """`--run-trf` is on by default and TRF-mod has no conda package, so a
    conda install that cannot compile dies on the first real run."""

    @pytest.fixture(scope="class")
    def environment_yml(self) -> str:
        return (REPO_ROOT / "environment.yml").read_text()

    @pytest.mark.parametrize("package", ["git", "make", "c-compiler"])
    def test_build_tools_are_present(self, environment_yml, package):
        assert any(
            line.strip().startswith(f"- {package}")
            for line in environment_yml.splitlines()
        ), f"{package} missing from environment.yml"

    def test_sdust_is_installed_rather_than_compiled(self, environment_yml):
        assert any(
            line.strip().startswith("- sdust")
            for line in environment_yml.splitlines()
        )


class TestKmer2LTRRuntimeDependencies:
    """`parasail` on Bioconda is the C library. A recipe naming it resolves,
    installs, and still has no `import parasail` -- so Kmer2LTR dies the first
    time it aligns, several hours into a run."""

    RECIPES = ["environment.yml",
               "recipe/meta.yaml",
               "modules/local/ltrquest/detect/environment.yml",
               "modules/local/ltrquest/cluster/environment.yml"]

    @pytest.mark.parametrize("recipe", RECIPES)
    def test_the_python_binding_is_what_is_named(self, recipe):
        entries = [line.strip().lstrip("- ").split("#")[0].strip()
                   for line in (REPO_ROOT / recipe).read_text().splitlines()
                   if line.strip().startswith("- ")]
        assert "parasail-python" in entries, f"{recipe} does not name parasail-python"
        assert "parasail" not in entries, f"{recipe} names the C library instead"
        assert "pywfa" in entries, f"{recipe} does not name pywfa"

    @pytest.mark.parametrize("recipe", RECIPES)
    def test_neither_arrives_through_pip(self, recipe):
        text = (REPO_ROOT / recipe).read_text()
        assert "- pip:" not in text, f"{recipe} still installs through pip"


class TestMissingInputMessage:
    """`ERROR: Genome not found: hg38.fa` inside a container sends people
    looking for a missing file when the real problem is a missing mount."""

    def _run(self, driver: Path, tmp_path: Path, *args: str, in_container: bool):
        env = dict(os.environ)
        for key in ("APPTAINER_CONTAINER", "SINGULARITY_CONTAINER"):
            env.pop(key, None)
        if in_container:
            env["APPTAINER_CONTAINER"] = "/fake/ltrquest.sif"
        return subprocess.run(
            ["bash", str(driver), *args],
            capture_output=True, text=True, cwd=str(tmp_path), env=env,
        )

    def test_outside_a_container_the_message_stays_short(self, driver, tmp_path):
        result = self._run(driver, tmp_path, "--genome", "nope.fa", in_container=False)
        assert result.returncode == 1
        assert "Genome not found: nope.fa" in result.stderr
        assert "mounted" not in result.stderr

    def test_a_relative_path_gets_the_symlink_advice(self, driver, tmp_path):
        result = self._run(driver, tmp_path, "--genome", "nope.fa", in_container=True)
        assert result.returncode == 1
        assert "Genome not found: nope.fa" in result.stderr
        assert "cd -P ." in result.stderr
        assert "ltrquest-container" in result.stderr

    def test_an_absolute_path_names_the_directory_to_mount(self, driver, tmp_path):
        result = self._run(
            driver, tmp_path, "--genome", "/scratch/genomes/hg38.fa", in_container=True
        )
        assert result.returncode == 1
        assert "-B /scratch/genomes" in result.stderr

    def test_missing_proteins_are_explained_too(self, driver, tmp_path):
        genome = tmp_path / "g.fa"
        genome.write_text(">chr1\nACGT\n")
        result = self._run(
            driver, tmp_path,
            "--genome", str(genome), "--proteins", "/scratch/prot.fa",
            in_container=True,
        )
        assert result.returncode == 1
        assert "Proteins not found" in result.stderr
        assert "-B /scratch" in result.stderr
