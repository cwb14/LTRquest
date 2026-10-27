"""A round's per-chunk work folders are removed in the background."""

from ltrquest import detect


def _tree(root, n_dirs=20, n_files=5):
    for d in range(n_dirs):
        sub = root / f"chunk{d:03d}"
        sub.mkdir(parents=True)
        for f in range(n_files):
            (sub / f"f{f}").write_text("x")
    (root / "top.txt").write_text("x")


def test_a_folder_leaves_its_name_at_once_and_is_gone_after_the_wait(tmp_path):
    work = tmp_path / "ltr_runs"
    _tree(work)
    (tmp_path / "keep.txt").write_text("x")
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "precious").write_text("x")
    (work / "chunk000" / "link").symlink_to(outside, target_is_directory=True)
    (work / "dirlink").symlink_to(outside, target_is_directory=True)

    detect._clean_dir_in_background(work)
    assert not work.exists()
    detect._wait_for_background_removals()
    assert sorted(p.name for p in tmp_path.iterdir()) == ["keep.txt", "outside"]
    assert (outside / "precious").exists()      # links are removed, never followed


def test_a_file_or_a_missing_path_is_removed_in_place(tmp_path):
    fa = tmp_path / "x.fa"
    fa.write_text(">a\nACGT\n")
    detect._clean_dir_in_background(fa)
    detect._clean_dir_in_background(tmp_path / "missing")
    assert not fa.exists()
    detect._wait_for_background_removals()
