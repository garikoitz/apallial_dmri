from pathlib import Path

from functions.pipeline_new_mask.exec_utils import check_file, run_parallel


def test_check_file_reports_missing(tmp_path: Path, capsys):
    missing = tmp_path / "missing.txt"
    ok = check_file(missing, "missing")
    out = capsys.readouterr().out
    assert ok is False
    assert "MISSING" in out


def test_check_file_reports_existing(tmp_path: Path, capsys):
    existing = tmp_path / "file.txt"
    existing.write_text("abc", encoding="utf-8")
    ok = check_file(existing, "existing")
    out = capsys.readouterr().out
    assert ok is True
    assert "existing" in out


def _fake_worker_success(args):
    sid, _cmd, _timeout, _log_dir, _label = args
    return sid, 0, "", ""


def _fake_worker_mixed(args):
    sid, _cmd, _timeout, _log_dir, _label = args
    if sid == "sub-b":
        return sid, 1, "", "boom"
    return sid, 0, "", ""


def test_run_parallel_all_success(tmp_path: Path):
    tasks = [("sub-a", "echo a", 1), ("sub-b", "echo b", 1)]
    result = run_parallel(
        tasks,
        label="Unit",
        max_workers=2,
        log_dir=tmp_path,
        worker_fn=_fake_worker_success,
    )
    assert sorted(result["success"]) == ["sub-a", "sub-b"]
    assert result["failed"] == []


def test_run_parallel_mixed_results(tmp_path: Path):
    tasks = [("sub-a", "echo a", 1), ("sub-b", "echo b", 1)]
    result = run_parallel(
        tasks,
        label="Unit",
        max_workers=2,
        log_dir=tmp_path,
        worker_fn=_fake_worker_mixed,
    )
    assert result["success"] == ["sub-a"]
    assert result["failed"] == ["sub-b"]
