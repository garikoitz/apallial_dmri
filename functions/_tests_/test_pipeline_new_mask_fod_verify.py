from pathlib import Path

from functions.pipeline_new_mask.fod_verify import collect_fod_verification, collect_group_response_status


class _Result:
    def __init__(self, stdout: str):
        self.stdout = stdout


def _fake_run(cmd: str, timeout=30, verbose=False):
    if "mask_anat" in cmd:
        return _Result("28 100 107\n")
    return _Result("28 100 107 45\n")


def _subject(base: Path, sid: str = "sub-01") -> dict:
    dwi = base / sid / "ses-1" / "dwi"
    dwi.mkdir(parents=True, exist_ok=True)
    return {"id": sid, "group": "CTRL", "dwi_dir": dwi}


def test_collect_fod_verification_missing_norm_group(tmp_path: Path):
    s = _subject(tmp_path)
    rows, issues = collect_fod_verification([s], run_cmd=_fake_run, modules="module load mrtrix")
    assert rows[0]["status"] == "MISSING"
    assert any("wm_fod_norm_group.mif MISSING" in x for x in issues)


def test_collect_fod_verification_ok_row(tmp_path: Path):
    s = _subject(tmp_path)
    dwi = s["dwi_dir"]
    for f in ["wm_fod_norm_group.mif", "wm_fod_group.mif", "mask_anat.mif", "response_wm.txt", "response_gm.txt", "response_csf.txt"]:
        (dwi / f).touch()

    rows, issues = collect_fod_verification([s], run_cmd=_fake_run, modules="module load mrtrix")
    assert rows[0]["status"] == "OK"
    assert rows[0]["grid_match"] == "✅"
    assert rows[0]["responses"] == "✅"
    assert not any("GRID MISMATCH" in x for x in issues)


def test_collect_group_response_status(tmp_path: Path):
    grp = tmp_path / "group_responses"
    grp.mkdir(parents=True)
    (grp / "group_response_wm.txt").touch()

    status = collect_group_response_status(tmp_path)
    wm = next(x for x in status if x["name"] == "group_response_wm.txt")
    gm = next(x for x in status if x["name"] == "group_response_gm.txt")
    assert wm["exists"] is True
    assert gm["exists"] is False
