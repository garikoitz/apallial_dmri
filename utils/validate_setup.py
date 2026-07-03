#!/usr/bin/env python
"""Validate that MBRAIN repository setup files are present and coherent."""

from __future__ import annotations

from pathlib import Path
import sys


REPO_DIR = Path(__file__).resolve().parent.parent


def check(condition: bool, name: str) -> bool:
    if condition:
        print(f"✅ {name}")
        return True
    print(f"❌ {name}")
    return False


def main() -> int:
    print("MBRAIN Setup Validation")
    passed = 0
    failed = 0

    checks = [
        (REPO_DIR.exists(), "Repository root exists"),
        ((REPO_DIR / "notebooks").exists(), "notebooks/ exists"),
        ((REPO_DIR / "functions").exists(), "functions/ exists"),
        ((REPO_DIR / "utils").exists(), "utils/ exists"),
        ((REPO_DIR / "config1.yaml").exists(), "config1.yaml exists"),
        ((REPO_DIR / "utils" / "setup_data.py").exists(), "utils/setup_data.py exists"),
        ((REPO_DIR / "utils" / "create_data_package.py").exists(), "utils/create_data_package.py exists"),
        ((REPO_DIR / "utils" / "validate_setup.py").exists(), "utils/validate_setup.py exists"),
        ((REPO_DIR / "notebooks" / "pipeline_groups_new_mask.py").exists(), "pipeline notebook exists"),
        ((REPO_DIR / "SETUP.md").exists(), "SETUP.md exists"),
        ((REPO_DIR / "QUICKSTART.md").exists(), "QUICKSTART.md exists"),
    ]

    for ok, name in checks:
        if check(ok, name):
            passed += 1
        else:
            failed += 1

    nb = REPO_DIR / "notebooks" / "pipeline_groups_new_mask.py"
    if nb.exists():
        text = nb.read_text(encoding="utf-8", errors="ignore")
        for token in ["PRE-STEP -2", "PRE-STEP -1a", "PRE-STEP -1b", "PRE-STEP -1c", "platform.system()"]:
            if check(token in text, f"Notebook contains '{token}'"):
                passed += 1
            else:
                failed += 1

    print(f"\nPassed: {passed}")
    print(f"Failed: {failed}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
