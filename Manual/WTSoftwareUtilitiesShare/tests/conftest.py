from __future__ import annotations

import sys
from pathlib import Path


def pytest_configure() -> None:
    """
    Keep imports simple for students: tests run against the in-repo `windtunnel/`
    package without requiring a packaging step.
    """
    project_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(project_root))
    # Allow imports from `tests/support` without making `tests/` a public package.
    sys.path.insert(0, str(project_root / "tests"))


