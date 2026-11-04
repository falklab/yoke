"""Run the fast tests related to files changed in a commit.

Used by the ``pytest-changed`` pre-commit hook, which passes the staged
Python files as arguments. Each file is mapped to tests by naming
convention:

- ``src/yoke/<module>.py`` runs ``tests/test_<module>.py``
- ``src/yoke/data/**`` runs ``tests/test_data.py``
- ``tests/test_*.py`` runs itself
- ``tests/conftest.py`` runs all tests

Other files run no tests. Tests marked ``slow``, ``local_data`` or
``requires_rtk`` are deselected. CI runs the full suite, so this is only a
quick local check.
"""

import pathlib
import subprocess
import sys

TESTS: pathlib.Path = pathlib.Path("tests")
SRC: pathlib.Path = pathlib.Path("src") / "yoke"
FAST_MARKERS: str = "not slow and not local_data and not requires_rtk"
NO_TESTS_COLLECTED: int = 5


def _tests_for(pth: pathlib.Path) -> list[pathlib.Path]:
    """Return the test files which cover a changed file."""
    if pth == TESTS / "conftest.py":
        return sorted(TESTS.glob("test_*.py"))
    if pth.parent == TESTS and pth.name.startswith("test_"):
        return [pth]
    if pth.is_relative_to(SRC / "data"):
        return [TESTS / "test_data.py"]
    if pth.parent == SRC:
        return [TESTS / f"test_{pth.stem}.py"]
    return []


def main(changed: list[str]) -> int:
    """Run pytest on the tests for changed files, returning the exit code."""
    tests: list[pathlib.Path] = sorted(
        {t for f in changed for t in _tests_for(pathlib.Path(f))}
    )
    tests = [t for t in tests if t.exists()]
    if not tests:
        return 0
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-q",
            "--no-cov",
            "-m",
            FAST_MARKERS,
            *(str(t) for t in tests),
        ]
    )
    if result.returncode == NO_TESTS_COLLECTED:
        return 0
    return result.returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
