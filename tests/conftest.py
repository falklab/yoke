import os
import shutil
from pathlib import Path

import pytest

LOCAL_DATA_ENV: str = "YOKE_LOCAL_DATA"
"""Environment variable giving the directory of large local test data."""


def pytest_collection_modifyitems(
    config: pytest.Config, items: list[pytest.Item]
) -> None:
    """Skip tests whose external requirements are not available."""
    local_data = os.environ.get(LOCAL_DATA_ENV)
    has_local_data = local_data is not None and Path(local_data).is_dir()
    has_rtk = shutil.which("rtk") is not None

    skip_local_data = pytest.mark.skip(
        reason=f"set {LOCAL_DATA_ENV} to a directory of local test data"
    )
    skip_rtk = pytest.mark.skip(reason="rtk not found on PATH")
    for item in items:
        if "local_data" in item.keywords and not has_local_data:
            item.add_marker(skip_local_data)
        if "requires_rtk" in item.keywords and not has_rtk:
            item.add_marker(skip_rtk)


@pytest.fixture(autouse=True)
def _run_from_repo_root(
    request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Run every test from the repository root.

    Some test data (e.g. assembly.txt) contains paths relative to the root.
    """
    monkeypatch.chdir(request.config.rootpath)


@pytest.fixture
def data_dir() -> Path:
    """Directory of small test data files committed to the repository."""
    return Path(__file__).parent / "data"


@pytest.fixture
def expected_dir() -> Path:
    """Directory of expected (golden) test outputs."""
    return Path(__file__).parent / "expected"


@pytest.fixture
def local_data_dir() -> Path:
    """Directory of large local test data, not committed to the repository.

    Tests using this should be marked `local_data`, so they are skipped when
    the directory is not available.
    """
    return Path(os.environ[LOCAL_DATA_ENV])
