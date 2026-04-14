from pathlib import Path

from _pytest.fixtures import fixture

from yoke import MAG
from yoke.extract import MAGSequences


@fixture
def example_mag() -> MAG:
    """A single MAG for tests."""
    return MAG(
        name = "Sample1__SemiBin_177.fa.gz",
        representative = False,
        completeness = 100,
        contamination = 0,
        contigs = []
    )

@fixture
def example_assembly_dir() -> Path:
    """Assembly directory for tests."""
    return (Path("tests") / "data").resolve()

def test_get_contig_link_file():
    assert False


def test_contig_fasta():
    assert False


def test_from_assembly_dir(
        example_mag: MAG,
        example_assembly_dir
):
    mag_seq = MAGSequences.from_assembly_dir(
        example_mag,
        example_assembly_dir
    )


def test_group_mags():
    assert False


def test_add_contig_ids():
    assert False


def test_extract_contigs():
    assert False


def test__extract_contigs_sample():
    assert False
