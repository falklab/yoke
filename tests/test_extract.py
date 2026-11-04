import gzip
from pathlib import Path
from typing import List

import pyfastx
from _pytest.fixtures import fixture

from yoke import MAG
from yoke.extract import MAGSequences, cli_extract_contigs


@fixture
def example_mag() -> MAG:
    """A single MAG for tests."""
    return MAG(
        name="Sample1__SemiBin_177.fa.gz",
        representative=False,
        completeness=100,
        contamination=0,
        contigs=[],
    )


@fixture
def example_mags() -> List[MAG]:
    """A list of MAGs originating from two different assemblies for tests."""

    return [
        MAG(
            name="Sample1__SemiBin_177.fa.gz",
            representative=False,
            completeness=100,
            contamination=0,
            contigs=[],
        ),
        MAG(
            name="Sample1__SemiBin_1.fa.gz",
            representative=False,
            completeness=100,
            contamination=0,
            contigs=[],
        ),
        MAG(
            name="Sample2__SemiBin_666.fa.gz",
            representative=True,
            completeness=100,
            contamination=100,
            contigs=[],
        ),
    ]


@fixture
def example_assembly_dir() -> Path:
    """Assembly directory for tests."""
    return (Path("tests") / "data").resolve()


def test_from_assembly_dir(example_mag: MAG, example_assembly_dir):
    mag_seq = MAGSequences.from_assembly_dir(example_mag, example_assembly_dir)
    assert mag_seq.mag is example_mag
    resolved_path = Path("tests", "data", "Assmbl_Exmpl", "metag").absolute()
    # Contig ID list should start empty
    assert isinstance(mag_seq.contig_ids, list)
    assert len(mag_seq.contig_ids) == 0


def test_group_mags(example_mags: List[MAG], example_assembly_dir: Path):
    mag_seqs = [
        MAGSequences.from_assembly_dir(x, example_assembly_dir)
        for x in example_mags
    ]
    grouped = MAGSequences.group_mags(mag_seqs)

    # Test expected groups identified
    expected_mags = dict(
        zip(
            [
                Path("tests", "data", "Assmbl_Exmpl", "metag").absolute(),
                Path(
                    "tests", "data", "Sample2", "assemblies", "metag"
                ).absolute(),
            ],
            [
                {"Sample1__SemiBin_177.fa.gz", "Sample1__SemiBin_1.fa.gz"},
                {"Sample2__SemiBin_666.fa.gz"},
            ],
        )
    )
    # Test groups have expected contents
    for grp, mags in expected_mags.items():
        assert any(x == grp for x in grouped.keys())
        assert set(x.mag.name for x in grouped[grp]) == mags


def test_add_contig_ids(example_mags, example_assembly_dir: Path):
    same_origin = [
        MAGSequences.from_assembly_dir(x, example_assembly_dir)
        for x in example_mags[:2]
    ]
    mags = MAGSequences.add_contig_ids(same_origin)
    # Test have expected (simple) list of contigs
    assert mags[0].contig_ids == ["Is_in_MAG"]
    assert mags[1].contig_ids == ["Is_in_MAG2"]


def test_extract_contigs(
    example_mags: List[MAG], example_assembly_dir: Path, tmp_path: Path
):
    mag_seqs = [
        MAGSequences.from_assembly_dir(x, example_assembly_dir)
        for x in example_mags
    ]
    created = MAGSequences.extract_contigs(
        mag_seqs, tmp_path, processes=2, progress_bar=True
    )
    assert len(created) == len(example_mags)


def test__extract_contigs_sample(
    example_mags: List[MAG], example_assembly_dir: Path, tmp_path: Path
):
    """Test that we can extract contigs from a single assembly."""
    tmp_out = tmp_path / "fasta"
    tmp_out.mkdir(parents=True, exist_ok=True)
    mags = [
        MAGSequences.from_assembly_dir(x, example_assembly_dir)
        for x in example_mags[:2]
    ]
    mags = MAGSequences.add_contig_ids(mags)
    fastas = MAGSequences._extract_contigs_sample(mags, tmp_out)

    # Check output is as expected - these are very small FASTA files
    # for which expected version have been created
    expected_path = Path(
        "tests", "expected", "test_extract_contigs_sample"
    ).absolute()
    assert expected_path.exists()
    for exp_file in expected_path.glob("*.fna"):
        exp_out = exp_file.name + ".gz"
        matched = [x for x in fastas if x.name == exp_out]
        assert len(matched) == 1
        # Could instead hash... but I only thought of it after
        with gzip.open(matched[0], "rt") as r, open(exp_file, "rt") as e:
            r_lines, e_lines = r.readlines(), e.readlines()
            assert len(r_lines) == len(e_lines)
            for a, b in zip(r_lines, e_lines):
                assert a == b
