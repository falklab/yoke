import gzip
from pathlib import Path

import pytest

from yoke import MAG
from yoke.extract import MAGSequences


@pytest.fixture
def example_mag() -> MAG:
    """A single MAG for tests."""
    return MAG(
        name="Sample1__SemiBin_177.fa.gz",
        representative=False,
        completeness=100,
        contamination=0,
        contigs=[],
    )


@pytest.fixture
def example_mags() -> list[MAG]:
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


@pytest.fixture
def example_assembly_dir(data_dir: Path) -> Path:
    """Assembly directory for tests."""
    return data_dir.resolve()


def test_from_assembly_dir(
    example_mag: MAG, example_assembly_dir: Path, data_dir: Path
):
    mag_seq = MAGSequences.from_assembly_dir(example_mag, example_assembly_dir)
    assert mag_seq.mag is example_mag
    # Sample1 is part of a coassembly, so is resolved via assembly.txt
    assert mag_seq.sample_assembly == data_dir / "Assmbl_Exmpl" / "metag"
    # Contig ID list should start empty
    assert isinstance(mag_seq.contig_ids, list)
    assert len(mag_seq.contig_ids) == 0


def test_group_mags(
    example_mags: list[MAG], example_assembly_dir: Path, data_dir: Path
):
    mag_seqs = [
        MAGSequences.from_assembly_dir(x, example_assembly_dir)
        for x in example_mags
    ]
    grouped = MAGSequences.group_mags(mag_seqs)

    # Test expected groups identified
    expected_mags = dict(
        zip(
            [
                data_dir / "Assmbl_Exmpl" / "metag",
                data_dir / "Sample2" / "assemblies" / "metag",
            ],
            [
                {"Sample1__SemiBin_177.fa.gz", "Sample1__SemiBin_1.fa.gz"},
                {"Sample2__SemiBin_666.fa.gz"},
            ],
            strict=True,
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
    example_mags: list[MAG], example_assembly_dir: Path, tmp_path: Path
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
    example_mags: list[MAG],
    example_assembly_dir: Path,
    expected_dir: Path,
    tmp_path: Path,
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
    expected_path = expected_dir / "test_extract_contigs_sample"
    assert expected_path.exists()
    for exp_file in expected_path.glob("*.fna"):
        exp_out = exp_file.name + ".gz"
        matched = [x for x in fastas if x.name == exp_out]
        assert len(matched) == 1
        # Could instead hash... but I only thought of it after
        with gzip.open(matched[0], "rt") as r, open(exp_file) as e:
            r_lines, e_lines = r.readlines(), e.readlines()
            assert len(r_lines) == len(e_lines)
            for a, b in zip(r_lines, e_lines, strict=True):
                assert a == b
