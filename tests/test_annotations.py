import pathlib
from collections import namedtuple
from typing import NamedTuple, List, Dict

import numpy as np
import pandas as pd
import pytest
from _pytest.fixtures import fixture
from scipy.stats import fisher_exact

from yoke import create
from yoke.create import MGS
from yoke.annotations import (
    HierarchicalSchema, RTKModules, cli_function_modules, function_enrichment,
    cli_function_enrichment
)

AnnPaths = namedtuple("AnnPaths", "annotation description")

@fixture
def kegg_annotation_path() -> NamedTuple:
    """Path to a KEGG and description for initial tests."""

    annotation_path: pathlib.Path = pathlib.Path(
        "/Volumes/DIME6/Anno/Func/DIAass_KGM.srt.gzgeneAss.gz",
    )
    hier_path: pathlib.Path = pathlib.Path(
        "/Volumes/hpc-home/downloads/ko_KO.descr"
    )

    return AnnPaths(
        annotation=annotation_path,
        description=hier_path
    )

@fixture
def magvsgc_path() -> pathlib.Path:
    """Path to a MAGvsGC file for testing"""

    magvsgc: pathlib.Path = pathlib.Path(
        "/Volumes/DIME6/Bin3SB/LOGandSUB/MAGvsGC.txt.gz"
    )
    return magvsgc

@fixture
def dime_magvsgc_path() -> pathlib.Path:
    """Path to the MAGvsGC for DIME, a mid-size study."""

    magvsgc: pathlib.Path = pathlib.Path(
        "data/DIME/MAGvsGC.txt.gz"
    )
    return magvsgc

@fixture
def dime_nog_annotation_path() -> NamedTuple:
    """Path for eggNOG annotation and description for DIME."""

    annotation_path: pathlib.Path = pathlib.Path(
        "data/DIME/DIAass_NOG.srt.gzgeneAss.gz"
    )
    hier_path: pathlib.Path = pathlib.Path(
        "data/DIME/NOG.descr"
    )

    return AnnPaths(
        annotation=annotation_path,
        description=hier_path
    )

@fixture
def dime_kegg_annotation_path() -> NamedTuple:
    """Path for eggNOG annotation and description for DIME."""

    annotation_path: pathlib.Path = pathlib.Path(
        "data/DIME/DIAass_KGM.srt.gzgeneAss.gz"
    )
    hier_path: pathlib.Path = pathlib.Path(
        "data/DIME/KEGG.descr"
    )

    return AnnPaths(
        annotation=annotation_path,
        description=hier_path
    )

def test_hs_kegg_init(dime_kegg_annotation_path):
    """Initialise a KEGG functional annotation from DataFrames."""

    annotation: pd.DataFrame = pd.read_csv(
        dime_kegg_annotation_path.annotation,
        sep="\t",
        header=None,
        index_col=0
    )
    description: pd.DataFrame = pd.read_csv(
        dime_kegg_annotation_path.description,
        sep="\t",
        header=None,
        index_col=0
    )

    s: HierarchicalSchema = HierarchicalSchema(
        name="KEGG",
        hierarchy=description,
        gene_map=annotation
    )

def test_hs_mgs_count(
        dime_kegg_annotation_path,
        dime_magvsgc_path
):
    mgs: Dict[str, MGS] = create.read_mag_vs_gc(
        dime_magvsgc_path,
        min_completeness=0,
        max_contamination=100
    )

    annotation: pd.DataFrame = pd.read_csv(
        dime_kegg_annotation_path.annotation,
        sep="\t",
        header=None,
        index_col=0
    )
    description: pd.DataFrame = pd.read_csv(
        dime_kegg_annotation_path.description,
        sep="\t",
        header=None,
        index_col=0
    )

    s: HierarchicalSchema = HierarchicalSchema(
        name="KEGG",
        hierarchy=description,
        gene_map=annotation
    )
    mgs_list: List[MGS] = list(mgs.values())[0:10]
    mgs: MGS = mgs_list[0]

    # f = s.function_counts(mgs)
    tbl: pd.DataFrame = s.function_counts(mgs_list)

@fixture
def pth_kegg_abd() -> pathlib.Path:
    return pathlib.Path("data/DIME/KGML0.txt")

def test_score_gbm(
        pth_kegg_abd: pathlib.Path,
        tmp_path: pathlib.Path
):
    gbm = RTKModules.from_bundled("GBM")
    pref: pathlib.Path = tmp_path / "test_gbm"
    gbm.score_table(
        pth_mgs_tbl=pth_kegg_abd,
        output_prefix=pref
    )

    for xtn in ["KOused", "mat", "descr", "MODscore"]:
        assert pref.with_suffix(f".{xtn}").exists()

@fixture
def pth_mgs_kegg() -> pathlib.Path:
    return pathlib.Path("output/KEGG")

def test_score_cli(
        pth_mgs_kegg: pathlib.Path,
        tmp_path: pathlib.Path
):
    """Attempt to score multiple tables using the CLI function."""

    kegg_mgs: List[pathlib.Path] = list(pth_mgs_kegg.glob("MGS.1*.tsv.gz"))
    cli_function_modules(
        kegg_mgs, "GBM", tmp_path, "GBM"
    )
    cli_function_modules(
        kegg_mgs, "BSB", tmp_path, "BSB"
    )
    for tbl in kegg_mgs:
        assert (
                tmp_path /
                f"{tbl.name.replace(".tsv.gz", "")}.GBM.mat"
        ).exists()
        assert (
                tmp_path /
                f"{tbl.name.replace(".tsv.gz", "")}.BSB.mat"
        ).exists()
        # Ensure uncompressed versions tidied up
        assert not tbl.with_suffix("").exists()

@fixture
def enrich_counts() -> pd.DataFrame:
    """Function counts for 8 MAGs with metadata, plus one MAG (M8) without."""

    return pd.DataFrame(
        {
            "K_a": [1, 2, 1, 3, 0, 0, 0, 0, 5],
            "K_all": [1, 1, 1, 1, 1, 1, 1, 1, 5],
            "K_mix": [1, 0, 1, 0, 1, 1, 0, 0, 5],
        },
        index=[f"M{i}" for i in range(9)]
    ).T

@fixture
def enrich_metadata() -> pd.DataFrame:
    """Metadata with a two-level, a boolean with NA, and a three-level
    column."""

    return pd.DataFrame(
        {
            "grp": ["a", "a", "a", "a", "b", "b", "b", "b"],
            "flag": [True, True, True, False, False, False, None, False],
            "three": ["x", "y", "z", "x", "y", "z", "x", "y"],
        },
        index=[f"M{i}" for i in range(8)]
    )

def test_enrichment_matches_scipy(enrich_counts, enrich_metadata):
    """Statistics match scipy on the equivalent 2x2 table."""

    res: pd.DataFrame = function_enrichment(
        enrich_counts, enrich_metadata, tests=["fisher"]
    )
    # K_a present in all 4 of group a, none of group b
    expected = fisher_exact(np.array([[4, 0], [0, 4]]))
    assert res.loc[("grp", "K_a"), "fisher_stat"] == expected.statistic
    assert res.loc[("grp", "K_a"), "fisher_p"] == pytest.approx(
        expected.pvalue
    )

def test_enrichment_columns_and_index(enrich_counts, enrich_metadata):
    """Output has one row per tested column and function, three level
    columns are skipped, and each test gives stat, p and adjusted p."""

    tests: List[str] = ["fisher", "barnard", "boschloo"]
    res: pd.DataFrame = function_enrichment(
        enrich_counts, enrich_metadata, tests=tests
    )
    assert list(res.index.names) == ["metadata", "function"]
    assert set(res.index.get_level_values("metadata")) == {"grp", "flag"}
    assert set(res.index.get_level_values("function")) == set(
        enrich_counts.index
    )
    assert list(res.columns) == [
        f"{t}_{s}" for t in tests for s in ["stat", "p", "p_adj"]
    ]
    # Undefined p-values (e.g. boschloo when a function is in every MAG) are
    # treated as 1 for correction
    assert (res.filter(like="_p_adj").values >= res.filter(
        regex="_p$").fillna(1.0).values).all()

def test_enrichment_drops_na(enrich_counts, enrich_metadata):
    """NA metadata values are excluded, and True is the first level."""

    res: pd.DataFrame = function_enrichment(
        enrich_counts, enrich_metadata, tests=["fisher"]
    )
    # M6 is NA so excluded. True: M0-M2, K_mix in 2. False: M3, M4, M5, M7,
    # K_mix in 2.
    expected = fisher_exact(np.array([[2, 1], [2, 2]]))
    assert res.loc[("flag", "K_mix"), "fisher_stat"] == pytest.approx(
        expected.statistic
    )
    assert res.loc[("flag", "K_mix"), "fisher_p"] == pytest.approx(
        expected.pvalue
    )

def test_enrichment_unmatched_samples(enrich_counts, enrich_metadata):
    """MAGs not in metadata do not affect results."""

    with_extra: pd.DataFrame = function_enrichment(
        enrich_counts, enrich_metadata
    )
    without: pd.DataFrame = function_enrichment(
        enrich_counts.drop(columns="M8"), enrich_metadata
    )
    pd.testing.assert_frame_equal(with_extra, without)

def test_enrichment_unknown_test(enrich_counts, enrich_metadata):
    with pytest.raises(ValueError):
        function_enrichment(enrich_counts, enrich_metadata, tests=["chi2"])

def test_enrichment_no_valid_columns(enrich_counts, enrich_metadata):
    """No two-level columns gives an empty result."""

    res: pd.DataFrame = function_enrichment(
        enrich_counts, enrich_metadata[["three"]]
    )
    assert res.empty

def test_enrichment_cli(
        enrich_counts: pd.DataFrame,
        enrich_metadata: pd.DataFrame,
        tmp_path: pathlib.Path
):
    """Count tables are combined, metadata columns subset, and results
    written."""

    pth_a: pathlib.Path = tmp_path / "a.tsv.gz"
    pth_b: pathlib.Path = tmp_path / "b.tsv"
    pth_md: pathlib.Path = tmp_path / "metadata.tsv"
    pth_out: pathlib.Path = tmp_path / "out" / "enrichment.tsv"
    enrich_counts.iloc[:, :4].to_csv(pth_a, sep="\t", compression="gzip")
    # Drop a function from the second table, should be treated as absent
    enrich_counts.iloc[:, 4:].drop(index="K_a").to_csv(pth_b, sep="\t")
    enrich_metadata.to_csv(pth_md, sep="\t")

    cli_function_enrichment(
        counts=[pth_a, pth_b],
        metadata=pth_md,
        output=pth_out,
        tests=["fisher"],
        columns=["grp"]
    )

    res: pd.DataFrame = pd.read_csv(pth_out, sep="\t", index_col=[0, 1])
    expected: pd.DataFrame = function_enrichment(
        enrich_counts, enrich_metadata[["grp"]], tests=["fisher"]
    )
    pd.testing.assert_frame_equal(
        res.sort_index(), expected.sort_index(), check_names=False
    )

    with pytest.raises(ValueError):
        cli_function_enrichment(
            counts=[pth_a], metadata=pth_md, output=pth_out,
            columns=["missing"]
        )

def test_enrichment_bad_adj_method(enrich_counts, enrich_metadata):
    with pytest.raises(ValueError, match="correction method"):
        function_enrichment(enrich_counts, enrich_metadata, adj_method="no")

def test_enrichment_no_overlap(enrich_counts, enrich_metadata):
    """No shared identifiers is an error rather than an empty result."""

    md: pd.DataFrame = enrich_metadata.set_axis(
        [f"X{i}" for i in range(8)]
    )
    with pytest.raises(ValueError, match="No MGS/MAGs"):
        function_enrichment(enrich_counts, md)

@pytest.mark.parametrize("where", ["metadata", "columns", "functions"])
def test_enrichment_duplicates(enrich_counts, enrich_metadata, where):
    if where == "metadata":
        enrich_metadata = pd.concat(
            [enrich_metadata, enrich_metadata.iloc[:1]]
        )
    elif where == "columns":
        enrich_counts = pd.concat(
            [enrich_counts, enrich_counts[["M0"]]], axis=1
        )
    else:
        enrich_counts = pd.concat([enrich_counts, enrich_counts.iloc[:1]])
    with pytest.raises(ValueError, match="Duplicate"):
        function_enrichment(enrich_counts, enrich_metadata)

def test_enrichment_non_numeric(enrich_counts, enrich_metadata):
    with pytest.raises(ValueError, match="non-numeric"):
        function_enrichment(enrich_counts.astype(str), enrich_metadata)

def test_enrichment_no_functions(enrich_counts, enrich_metadata):
    res: pd.DataFrame = function_enrichment(
        enrich_counts.iloc[:0], enrich_metadata
    )
    assert res.empty

def test_enrichment_mixed_type_levels(enrich_counts, enrich_metadata):
    """Levels of different types are still tested."""

    enrich_metadata["grp"] = [1, 1, 1, 1, "b", "b", "b", "b"]
    res: pd.DataFrame = function_enrichment(
        enrich_counts, enrich_metadata[["grp"]], tests=["fisher"]
    )
    assert res.loc[("grp", "K_a"), "fisher_p"] == pytest.approx(
        fisher_exact(np.array([[4, 0], [0, 4]])).pvalue
    )

def test_enrichment_level_lost_after_matching(
        enrich_counts, enrich_metadata, caplog
):
    """A column is skipped, with the levels named, when MGS/MAGs for one
    level are not in the counts."""

    res: pd.DataFrame = function_enrichment(
        enrich_counts[["M0", "M1", "M2", "M3"]], enrich_metadata[["grp"]]
    )
    assert res.empty
    assert "1 levels (a)" in caplog.text

def test_enrichment_nan_counts_absent(enrich_counts, enrich_metadata):
    """NaN counts are treated as the function being absent."""

    with_nan: pd.DataFrame = function_enrichment(
        enrich_counts.where(enrich_counts > 0), enrich_metadata
    )
    pd.testing.assert_frame_equal(
        with_nan, function_enrichment(enrich_counts, enrich_metadata)
    )

def test_enrichment_absent_everywhere(enrich_counts, enrich_metadata):
    """A function in no MAGs gives p of 1."""

    enrich_counts.loc["K_none"] = 0
    res: pd.DataFrame = function_enrichment(
        enrich_counts, enrich_metadata, tests=["fisher"]
    )
    assert res.loc[("grp", "K_none"), "fisher_p"] == 1.0
    assert res.loc[("grp", "K_none"), "fisher_p_adj"] == 1.0

def test_enrichment_cli_numeric_ids(
        enrich_counts: pd.DataFrame,
        enrich_metadata: pd.DataFrame,
        tmp_path: pathlib.Path
):
    """Numeric-looking identifiers match between counts and metadata,
    including leading zeros."""

    ids: List[str] = [f"00{i}" for i in range(9)]
    pth_c: pathlib.Path = tmp_path / "counts.tsv"
    pth_md: pathlib.Path = tmp_path / "metadata.tsv"
    pth_out: pathlib.Path = tmp_path / "out.tsv"
    enrich_counts.set_axis(ids, axis=1).to_csv(pth_c, sep="\t")
    enrich_metadata.set_axis(ids[:8]).to_csv(pth_md, sep="\t")

    cli_function_enrichment(
        counts=[pth_c], metadata=pth_md, output=pth_out, tests=["fisher"]
    )
    res: pd.DataFrame = pd.read_csv(pth_out, sep="\t", index_col=[0, 1])
    assert set(res.index.get_level_values(0)) == {"grp", "flag"}

def test_enrichment_cli_duplicates(
        enrich_counts: pd.DataFrame,
        enrich_metadata: pd.DataFrame,
        tmp_path: pathlib.Path
):
    """Duplicate functions within a table, or MAGs across tables, are
    reported."""

    pth_c: pathlib.Path = tmp_path / "counts.tsv"
    pth_dup: pathlib.Path = tmp_path / "dup.tsv"
    pth_md: pathlib.Path = tmp_path / "metadata.tsv"
    pth_out: pathlib.Path = tmp_path / "out.tsv"
    enrich_counts.to_csv(pth_c, sep="\t")
    pd.concat([enrich_counts, enrich_counts.iloc[:1]]).to_csv(
        pth_dup, sep="\t"
    )
    enrich_metadata.to_csv(pth_md, sep="\t")

    with pytest.raises(ValueError, match="function IDs in .*dup.tsv"):
        cli_function_enrichment(
            counts=[pth_dup], metadata=pth_md, output=pth_out
        )
    with pytest.raises(ValueError, match="across function count tables"):
        cli_function_enrichment(
            counts=[pth_c, pth_c], metadata=pth_md, output=pth_out
        )
