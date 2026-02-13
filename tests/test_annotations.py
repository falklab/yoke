import pathlib
from collections import namedtuple
from typing import NamedTuple, List, Dict

import pandas as pd
from _pytest.fixtures import fixture

import create
from create import MGS
from annotations import HierarchicalSchema, RTKModules, cli_function_modules

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
    for tbl in kegg_mgs:
        assert (
                tmp_path /
                f"{tbl.name.replace(".tsv.gz", "")}.GBM.mat"
        ).exists()
        # Ensure uncompressed versions tidied up
        assert not tbl.with_suffix("").exists()