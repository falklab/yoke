import pathlib
from collections import namedtuple
from typing import NamedTuple, List, Dict

import pandas as pd
from _pytest.fixtures import fixture

import create
from annotations import HierarchicalSchema


@fixture
def kegg_annotation_path() -> NamedTuple:
    """Path to a KEGG and description for initial tests."""
    AnnPaths = namedtuple("AnnPaths", "annotation description")
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
def magvsgc_path() -> NamedTuple:
    """Path to a KEGG and description for initial tests."""
    magvsgc: pathlib.Path = pathlib.Path(
        "/Volumes/DIME6/Bin3SB/LOGandSUB/MAGvsGC.txt.gz"
    )
    return magvsgc

def test_hs_kegg_init(kegg_annotation_path):
    """Initialise a KEGG functional annotation from DataFrames."""

    annotation: pd.DataFrame = pd.read_csv(
        kegg_annotation_path.annotation,
        sep="\t",
        header=None,
        index_col=0
    )
    description: pd.DataFrame = pd.read_csv(
        kegg_annotation_path.description,
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
        kegg_annotation_path,
        magvsgc_path
):
    mgs: Dict[str, create.MGS] = create.read_mag_vs_gc(
        magvsgc_path,
        min_completeness=0,
        max_contamination=100
    )

    annotation: pd.DataFrame = pd.read_csv(
        kegg_annotation_path.annotation,
        sep="\t",
        header=None,
        index_col=0
    )
    description: pd.DataFrame = pd.read_csv(
        kegg_annotation_path.description,
        sep="\t",
        header=None,
        index_col=0
    )

    s: HierarchicalSchema = HierarchicalSchema(
        name="KEGG",
        hierarchy=description,
        gene_map=annotation
    )

    tbl: pd.DataFrame = s.function_counts(list(mgs.values()))