"""Link between gene IDs and functional annotations

Genes may be annotated against several functional schema, such as KEGG
Orthologs, eggNOGs, GBM modules etc. This module provides methods to handle
the connection between gene identities and these annotation schema.
"""
from __future__ import annotations

import csv
import dataclasses
import gzip
import itertools
import logging
import pathlib
import sys
from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Any, Iterable, List, Set, Dict, Union

import numpy as np
import pandas as pd
import tqdm
from multimethod import multimethod
from pandas.core.computation.expr import intersection

from yoke import create
from yoke.create import MAG, MGS

logger: logging.Logger = logging.getLogger(__name__)


@dataclass
class GeneFunction:
    """Represents one function and gene IDs which have been annotated with
    it.

    Attributes:
        id: Function ID.
        desc: Description of function.
        color: Color to use when representing this gene.
    """
    id: str
    desc: Optional[str] = None
    color: Optional[str] = None
    gene_ids: Set[int] = dataclasses.field(default_factory=set)

    @property
    def sparse_series(self) -> pd.Series:
        a = pd.Series(
            [1] * len(self.gene_ids),
            dtype=pd.SparseDtype(np.dtype("int_"), fill_value=0),
        )
        return a

class FunctionSchema(Enum):
    """Functional annotation schema supported by default load methods"""
    KEGG = 1


class HierarchicalSchema:
    """Functional annotation schema which has one or more levels of
    hierarchy.

    Hierarchical annotation schema can be represented by a tree, allowing
    for multifurcation. Examples include GBM modules or BRITE Hierarchies."""

    def __init__(
            self,
            name: str,
            hierarchy: pd.DataFrame,
            gene_map: pd.DataFrame,
            **kwargs
    ):
        """Initialise a hierarchical annotation schema.

        Args:
            name: Name of the annotation schema.
            hierarchy_table: Table where the index is function IDs, and each
            column represents a level of the hierarchy.
            gene_map: Table with the index being gene ID, and the only column
            being function ID.
            **kwargs: Can add progress_bar=False to disable progress bar
        """
        self.__name: str = name
        self.__hierarchy: pd.DataFrame = hierarchy
        self.__gene_map: Dict = dict(zip(
            gene_map.index, gene_map.iloc[:, 0]
        ))
        self.__functions: Dict[str, GeneFunction] = (
            self.__make_gene_functions(hierarchy, gene_map)
        )

    @property
    def name(self) -> str:
        return self.__name

    @property
    def hierarchy(self) -> pd.DataFrame:
        return self.__hierarchy

    @property
    def functions(self) -> Dict[str, GeneFunction]:
        return self.__functions

    @property
    def gene_map(self) -> pd.DataFrame:
        return self.__gene_map

    @property
    def gene_matrix(self) -> pd.DataFrame:
        mat: pd.DataFrame = pd.concat(
            [f.sparse_series for f in self.functions.values()],
            axis=1
        )
        return mat

    def __make_hierarchy(
            self,
            hierarchy: Optional[pd.DataFrame],
            annotation: pd.DataFrame
    ) -> pd.DataFrame:
        """Properly format a hierarchy table, or make a dummy if is no table."""
        if hierarchy is not None:
            return hierarchy
        # Make a single level hierarchy with all functions from the annotation
        funcs: pd.Series = annotation.iloc[:, 0].unique()

    def get_gene_function(self, gene_id) -> Optional[GeneFunction]:
        """Get the annotation for a given gene. Returns None if no
        annotation."""
        if gene_id not in self.gene_map:
            return None
        function_id = self.gene_map[gene_id]
        return self.functions[function_id]


    @multimethod
    def function_counts(
            self,
            mag: MAG
    ):
        """Count number of genes in this MAG with each function in this
        schema.

        Args:
            mag: MAG to produce counts for.
        """
        f_counts: Counter = Counter(
            (
                gf.id for gf in
                (self.get_gene_function(gene) for gene in mag.genes)
                if gf is not None
            )
        )
        counts: pd.Series = pd.Series(
            f_counts,
            name=mag.name
        )
        return counts

    @multimethod
    def function_counts(
            self,
            mgs: MGS
    ):
        """Count number of genes in each MAG in this MGS with each function
        in this schema.

        Args:
            mgs: MGS to produce counts for.
        """

        # Make a multiindex with MAG and MGS - MGS is redundant, but it
        # makes columnwise concating easier
        cols: pd.MultiIndex = pd.MultiIndex.from_tuples(
            zip(
                [mag.name for mag in mgs.mags],
                itertools.cycle([mgs.name])
            ), names=["mag", "mgs"]
        )
        df = pd.concat(
            [self.function_counts(mag) for mag in
             tqdm.tqdm(mgs.mags, desc=mgs.name)],
            axis=1
        )
        df.columns = cols
        return df

    @multimethod
    def function_counts(
            self,
            mgs: List[MGS]
    ) -> pd.DataFrame:
        """Count number of genes in each MAG in these MGS with each functional
        in this schema.

        Args:
            mgs: Collection of MGS to produce counts for.
        """
        return pd.concat(
            [self.function_counts(m) for m in mgs],
            axis=1
        )

    @staticmethod
    def __make_gene_functions(
            hierarchy: pd.DataFrame,
            gene_map: pd.DataFrame
    ) -> Dict[str, GeneFunction]:
        """Create GeneFunction objects for this hierarchy."""

        # Refactor for speed
        # Try group by apply
        ret_dict: Dict[str, GeneFunction] = dict()
        for function_id, ann_group in gene_map.groupby(by=1):
            ret_dict[function_id] = GeneFunction(
                id=function_id,
                desc=(
                    None if function_id not in hierarchy.index else
                    hierarchy.loc[function_id].iloc[-1]
                ),
                gene_ids=set(ann_group.index)
            )
        return ret_dict

    @staticmethod
    def from_file(
            name: str,
            hierarchy: Union[str, pathlib.Path],
            annotation: Union[str, pathlib.Path],
            delim: str = "\t"
    ) -> HierarchicalSchema:
        """Read schema from files.

        Create an annotation schema from a definition of the hierarchy,
        and an annotation file which connects genes and function IDS.
        Annotation files should be only two columns, with no headers.
        """

        h_df: pd.DataFrame = pd.read_csv(
            hierarchy,
            sep=delim,
            index_col=0
        )
        a_df: pd.DataFrame = pd.read_csv(
            annotation,
            sep=delim,
            index_col=0
        )

        return HierarchicalSchema(
            name=name,
            hierarchy=h_df,
            gene_map=a_df
        )

def _combine_mgs(
        mag_counts: pathlib.Path,
        mgs: List[str]
) -> pd.DataFrame:
    """Combine function counts of multiple MGS into one table. The returned
    table will have one column for each MAG in any of the listed MGS. This
    can be used to get a table of all MGS with a certain phenotype, or from
    a given genus etc.
    """
    mgs_dfs: List[pd.DataFrame] = [
        pd.read_csv(mag_counts / f"{x}.tsv", sep="\t", index_col=0) for x in mgs
    ]
    comb_df: pd.DataFrame = pd.concat(mgs_dfs, axis=1)
    return comb_df

def _csv_field_adjust():
    """Increase max field size for csv to help parse larger entries in
    MAGvsGC"""
    logger.debug("Increasing max field size to %s", str(sys.maxsize))
    logger.debug("This has been needed in the past to parse large entries in "
                 "MAGvsGC")
    csv.field_size_limit(sys.maxsize)


def cli_function_count(
        magvsgc: pathlib.Path,
        annotations: pathlib.Path,
        descriptions: Optional[pathlib.Path],
        output: pathlib.Path,
        name: str,
        completeness: float = 0,
        contamination: float = 100,
        compress: bool = True
) -> None:
    """Count how many genes have each functional annotation in all MAGS."""

    _csv_field_adjust()
    logger.setLevel(logging.INFO)

    logger.info("function_count")
    logger.info("|- Parameters -|")
    logger.info("MAGvsGC file: %s", magvsgc)
    logger.info("Gene annotations: %s", annotations)
    logger.info("Annotation descriptions: %s", descriptions)
    logger.info("Output: %s", output)
    logger.info("Name: %s", name)
    logger.info("Completeness: %s", completeness)
    logger.info("Contamination: %s", contamination)
    logger.info("Compress: %s", compress)

    pth_out: pathlib.Path = output / name
    pth_out.mkdir(parents=False, exist_ok=True)

    logger.info("Reading annotation descriptions")
    description: pd.DataFrame = pd.read_csv(
        descriptions,
        sep="\t",
        header=None,
        index_col=0
    )
    logger.info("Reading gene annotations")
    annotation: pd.DataFrame = pd.read_csv(
        annotations,
        sep="\t",
        header=None,
        index_col=0
    )
    logger.info("Making annotation object")
    s: HierarchicalSchema = HierarchicalSchema(
        name=name,
        hierarchy=description,
        gene_map=annotation
    )

    logger.info("Creating function counts for each MGS")
    mgs_gen = create.generate_mgs_from_magvsgc(
        magvsgc,
        min_completeness=completeness,
        max_contamination=contamination,
        compressed=True
    )

    # Collect means to make a summary DF also
    means: List[pd.Series] = []
    comp_args: Dict = dict(compression="gzip") if compress else dict()
    for mgs in mgs_gen:
        pth_mgs: pathlib.Path = (
                pth_out / f"{mgs.name}.tsv.gz"
                if compress else
                pth_out / f"{mgs.name}.tsv"
        )
        tbl: pd.DataFrame = s.function_counts(mgs)
        tbl.columns = tbl.columns.droplevel("mgs")
        logger.debug("Write to %s", pth_mgs)
        tbl.to_csv(pth_mgs, sep="\t", **comp_args)
        mean: pd.Series = tbl.mean(axis=1)
        mean.index.name = mgs.name
        means.append(mean)
    pth_all: pathlib.Path = output / f"{name}.tsv"
    logger.debug("Producing per MGS mean table")
    df_means: pd.DataFrame = pd.concat(means, axis=1)
    logger.info("MGS means saved to %s", pth_all)
    df_means.to_csv(pth_all, sep="\t")

def cli_function_hoist(
        mag_counts: pathlib.Path,
        magvsgc_file: pathlib.Path,
        output: pathlib.Path,
        ranks: List[str] = ['Domain', 'Phylum', 'Class', 'Order', 'Family',
                            'Genus', 'Species'],
        compressed: bool = True,
        cols: dict[str, str] = None
) -> None:
    """Summarise previously produced function counts per MAG at different
    taxonomic levels."""

    _csv_field_adjust()
    logger.setLevel(logging.INFO)

    logger.info("function_summarise")
    logger.info("|- Parameters -|")
    logger.info("MAG function counts: %s", mag_counts)
    logger.info("MAGvsGC file: %s", magvsgc_file)
    logger.info("Ranks: %s", ",".join(ranks))
    logger.info("Output: %s", output)
    logger.info("Compress: %s", compressed)

    # Take per MGS taxonomy from MAGvsGC file
    if cols is None:
        cols = create.DEF_MAGCVSGC_COLS

    logger.info("Reading MAGvsGC file %s", magvsgc_file)
    openfn = (lambda x: gzip.open(x, 'rt')) if compressed else open
    with openfn(magvsgc_file) as magvsgc:
        reader = csv.DictReader(magvsgc, delimiter="\t")
        mgs_tax: Dict[str, pd.Series] = dict()
        for row in reader:
            mgs_id: str = row[cols['mgs']]
            if (mgs_id not in mgs_tax) and (len(mgs_id) > 0):
                tdict = {r: row[r] for r in ranks}
                tdict['mgs'] = mgs_id
                mgs_tax[mgs_id] = tdict
    tax_df: pd.DataFrame = pd.DataFrame(
        [pd.Series(x) for x in mgs_tax.values()]
    )

    # Group MGS by each level of each requested rank
    for i, rank in enumerate(ranks):
        logger.info("Summarise to rank %s", rank)
        pth_rank = output / rank
        pth_rank.mkdir(parents=True, exist_ok=True)
        means: List[pd.Series] = []
        tax: List[str] = list(tax_df[rank].unique())
        for val in (pbar := tqdm.tqdm(tax)):
            pbar.set_description(f"{rank}: {val}")
            pth_full = pth_rank / f"{val}.tsv"
            df = _combine_mgs(
                mag_counts,
                tax_df.loc[tax_df[rank] == val, 'mgs'].tolist()
            )
            logger.debug("Write to %s", pth_full)
            df.to_csv(pth_full, sep="\t")
            logger.debug("Calculate mean for taxon")
            mean_series = df.mean(axis=1)
            mean_series.name = val
            means.append(mean_series)
        logger.info("Writing taxon mean counts")
        df_mean = pd.DataFrame(means)
        df_mean.T.to_csv(output / f"{rank}.tsv", sep="\t")