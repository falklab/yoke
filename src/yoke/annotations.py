"""Link between gene IDs and functional annotations

Genes may be annotated against several functional schema, such as KEGG
Orthologs, eggNOGs, GBM modules etc. This module provides methods to handle
the connection between gene identities and these annotation schema.
"""
from __future__ import annotations
import dataclasses
import itertools
import pathlib
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Any, Iterable, List, Set, Dict, Union

import numpy as np
import pandas as pd
import tqdm
from multimethod import multimethod

from create import MAG, MGS


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
        self.__gene_map: pd.DataFrame = gene_map
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
        counts: pd.Series = pd.Series(
            [len(mag.genes.intersection(f.gene_ids)) for f in self.functions],
            index=mag.name
        )

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
                [mag.id for mag in mgs.mags],
                itertools.cycle([mgs.id])
            ), names=["mag", "mgs"]
        )
        df = pd.concat(
            [self.function_counts(mag) for mag in mgs.mags],
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
            gene_map: pd.DataFrame,
            progress_bar: bool = True
    ) -> Dict[str, GeneFunction]:
        """Create GeneFunction objects for this hierarchy."""

        # TODO: Refactor, this is too slow
        # Subset to only the annotations which are in both the annotation and
        # the description
        intersection: np.ndarray = (
            np.intersect1d(
                np.unique(gene_map.iloc[:,-1]),
                hierarchy.index.unique()
            )
        )
        intersection_iterator = tqdm.tqdm(
            intersection,
            desc="Link annotations to genes",
            disable=not progress_bar
        )

        # Working with ndarray objects only rather than DataFrames for speed
        map_arr: np.ndarray = gene_map.reset_index().values
        return {
            function_id: GeneFunction(
                id=function_id,
                desc=hierarchy.loc[function_id].iloc[-1],
                gene_ids=set(map_arr[map_arr[:, 1] == function_id, 0])
            )
            for function_id in intersection_iterator
        }



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