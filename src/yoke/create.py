"""Create representations of the gene content of MGS for functional or pangenome
analysis.

Each Metagenomic Species (MGS) is a collection of sufficiently similar
metagenome assembled genomes (MAGs). Each MAG contains many genes, arranged in
order on contigs. We want to analyse and understand the gene content of MGS
in a variety of ways, such as estimate which pathways they contain,
or rearrangements of genes within an MGS.

The package was written primarily to handle output of the MATAFILER pipeline.
This is a gene-centric approach, based on a clustered gene-cataolgue.
Consequently the starting point is at a gene level, representing a MAG as a
collection of contigs, each of which is an ordered list of gene identifiers
from the gene catalogue. Generalisation to other outputs has been kept in
mind during development, but some conversion of data might be necessary to
use the package for other pipelines.
"""

from __future__ import annotations

import itertools
import pathlib
from typing import Dict, List, Optional, Set, Tuple, Union

import pandas as pd

DEF_MAGCVSGC_COLS = dict(
    completeness="Completeness",
    contamination="Contamination",
    representative="Representative4MGS",
    mag="MAG",
    mgs="MGS",
    contigs="other_Genes"
)
"""Column names in the MATAFILER MAGvsGC file."""
# Made this as a dictionary so if output format changes can just change
# column headers here

class MGS:
    """A metagenomic species."""

    def __init__(
        self,
        name: str,
        mags: List[MAG]
    ) -> None:
        """A metagenomic species.

        Args:
            name (str): Name to identify
            mags (List[MAG]): MAGs which make up this MGS
        """
        self.name: str = name
        self.mags: List[MAG] = mags

    @staticmethod
    def from_dataframe_group(
        group: pd.DataFrame
    ) -> MGS:
        """Produce single MGS object from a DataFrame, typically MAGvsGC from
        MG-TK.

        Pass a section of the MAGvsGC DataFrame containing all MAG from one
        MGS. An MGS class object will be created from these MGS.

        Args:
            group (pd.DataFrame): DataFrame with each row being a MAG from
                the same MGS. Should have all the columns in MAGVGC_COLS

        Returns:
            MGS: An MGS containing a MAG for each row
        """
        # Make each MAG from a single row of the DataFrame
        mags: List[MAG] = list(group.apply(MAG.from_series, axis=1))
        name: str = group.iloc[0][DEF_MAGCVSGC_COLS['mgs']]
        return MGS(name=name, mags=mags)

    def __repr__(self) -> str:
        return f'MGS[{self.name}, {len(self.mags)} MAGs]'

    @property
    def genes(self) -> Set[int]:
        """Set of all genes in this MGS."""

        return set().union(*(x.genes for x in self.mags))

    @property
    def edges(self) -> List[Tuple[int, int, Dict[str, int]]]:
        """Edges between genes in this MGS.
        
        The return tuples are a, b, total_count, {MAG: count}
        where a, b are incident vertices of the edge, and MAG: count is a 
        dictionary with the number of times the edge is observed in each MAG.

        Returns:
            List[Tuple[int, int, Dict[str, int]]]: List with tuples in the
                format (a, b, count) where a and b are gene ids, and count
                is a dictionary giving how many time this edge was observed
                in each MAG.

        """

        count_dict: Dict[Tuple[int, int], Dict[str, int]] = dict()
        # Got to be a smarter way of doing this, but this works
        # At least rewrite with defaultdict or something
        for mag in self.mags:
            for edge in mag.edges:
                if edge not in count_dict:
                    count_dict[edge] = dict()
                if mag.name not in count_dict[edge]:
                    count_dict[edge][mag.name] = 0
                count_dict[edge][mag.name] += 1
        # Convert to a list
        return list(
            (e[0], e[1], c) for e, c in count_dict.items()
        )


class MAG:
    """A MAG based on gene catalogue content."""

    def __init__(
        self,
        name: str,
        representative: bool,
        completeness: float,
        contamination: float,
        contigs: List[Contig]
    ) -> None:
        """A MAG based on gene catalogue content.

        Args:
            name (str): Name
            representative (bool): Is this MAG the representative for the MGS?
            completeness (float): Completeness
            contamination (float): Contamination
            contigs (List[Contig]): List of contig objects
        """
        self.name: str = name
        self.representative: bool = representative
        self.completeness: float = completeness
        self.contamination: float = contamination
        self.contigs: List[Contig] = contigs

    @staticmethod
    def from_series(
            mag_series: pd.Series,
            cols: Dict[str, str] = None
    ) -> MAG:
        """_summary_

        Args:
            mag_series (pd.Series): Series containing required fields
            cols (Dict[str, str]): Names of columns in the input

        Returns:
            MAG: MAG based on the series provided
        """
        if cols is None:
            cols = DEF_MAGCVSGC_COLS
        name: str = mag_series[cols['mag']]
        representative: bool = (
                len(mag_series[cols['representative']]) > 0)
        completeness: float = float(
            mag_series[cols['completeness']])
        contamination: float = float(
            mag_series[cols['contamination']])
        contigs: List[Contig] = Contig.from_str(
            mag_series[cols['contigs']],
            name_prefix=f'{name}_')

        return MAG(
            name=name,
            representative=representative,
            completeness=completeness,
            contamination=contamination,
            contigs=contigs
        )

    @property
    def genes(self) -> Set[int]:
        """Set of all genes in this MAG."""
        return set(itertools.chain.from_iterable(
            x.genes for x in self.contigs
        ))

    @property
    def edges(self) -> List[Tuple[int, int]]:
        """All edges in this MAG.

        Returns:
            List[Tuple[int, int]]: List of all edges as a tuple
        """

        edges: List[Tuple[int, int]] = list(itertools.chain.from_iterable(
            x.edges for x in self.contigs
        ))
        return edges


    def __repr__(self) -> str:
        return f'MAG[{self.name}, {len(self.contigs)} contigs]'


class Contig:
    """One contig, which is an ordered list of genes."""

    def __init__(
        self,
        genes: List[int],
        name: Optional[str] = None
    ) -> None:
        """Representation of a single contig.

        Args:
            genes (List[int]): Ordered list of genes on contig
            name (Optional[str], optional): Name of contig. Defaults to None.
        """

        self.genes: List[int] = [int(i) for i in genes if i != "?"]
        self.name: str = "" if name is None else name

    @staticmethod
    def from_str(
        contig_str: str,
        gene_delim: str = ",",
        contig_delim: str = ",,",
        name_prefix: str = ""
    ) -> List[Contig]:
        """Create Contig from string representation of sequence of genes.

        Each contig will be given a number, prefixed with name_prefix.

        Args:
            contig_str (str): String giving the order of genes on all contigs.
            gene_delim (str, optional): Character which delimits genes. Defaults
                to ",".
            contig_delim (str, optional): String which indicates the end of a
                contig. Defaults to ",,".
            name_prefix (str, optional): String to place before contig number. 
                Defaults to "".

        Returns:
            List[Contig]: _description_
        """
        # Split to contigs
        contigs: List[Contig] = [
            Contig(ctg.split(gene_delim), name=f'{name_prefix}{i}')
            for i, ctg in enumerate(contig_str.split(contig_delim))
        ]
        return contigs

    def __repr__(self) -> str:
        return f'Contig[{self.name}, {len(self.genes)} genes]'

    @property
    def edges(self) -> List[Tuple[int, int]]:
        """All the edges supported by this contig (consecutive genes).
        

        Returns:
            List[Tuple[int, int]]: List of adjacent genes on this contig.
        """

        return [
            (self.genes[i], self.genes[i+1]) for i in range(len(self.genes)-1)
        ]


def read_mag_vs_gc(
    file: Union[str, pathlib.Path],
    min_completeness: float = 0,
    max_contamination: float = 100,
    cols: Dict[str, str] = None,
    **kwargs
) -> Dict[str, MGS]:
    """Create MGS from file connecting MAGs to gene catalogue.

    The MAGvsGC file in MG-TK provides an ordered lists of genes found
    on each contig in MAGs, along with metdata on each MAG and MGS. This
    will be parsed to make a set of objects represent all these elements, which
    can be used to create a pangenome graph and other analyses.

    :param file: Path to MAGvsGC
    :param min_completeness: Exclude any MAGs below minimum completeness
    :param max_contamination: Exclude any MAGs above this contamination
    :param cols: Column names in the input file
    :param **kwargs: Passed to pandas read_csv
    """

    if cols is None:
        cols = DEF_MAGCVSGC_COLS

    df: pd.DataFrame = pd.read_csv(
        file,
        sep="\t",
        low_memory=True
    )

    # Filter MAGs
    df = (df[(df[cols['completeness']] > min_completeness) &
             (df[cols['contamination']] < max_contamination)]
            .astype({cols['representative']: "str"})
        )

    # Group by MGS
    mgs: List[MGS] = list(
        df.groupby(by=cols['mgs']).
        apply(MGS.from_dataframe_group)
    )

    return dict(zip((x.name for x in mgs), mgs))