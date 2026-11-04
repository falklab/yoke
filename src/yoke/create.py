"""Represent the gene content of MGS for functional or pangenome analysis.

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

import csv
import gzip
import itertools
import pathlib
from collections import Counter
from collections.abc import Iterator

import numpy as np
import pandas as pd

DEF_MAGCVSGC_COLS = dict(
    completeness="Completeness",
    contamination="Contamination",
    representative="Representative4MGS",
    mag="MAG",
    mgs="MGS",
    contigs="other_Genes",
)
"""Column names in the MATAFILER MAGvsGC file."""
# Made this as a dictionary so if output format changes can just change
# column headers here


class MGS:
    """A metagenomic species."""

    def __init__(self, name: str, mags: list[MAG]) -> None:
        """A metagenomic species.

        Args:
            name (str): Name to identify
            mags (List[MAG]): MAGs which make up this MGS
        """
        self.name: str = name
        self.mags: list[MAG] = mags

    @staticmethod
    def from_dataframe_group(group: pd.DataFrame) -> MGS:
        """Produce one MGS object from a DataFrame, typically MG-TK MAGvsGC.

        Pass a section of the MAGvsGC DataFrame containing all MAG from one
        MGS. An MGS class object will be created from these MGS.

        Args:
            group (pd.DataFrame): DataFrame with each row being a MAG from
                the same MGS. Should have all the columns in MAGVGC_COLS

        Returns:
            MGS: An MGS containing a MAG for each row
        """
        # Make each MAG from a single row of the DataFrame
        mags: list[MAG] = list(group.apply(MAG.from_series, axis=1))
        name: str = group.iloc[0][DEF_MAGCVSGC_COLS["mgs"]]
        return MGS(name=name, mags=mags)

    def __repr__(self) -> str:
        return f"MGS[{self.name}, {len(self.mags)} MAGs]"

    @property
    def genes(self) -> set[int]:
        """Set of all genes in this MGS."""
        return set().union(*(x.genes for x in self.mags))

    @property
    def gene_matrix(self) -> pd.DataFrame:
        """Matrix of gene count for each MAG in the MGS.

        The matrix will have genes on rows, with MAGs on columns. Each entry
        will be the number of times the gene occurs.
        """
        return pd.concat([mag.gene_series for mag in self.mags])

    @property
    def edges(
        self,
    ) -> list[tuple[int, int, dict[str, int], dict[str, set[Contig]]]]:
        """Edges between genes in this MGS.

        The return tuples are a, b, {MAG: count}, {MAG: contigs}
        where a, b are incident vertices of the edge, and MAG: count is a
        dictionary with the number of times the edge is observed in each MAG.

        Returns:
            List[Tuple[int, int, Dict[str, int], Dict[str, Contig]]: List with
            tuples in the format (a, b, count) where a and b are gene ids, and
            count is a dictionary giving how many times this edge was observed
            in each MAG.
        """
        count_dict: dict[tuple[int, int], dict[str, int]] = dict()
        contig_dict: dict[tuple[int, int], dict[str, set[Contig]]] = dict()
        # Got to be a smarter way of doing this, but this works
        # At least rewrite with defaultdict or something
        for mag in self.mags:
            for x in mag.edges:
                edge = x[0], x[1]
                if edge not in count_dict:
                    count_dict[edge] = dict()
                    contig_dict[edge] = dict()
                if mag.name not in count_dict[edge]:
                    count_dict[edge][mag.name] = 0
                    contig_dict[edge][mag.name] = set()
                count_dict[edge][mag.name] += 1
                contig_dict[edge][mag.name].add(x[-1])
        # Convert to a list
        return list(
            (e[0], e[1], c, contig_dict[e]) for e, c in count_dict.items()
        )


class MAG:
    """A MAG based on gene catalogue content."""

    def __init__(
        self,
        name: str,
        representative: bool,
        completeness: float,
        contamination: float,
        contigs: list[Contig],
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
        self.contigs: list[Contig] = contigs

    @staticmethod
    def from_series(
        mag_series: pd.Series, cols: dict[str, str] | None = None
    ) -> MAG:
        """Create a MAG from one row of a MAGvsGC table.

        Args:
            mag_series (pd.Series): Series containing required fields
            cols (Dict[str, str]): Names of columns in the input

        Returns:
            MAG: MAG based on the series provided
        """
        if cols is None:
            cols = DEF_MAGCVSGC_COLS
        name: str = mag_series[cols["mag"]]
        representative: bool = len(mag_series[cols["representative"]]) > 0
        completeness: float = float(mag_series[cols["completeness"]])
        contamination: float = float(mag_series[cols["contamination"]])
        contigs: list[Contig] = Contig.from_str(
            mag_series[cols["contigs"]], name_prefix=f"{name}_"
        )

        return MAG(
            name=name,
            representative=representative,
            completeness=completeness,
            contamination=contamination,
            contigs=contigs,
        )

    @property
    def genes(self) -> set[int]:
        """Set of all genes in this MAG."""
        return set(
            itertools.chain.from_iterable(x.genes for x in self.contigs)
        )

    @property
    def gene_series(self) -> pd.Series:
        """Series with count of observed genes."""
        gene_count: Counter = Counter(self.genes)
        return pd.Series(gene_count, name=self.name)

    @property
    def edges(self) -> list[tuple[int, int, Contig]]:
        """All edges in this MAG.

        Return every edge implied by any contig in this MAG. This will include
        duplicate edges where the same edges appears in the MAG multiple times,
        or even where the same edges appears on a contig multiple times.

        Returns:
            List[Tuple[int, int, Contig]]: List of all edges as a tuple
        """
        edges: list[tuple[int, int, Contig]] = [
            (x[0], x[1], contig)
            for contig in self.contigs
            for x in contig.edges
        ]
        return edges

    def __repr__(self) -> str:
        return f"MAG[{self.name}, {len(self.contigs)} contigs]"


class Contig:
    """One contig, which is an ordered list of genes."""

    def __init__(self, genes: list[int], name: str | None = None) -> None:
        """Representation of a single contig.

        Args:
            genes (List[int]): Ordered list of genes on contig
            name (Optional[str], optional): Name of contig. Defaults to None.
        """
        self.genes: list[int] = [int(i) for i in genes if i != "?"]
        self.name: str = "" if name is None else name

    @staticmethod
    def from_str(
        contig_str: str,
        gene_delim: str = ",",
        contig_delim: str = ",,",
        name_prefix: str = "",
    ) -> list[Contig]:
        """Create Contig from string representation of sequence of genes.

        Each contig will be given a number, prefixed with name_prefix.

        Args:
            contig_str (str): String giving the order of genes on all contigs.
            gene_delim (str, optional): Character which delimits genes.
                Defaults to ",".
            contig_delim (str, optional): String which indicates the end of a
                contig. Defaults to ",,".
            name_prefix (str, optional): String to place before contig number.
                Defaults to "".

        Returns:
            List[Contig]: _description_
        """
        # Split to contigs
        contigs: list[Contig] = [
            Contig(ctg.split(gene_delim), name=f"{name_prefix}{i}")
            for i, ctg in enumerate(contig_str.split(contig_delim))
        ]
        return contigs

    def __repr__(self) -> str:
        return f"Contig[{self.name}, {len(self.genes)} genes]"

    @staticmethod
    def reverse(contig: Contig) -> Contig:
        """Copy a contig with gene order reversed.

        Args:
            contig: Contig to be reversed
        Returns:
            Contig: Copy of contig with gene order reversed. The contig
            name will have "_r" appended to indicate this has been reversed.
        """
        rev_genes = [int(x) for x in contig.genes[::-1]]
        return Contig(genes=rev_genes, name=f"{contig.name}_r")

    @property
    def edges(self) -> list[tuple[int, int]]:
        """All the edges supported by this contig (consecutive genes).

        Returns:
            List[Tuple[int, int]]: List of adjacent genes on this contig.
        """
        return [
            (self.genes[i], self.genes[i + 1])
            for i in range(len(self.genes) - 1)
        ]

    def get_indices(self, genes: list[int]) -> np.ndarray:
        """Positions on contig of a given sequence of genes.

        Args:
            genes: Sequence of genes identifiers to be found

        Returns:
            Optional[List[List[int]]: List where each item is the list of
                sequential indices for an occurrence of the pattern of genes.
                The indices correspond to the order of the genes in the input
                pattern. Returns None when the pattern is not found.
        """
        contig_genes: np.ndarray = np.array(self.genes)
        pattern: np.ndarray = np.array(genes)
        # Sliding window is inefficient, potentially an area to optimise
        # if running slow - but realistically I don't expect this function
        # to be the workhorse of the package
        window: np.ndarray = np.lib.stride_tricks.sliding_window_view(
            contig_genes, pattern.shape[0]
        )

        forward_starts: np.ndarray = (
            np.where(np.all(window == pattern, axis=1))
        )[0]
        forward_shifts = np.tile(
            np.arange(pattern.shape[0]), (forward_starts.shape[0], 1)
        )
        forward_indices = forward_shifts + forward_starts.reshape((-1, 1))

        # Palindromic patterns do not need reverse matches calculating
        if not np.array_equal(genes, genes[::-1]):
            reverse_starts: np.ndarray = np.where(
                np.all(window == pattern[::-1], axis=1)
            )[0]
            reverse_shifts = np.tile(
                np.arange(pattern.shape[0] - 1, -1, -1),
                (reverse_starts.shape[0], 1),
            )
            reverse_indices = reverse_shifts + reverse_starts.reshape((-1, 1))
            return np.concatenate([forward_indices, reverse_indices])
        else:
            return forward_indices

    @staticmethod
    def get_contexts(
        contigs: list[Contig], genes: list[int], flanking: int
    ) -> pd.DataFrame:
        """Get a sequence of genes plus flanking genes on multiple contigs.

        Get the genomic context of a pattern of genes, on multiple contigs.
        See `Contig.get_context()` for details of arguments and method.
        """
        ctxts: pd.DataFrame = pd.concat(
            [
                pd.DataFrame(c.get_context(genes, flanking)).assign(
                    contig=c.name
                )
                for c in contigs
            ]
        )
        return ctxts

    def get_context(self, genes: list[int], flanking: int) -> np.ndarray:
        """Get a sequence of genes plus some up or downstream.

        Finds each occurrence of a given pattern of genes, and returns
        those genes with a number of flanking genes either side. This allows
        a view of the occurrence of the genes in their genomic context. The
        query genes will be placed in the centre of the return sequence. If
        there are no genes up or down stream (due to end of contig) those
        positions will be filled with NaN.

        Args:
            genes: List with sequence of genes to be identified
            flanking: Number of genes either side to return

        Returns:
            Array where is row is an occurrence of the pattern, and each entry
            a gene ID. The central columns in the matrix are the pattern.
        """
        indices: np.ndarray = self.get_indices(genes=genes)
        coll = []
        for m in indices:
            start = min(m) - flanking
            end = max(m) + flanking + 1
            genes_nums = [
                self.genes[i] if (i >= 0) and (i < len(self.genes)) else None
                for i in range(start, end)
            ]
            coll.append(genes_nums)
        return np.array(coll)


def read_mag_vs_gc(
    file: str | pathlib.Path,
    min_completeness: float = 0,
    max_contamination: float = 100,
    cols: dict[str, str] | None = None,
    **kwargs,
) -> dict[str, MGS]:
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

    df: pd.DataFrame = pd.read_csv(file, sep="\t", low_memory=True)

    # Filter MAGs
    df = df[
        (df[cols["completeness"]] > min_completeness)
        & (df[cols["contamination"]] < max_contamination)
    ].astype({cols["representative"]: "str"})

    # Group by MGS
    mgs: list[MGS] = list(
        df.groupby(by=cols["mgs"]).apply(MGS.from_dataframe_group)
    )

    return dict(zip((x.name for x in mgs), mgs, strict=True))


def generate_mgs_from_magvsgc(
    file: str | pathlib.Path,
    compressed: bool = True,
    min_completeness: float = 0,
    max_contamination: float = 100,
    cols: dict[str, str] | None = None,
    **kwargs,
) -> Iterator[MGS]:
    """Yield MGS objects based on MAGvsGC content.

    Compared to read_mag_vs_gc, this avoids reading the whole table into
    memory, and yields individual MGS objects rather than collecting them.
    """
    if cols is None:
        cols = DEF_MAGCVSGC_COLS

    openfn = (lambda x: gzip.open(x, "rt")) if compressed else open
    with openfn(file) as magvsgc:
        reader = csv.DictReader(magvsgc, delimiter="\t", **kwargs)
        mags: list[dict] = []
        prev_mgs: str | None = None
        for row in reader:
            mgs: str = row[cols["mgs"]]
            if (prev_mgs != mgs) and (len(mags) > 0):
                # Create and yield an MGS
                df: pd.DataFrame = pd.DataFrame(mags).astype(
                    {
                        cols["representative"]: "str",
                        cols["completeness"]: "float",
                        cols["contamination"]: "float",
                    }
                )
                # Filter MAGs
                df = df[
                    (df[cols["completeness"]] > min_completeness)
                    & (df[cols["contamination"]] < max_contamination)
                ].astype({cols["representative"]: "str"})
                yield MGS.from_dataframe_group(df)
                mags = []

            # If this MAG is not part of an MGS skip it
            if len(mgs.strip()) < 1:
                continue

            # Select only requested columns
            mags.append({v: row[v] for k, v in DEF_MAGCVSGC_COLS.items()})

            prev_mgs = mgs

        # Yield final MGS
        if len(mags) > 0:
            # TODO: Tidy and remove this duplication
            df: pd.DataFrame = pd.DataFrame(mags).astype(
                {
                    cols["representative"]: "str",
                    cols["completeness"]: "float",
                    cols["contamination"]: "float",
                }
            )
            # Filter MAGs
            df = df[
                (df[cols["completeness"]] > min_completeness)
                & (df[cols["contamination"]] < max_contamination)
            ].astype({cols["representative"]: "str"})
            if df.shape[0] > 0:
                yield MGS.from_dataframe_group(df)
