"""Analysis of gene-level pangenome graphs of MGS.

Provides functions to create and analyse pangenome graphs for MGS objects.
"""
from __future__ import annotations
from typing import List, Set, Any, NamedTuple, Tuple, Dict

from multimethod import multimethod
import networkx as nx
import pandas as pd
import plotnine as pn

from yoke import MGS, MAG, Contig


class MGSGraph:
    """Graph representation of an MGS.

    A pangene graph representing the MAGs in an MGS.
    """

    def __init__(
            self,
            mgs: MGS
    ) -> None:
        self.__mgs = mgs
        self.__graph = self.__make_graph()

    @property
    def mgs(self) -> MGS:
        return self.__mgs

    @property
    def graph(self) -> nx.Graph:
        return self.__graph

    def __make_graph(self) -> nx.Graph:
        """Create a pangenome graph for this MGS.

        Returns:
            nx.Graph: Undirected pangenome graph. Each edge has the following
                properties:
                * weight: the number of times that edge occurs in all MAGs
                * mag_count: dictionary with number of times the edge found
                    per MAG
                * mag_set: the MAGs in which this edge was found
        """
        graph: nx.Graph = nx.Graph()
        # Add all vertices
        graph.add_nodes_from(self.mgs.genes)
        # Add all edges
        for a, b, count, contig in self.mgs.edges:
            graph.add_edge(a, b,
                           weight=sum(count.values()),
                           contigs=contig,
                           mag_count=count,
                           mag_set=count.keys())
        return graph


    def hg_paths(self) -> List[List[int]]:
        """Find all homogenomic paths in the graph.

        A path is homogenomic when the edges between each pair of points on
        the path have the same set support set of MAGs. Stated in terms of
        biology, a homogenomic path represents a sequence of genes which occur
        in the same order in a subset of MAGs for this MGS.
        """

        # All homogenomic stretches will necessarily be included in a spanning
        # tree of the graph. Currently, I don't think the nature of the spanning
        # tree matters too much (i.e it need not be the minimum).

        st: nx.Graph = nx.minimum_spanning_tree(pg)
        # Traverse tree building homogenomic paths as we go
        hg_paths: List[List[int]] = []
        curr_path: List[int] = []
        for v in nx.dfs_preorder_nodes(st):
            degree: int = st.degree(v)
            if degree == 0:
                curr_path = []
                continue
            mag_sets: List[Set] = [
                set(x['mag_set']) for x in st.adj[v].values()
            ]
            if degree < 3 and set.union(*mag_sets) == mag_sets[0]:
                curr_path.append(v)
            elif len(curr_path) > 1:
                hg_paths.append(curr_path)
                curr_path = []
            else:
                curr_path = []
        return hg_paths

    def tandem_duplication(self) -> List[TandemIndel]:
        """Find tandem duplications of genes.

        Tandem genes are when the same gene occurs next to each other without
        intervening genes. A tandem duplication is where some MAGs have
        a tandem gene where others have a single copy.
        """

        v_self = nx.nodes_with_selfloops(self.__graph)
        # If the weight of the self loop is lower than the other edges, then
        # this is a duplication.
        res: List[TandemIndel] = []
        for v in v_self:
            adj = self.graph.adj[v]
            # If 3 adjacent vertices, this is the simplest case we could hope
            # for
            if len(adj) == 3:
                dup_set: Set = set(adj[v]['mag_set'])
                neighbours = list(adj.keys())
                neighbours.remove(v)
                l_set = set(adj[neighbours[0]]['mag_set'])
                r_set = set(adj[neighbours[1]]['mag_set'])
                # Being very strict with criteria
                # neighbouring bits must be all the same taxa, loop must be a subset
                if ((l_set == r_set) and (dup_set.issubset(l_set)) and
                    (len(dup_set) < len(l_set))):
                    support_contigs: Dict[str, Set[Contig]] = (
                        merge_contig_dicts(
                            [adj[v]['contigs'] for v in adj]
                        )
                    )
                    res.append(TandemIndel(
                        gene=v,
                        vertex=v,
                        duplicated_in=dup_set,
                        single_in=l_set - dup_set,
                        contigs=support_contigs
                    ))
        return res

class TandemIndel(NamedTuple):
    gene: int
    """Gene which has a tandem duplication."""
    vertex: int
    """Vertex ID which has duplication."""
    duplicated_in: Set[str]
    """Genomes the duplication appears in."""
    single_in: Set[str]
    """Genomes which do not have the duplication."""
    contigs: Dict[str, Set[Contig]]
    """Contigs in each genome which support this variant. Key is MAG name,
    and value a set of Contigs."""

    @property
    def id(self) -> str:
        """ID for this variant. Unique only within the MGS."""
        return f"tr_{self.vertex}"

    @property
    def size(self) -> Tuple[int, int]:
        """Number of MAGs which have duplicate and single respectively."""

        return len(self.duplicated_in), len(self.single_in)

    @property
    def ratio(self) -> float:
        """Ratio of duplicate to single MAGs."""

        d, s = self.size
        return d / s

    def plot(
            self,
            window: int = 5
    ) -> pn.ggplot:
        """Plot showing the region of contigs containing the variant.

        Args:
            window: Number of genes up and downstream of the variant region
            to display
        """
        # Find variant regions in all contigs, and trim to window length
        duplicated_df = pd.concat(
            [
                Contig.get_contexts(
                    contigs=list(self.contigs[d]),
                    genes=[self.gene, self.gene],
                    flanking=window
                ).assign(
                    mag=d,
                    duplicated=True
                ) for d in self.duplicated_in
            ]
        )
        single_df = pd.concat(
            [
                Contig.get_contexts(
                    contigs=list(self.contigs[s]),
                    genes=[self.gene],
                    flanking=window
                ).assign(
                    mag=s,
                    duplicated=False
                ) for s in self.single_in
            ]
        )
        # Dummy empty column in single
        single_df = single_df.rename(columns=dict(zip(
            range(window + 1, (2 * window) + 2),
            range(window + 2, (2 * window) + 3)
        )))
        single_df.insert(
            window + 1,
            window + 1,
            pd.Series([None] * single_df.shape[0])
        )
        df = pd.concat([duplicated_df, single_df])

        # Convert to a format plottable by plotnine
        dfs = (
            df
            .set_index(["mag", "contig", "duplicated"])
            .stack()
            .to_frame("gene_id")
            .reset_index()
            .rename(columns=dict(level_3="position"))
        )

        contig_count: int = sum(len(x) for x in self.contigs.values())
        plt: pn.ggplot = (
            pn.ggplot(
                dfs,
                pn.aes(
                    x="position",
                    y="contig",
                    color="gene_id"
                )
            )
            + pn.geom_point()
            + pn.geom_text(
                mapping=pn.aes(label="gene_id"),
                nudge_y=0.5,
                size=5
            )
            + pn.guides(
                color="none"
            )
            + pn.theme(
                figure_size=(window + 2, contig_count / 2.5),
            )
        )
        return plt

    def as_series(self) -> pd.Series:
        """Information on variant as a Pandas series.

        Provides the ID, gene, size and ratio, along with a boolean for
        all MAGs in which the variant is or is not observed in. MAGs without
        the gene will not be included.
        """

        all_mags: list[str] = list(self.single_in.union(self.duplicated_in))
        res: pd.Series = pd.Series(
            index=all_mags,
            data=[x in self.duplicated_in for x in all_mags]
        )
        for p in ['id', 'gene', 'size', 'ratio']:
            res[p] = getattr(self, p)
        return res

    @multimethod
    def as_table(mgs: str, variants: List[TandemIndel]) -> pd.DataFrame:
        """Information on variants from a single MGS as a DataFrame."""
        df: pd.DataFrame = (
            pd.DataFrame([x.as_series() for x in variants])
            .assign(mgs=mgs)
        )
        front: List[str] = ["id", "gene", "size", "ratio"]
        df = df[
            front +
            sorted([x for x in df.columns if x not in front])
        ]
        return df

    @as_table.register
    def _(variants: Dict[str, List[TandemIndel]]) -> pd.DataFrame:
        """Information on variant from multiple MGS as a DataFrame."""
        return pd.concat(
            [TandemIndel.as_table(m, v) for m, v in variants.items()]
        )

    as_table = staticmethod(as_table)


def merge_contig_dicts(
        list_dicts: List[Dict[str, Set[Contig]]]
) -> Dict[str, Set[Contig]]:
    all_keys: Set[str] = set.union(*(set(x.keys()) for x in list_dicts))
    res: Dict[str, Set[Contig]] = dict()
    for k in all_keys:
        contigs: Set[Contig] = set.union(
            *(set() if k not in x else x[k] for x in list_dicts)
        )
        res[k] = contigs
    return res