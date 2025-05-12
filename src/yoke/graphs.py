"""Analysis of gene-level pangenome graphs of MGS.

Provides functions to create and analyse pangenome graphs for MGS objects.
"""
from __future__ import annotations
from typing import List, Set, Any, NamedTuple, Tuple

import networkx
import networkx as nx

from yoke import MGS, MAG


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
        for a, b, count in self.mgs.edges:
            graph.add_edge(a, b,
                           weight=sum(count.values()),
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
                    res.append(TandemIndel(
                        gene=v,
                        vertex=v,
                        duplicated_in=dup_set,
                        single_in=l_set
                    ))
        return res

class TandemIndel(NamedTuple):
    gene: int
    vertex: int
    duplicated_in: Set[str]
    single_in: Set[str]

    @property
    def size(self) -> Tuple[int, int]:
        """Number of MAGs which have duplicate and single respectively."""

        return len(self.duplicated_in), len(self.single_in)

    @property
    def ratio(self) -> float:
        """Ratio of duplicate to single MAGs."""

        d, s = self.size
        return d / s

    def get_contigs(self, graph: MGSGraph) -> Any:
        """Find the contigs which support this variant."""

        