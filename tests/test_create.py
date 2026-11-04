from typing import List

import numpy as np
import pytest

from yoke.create import Contig


@pytest.fixture
def simple_contig() -> Contig:
    contig: Contig = Contig.from_str(
        "1,2,3,4,5,4,5,8,3,5,4,1,10,32,42", name_prefix="C"
    )[0]
    return contig


def test_contig_from_str():
    contigs: List[Contig] = Contig.from_str(
        "1,2,3,4,5,,5,6,3,2,2,1", name_prefix="C"
    )
    assert len(contigs) == 2
    assert contigs[0].genes == list(range(1, 6))
    assert contigs[1].genes == [5, 6, 3, 2, 2, 1]


def test_contig_find_indices(simple_contig):
    contig: Contig = simple_contig
    indices = contig.get_indices([4, 5])

    expected: np.ndarray = np.array([[3, 4], [5, 6], [5, 4], [10, 9]])
    assert np.array_equal(expected, indices)

    # Test with a pattern which does not appear, expect empty array
    indices = contig.get_indices([666])
    assert indices.shape == (0, 1)


def test_contig_get_context(simple_contig):
    contig: Contig = simple_contig
    ctxt = contig.get_context(genes=[4, 5], flanking=int(2))
    expect = np.array(
        [
            [2, 3, 4, 5, 4, 5],
            [4, 5, 4, 5, 8, 3],
            [3, 4, 5, 4, 5, 8],
            [8, 3, 5, 4, 1, 10],
        ]
    )
    assert np.array_equal(ctxt, expect)


def test_contig_reverse(simple_contig):
    """Contig reverse creates a new contig."""
    rev: Contig = Contig.reverse(simple_contig)
    assert rev is not simple_contig
    assert rev.genes is not simple_contig.genes


def test_contig_get_contexts(simple_contig):
    """Get the context of a pattern of genes from multiple contigs."""
    ctgs: List[Contig] = [simple_contig, Contig.reverse(simple_contig)]
    Contig.get_contexts(ctgs, [4, 5], flanking=2)
