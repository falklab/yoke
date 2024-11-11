from typing import List

from create import Contig


def test_from_str():
    contigs: List[Contig] = Contig.from_str(
        "1,2,3,4,5,,5,6,3,2,2,1",
        name_prefix="C"
    )
    assert len(contigs) == 2
    assert contigs[0].genes == list(range(1,6))
    assert contigs[1].genes == [5, 6, 3, 2, 2, 1]


