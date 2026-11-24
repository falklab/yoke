import pathlib
from enum import Enum
from importlib.resources import files
from typing import NamedTuple


class ModuleInfo(NamedTuple):
    """Information about modules which are bundled with yoke."""

    name: str
    """Short name of the modules used in output files (e.g. BSB)"""
    about: str
    """Longer description of the modules provided for context"""
    citation: str
    """Citation for these modules"""
    definition: pathlib.Path
    """Path to the module definition (steps, enzymes, etc)"""
    description: pathlib.Path
    """Path to module descriptions (semantic meaning of module IDs)"""
    hierarchy: pathlib.Path
    """Path to module hierarchy"""


class BundledModules(Enum):
    GBM = ModuleInfo(
        name="GBM",
        about=(
            "Human gut-brain modules (MGB001-MGB056): microbial production "
            "and degradation of neuroactive compounds. Based on KO "
            "annotations."
        ),
        citation=(
            "Valles-Colomer et al. The neuroactivity potential of the human "
            "gut microbiota in quality of life and depression (2019). "
            "https://doi.org/10.1038/s41564-018-0337-x"
        ),
        definition=pathlib.Path(
            files("yoke.data").joinpath(pathlib.Path("GBM", "modGBM.list"))
        ),
        description=pathlib.Path(
            files("yoke.data").joinpath(pathlib.Path("GBM", "modGBM.descr"))
        ),
        hierarchy=pathlib.Path(
            files("yoke.data").joinpath(
                pathlib.Path("GBM", "modGBM_hiera.txt")
            )
        ),
    )
    BSB = ModuleInfo(
        name="BSB",
        about=(
            "Human gut metabolic modules (MF0001-MF0133). Based on KO "
            "annotations."
        ),
        citation=(
            "Vieira-Silva et al. Species–function relationships shape "
            "ecological properties of the human gut microbiome (2016). "
            "https://doi.org/10.1038/nmicrobiol.2016.88"
        ),
        definition=pathlib.Path(
            files("yoke.data").joinpath(pathlib.Path("BSB", "modg.list"))
        ),
        description=pathlib.Path(
            files("yoke.data").joinpath(pathlib.Path("BSB", "modg.descr"))
        ),
        hierarchy=None,
    )
