import pathlib
from enum import Enum
from importlib.resources import files
from typing import NamedTuple


class ModuleInfo(NamedTuple):
    """Information about modules which are bundled with yoke"""

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
        about="",
        citation="",
        definition=pathlib.Path(
            files("yoke.data")
            .joinpath(pathlib.Path("GBM", "modGBM.list"))
        ),
        description=pathlib.Path(
            files("yoke.data")
            .joinpath(pathlib.Path("GBM", "modGBM.descr"))
        ),
        hierarchy=pathlib.Path(
            files("yoke.data")
            .joinpath(pathlib.Path("GBM", "modGBM_hiera.txt"))
        )
    )
