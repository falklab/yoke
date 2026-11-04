"""
In MATAFILER output, MAGs are defined as a list of contig IDs, and as a list of gene IDs.
There are circumstances where we might desire a simple FASTA file containing these
sequences for analysis.
This module provides functions for generating these files either on the fly to stdout,
or in bulk writing to disk.
"""

from __future__ import annotations

import gzip
import itertools
import logging
import multiprocessing
import pathlib
from collections import defaultdict
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Optional, Dict, Iterable, List, Tuple, Any

import pyfastx
from tqdm import tqdm

from yoke.create import MAG
from yoke import create
from yoke.annotations import _csv_field_adjust

logger: logging.Logger = logging.getLogger(__name__)


@dataclass
class MAGSequences:
    """Information related to sequences and IDs for a MAG."""

    mag: MAG
    sample_assembly: Path
    contig_ids: List[str] = field(default_factory=list)

    def get_contig_link_file(self, binner: Optional[str] = None) -> Path:
        """Identify the file which links MAG name to contig ID.

        Typically an assembly directory has subdirectories `./Binning/{binner},
        and in this the file {samplename} links MAG to contig.

        Args:
            binner: Only use this results from this binner. Otherwise, the
            first directory found will be used.

        Returns:
            Path to the contig link file
        """

        binning_dir: Path = self.sample_assembly / "Binning"
        if binner is None:
            # Take first subdirectory
            binning_dir = binning_dir / next(
                d for d in binning_dir.iterdir() if d.is_dir()
            )
            logger.warning(
                "Binner not specified, using first directory: %s",
                str(binning_dir.parts[-1]),
            )
        else:
            binning_dir = binning_dir / binner

        # Contig link file is the one which matches the sample name
        contig_file = next(
            path
            for path in binning_dir.iterdir()
            if (path.is_file() and path.name in self.mag.name)
        )

        return contig_file

    @property
    def contig_fasta(self) -> Path:
        """Identify path to FASTA containing all sequences."""
        return self.sample_assembly / "scaffolds.fasta.filt"

    def __open_gzip(
        self,
        output_dir: Path,
        suffix=".contigs.fna.gz",
        **kwargs,
    ) -> gzip.GzipFile:
        """Open a gzip file for writing sequences to/

        Args:
            output_dir: Path to the output directory
            suffix: File suffix to use
            **kwargs: Passed to gzip.open

        Returns:
            Connection to Gzip file
        """
        opath = output_dir / f"{self.mag.name}{suffix}"
        return gzip.open(opath, **kwargs)

    @staticmethod
    def from_assembly_dir(mag: MAG, assembly_dir: Path) -> MAGSequences:
        """Create a MAGSequences object for a MAG looking in a project
        assembly directory."""

        return MAGSequences(
            mag=mag,
            sample_assembly=SequenceExtractor(assembly_dir).get_mag_assembly(
                mag
            ),
        )

    @staticmethod
    def group_mags(
        mags: Iterable[MAGSequences],
    ) -> Dict[Path, List[MAGSequences]]:
        """Group MAGs by assembly of origin.

        Uses MAG name to resolve which assembly the MAG originates from, and group these by

        Args:
            mags: A collection of MAGs from the same project.

        Returns:
           Dict[Path, List[MAG]]]: A dict in which the key is an assembly directory, and a
            value a list of all the MAGs in the input associated with that assembly.
        """

        collected = defaultdict(list)
        for mag in mags:
            collected[mag.sample_assembly].append(mag)

        return collected

    @staticmethod
    def add_contig_ids(
        mags: List[MAGSequences], **kwargs
    ) -> List[MAGSequences]:
        """Add contig IDs to a group of MAGs from the same origin assembly. This is
        more efficient than handling each MAG individually.

        Args:
              mags: A collection of MAGs from the same project.
              **kwargs: Passed to MAGSequences.get_contig_link_file

        Returns:
              Input mags with contig IDs added.
        """

        logger.debug("Attach contig IDs to %s MAGs", len(mags))
        logger.debug("Validate all from sample sample assembly")
        sample_paths = set(mag.sample_assembly for mag in mags)
        if len(sample_paths) != 1:
            logger.fatal(
                "All MAGs must originate from same assembly when extracting IDs. "
                "Found %s for MAGs %s",
                sample_paths,
                [mag.mag.name for mag in mags],
            )
            raise Exception(
                "MAGs passed to get_contig_ids must have the same "
                "assembly path"
            )

        ids = {mag.mag.name.split("__")[1]: mag for mag in mags}
        with open(mags[0].get_contig_link_file(**kwargs), "rt") as f:
            for line in f:
                parts = line.rstrip().split("\t")
                # Do not append if not a MAG we care about
                if parts[1] in ids.keys():
                    ids[parts[1]].contig_ids.append(parts[0])

        return mags

    @staticmethod
    def extract_contigs(
        mags: List[MAGSequences],
        output_dir: Path,
        binner: Optional[str] = None,
        processes: int = 1,
        progress_bar: bool = False,
        **kwargs,
    ) -> List[Path]:
        """Extract contigs from MAGs to compressed FASTA files

        Write individual FASTA files representing the MAGs to a directory.

        Args:
            mags: List of MAGs to write. Can originate from different assemblies.
            output_dir: Path to write output in
            processes: Number of multiprocessing processes to use
            **kwargs: Passed to multiprocessing.Pool construtor

        Return:
            Paths for all the created files
        """

        grouped_mags = MAGSequences.group_mags(mags)
        created = []
        if processes > 1:
            with multiprocessing.Pool(processes=processes, **kwargs) as pool:
                val_iter = (
                    (x, output_dir, binner) for x in grouped_mags.values()
                )
                if progress_bar:
                    val_iter = tqdm(val_iter, total=len(grouped_mags))
                c = pool.starmap(_multiproc_run, val_iter, **kwargs)
            created = list(itertools.chain.from_iterable(c))
        else:
            val_iter = grouped_mags.items()
            if progress_bar:
                val_iter = tqdm(val_iter, total=len(grouped_mags))
            for ass_dir, mags in val_iter:
                mags = MAGSequences.add_contig_ids(mags, binner=binner)
                c = MAGSequences._extract_contigs_sample(
                    mags, output_dir=output_dir
                )
                created = [*created, *c]
        return created

    @staticmethod
    def _extract_contigs_sample(
        mags: List[MAGSequences], output_dir: Path, **kwargs
    ) -> List[Path]:
        """Extract contigs from MAGs from a single assembly to compressed FASTA files.

        Args:
            mags: List of MAGs to write, must originate from a single assembly
            output_dir: Path to write output in

        Return:
            Paths for all the created files
        """

        # Create a dict mapping contig id to MAG
        # Probably a more elegant way to do this
        contig_map: Dict = dict()
        for mag in mags:
            for contig in mag.contig_ids:
                contig_map[contig] = mag

        with ExitStack() as stack:
            # Will have a dynamic number of file contexts to manage
            output_files: Dict = {
                m.mag.name: stack.enter_context(
                    m.__open_gzip(output_dir, mode="wt")
                )
                for m in mags
            }

            contig_fasta = pyfastx.Fasta(
                str(mags[0].contig_fasta), build_index=False
            )

            for name, seq in contig_fasta:
                try:
                    mag_name = contig_map[name].mag.name
                    output_files[mag_name].write(f">{name}\n")
                    output_files[mag_name].write(f"{seq}\n")
                except KeyError as e:
                    logger.debug("%s not found in contig map", name)

        return [Path(o.name) for o in output_files.values()]


class SequenceExtractor:
    """Extract sequences for MAGs from their respective assemblies."""

    SAMPLE_SEPARATOR: str = "__"
    """String which delimits sample name at start of MAG name."""
    ASSEMBLY_PATH: Path = Path("assemblies") / "metag"
    """Location of the assembly within a sample or assembly group directory."""
    CONTIG_FASTA: Path = Path("scaffolds.fasta.filt")
    """Location of the contig FASTA within a sample assembly directory."""
    PROTEIN_FASTA: Path = Path("genePred") / "proteins.shrtHD.faa.gz"
    """Location of the gene protein FASTA within a sample assembly directory."""
    TRIM_NAME_REGEX: re.Pattern = re.compile(r"^(.*)M\d{1,2}$")

    def __init__(self, assembly_path: Path) -> None:
        self.assembly_path: Path = assembly_path

    def get_mag_assembly(self, mag: MAG) -> Path:
        """Find the path for the assembly of a MAG.

        Args:
            mag(MAG): the MAG for find the assembly for

        Returns:
            Path: Path to the assembly directory for this MAG
        """

        sample_name = mag.name.split(self.SAMPLE_SEPARATOR)[0]
        sample_dir = self.assembly_path / sample_name
        if not sample_dir.is_dir():
            # Try dropping the trailing M\d{1,2}
            logger.debug("%s is not a directory, will try variant", sample_dir)
            matches: Optional[re.Match] = self.TRIM_NAME_REGEX.search(
                sample_name
            )
            if matches is not None:
                sample_name = str(matches.group(1))
                logger.debug("Trimmed sample name to %s", sample_name)
                sample_dir = self.assembly_path / sample_name
            else:
                logger.debug("Unable to trim")
            if not sample_dir.is_dir():
                logger.fatal(
                    "Directory for MAG %s not found (%s)", mag.name, sample_dir
                )
                raise NotADirectoryError(
                    f"{sample_dir} not found or is is not a directory"
                )

        # If a sample is part of a coassembly, it will have an `assembly.txt` file with the
        # location of the assembly
        assembly_txt = sample_dir / self.ASSEMBLY_PATH / "assembly.txt"
        assembly_path = None
        if assembly_txt.is_file():
            assembly_path = Path(
                assembly_txt.read_text().splitlines()[0].strip()
            )
            logger.debug(f"Assembly path from txt: {assembly_path}")
            if not assembly_path.is_dir():
                logger.warning(
                    (
                        "MAG %s: assembly.txt found, but assembly "
                        "directory does not exist (%s, %s)"
                    ),
                    mag.name,
                    str(assembly_txt),
                    str(assembly_path),
                )
                assembly_path = None

        if assembly_path is None:
            logger.debug(
                (
                    "MAG %s: assembly.txt not found (%s) or referred to directory did"
                    "not exist"
                ),
                mag.name,
                assembly_txt,
            )
            # Assume this is a single sample assembly, check for existence of scaffolds
            assembly_path = sample_dir / self.ASSEMBLY_PATH

        if not (assembly_path / self.CONTIG_FASTA).is_file():
            raise FileNotFoundError(
                "MAG %s: No scaffold FASTA found in directory %s",
                mag.name,
                str(assembly_path),
            )

        return assembly_path.absolute()


def _multiproc_run(
    mags: List[MAGSequences], output_dir: Path, binner: Optional[str]
):
    mags = MAGSequences.add_contig_ids(mags, binner=binner)
    return MAGSequences._extract_contigs_sample(mags, output_dir)


def cli_extract_contigs(
    magvsgc: pathlib.Path,
    assembly: pathlib.Path,
    output: pathlib.Path,
    completeness: float = 0,
    contamination: float = 100,
    binner: Optional[str] = None,
    processes: int = 1,
) -> None:
    """Count how many genes have each functional annotation in all MAGS."""

    _csv_field_adjust()
    logger.setLevel(logging.INFO)

    logger.info("extract_contigs")
    logger.info("|- Parameters -|")
    logger.info("MAGvsGC file: %s", magvsgc)
    logger.info("Assembly Directory: %s", assembly)
    logger.info("Bin Directory: %s", binner)
    logger.info("Completeness: %s", completeness)
    logger.info("Contamination: %s", contamination)
    logger.info("Processes: %s", processes)
    logger.info("Output: %s", output)

    pth_out: pathlib.Path = output
    pth_out.mkdir(parents=False, exist_ok=True)

    logger.info("Reading MAGs from MAGvsGC")
    mgs_gen = [
        x
        for x in create.generate_mgs_from_magvsgc(
            magvsgc,
            min_completeness=completeness,
            max_contamination=contamination,
            compressed=True,
        )
    ]

    logger.info("Found %s MGS", len(mgs_gen))
    logger.debug("Flattening MGS to MAGs")
    mags = list(itertools.chain.from_iterable(x.mags for x in mgs_gen))
    logger.info("Found %s MAGs", len(mags))
    mag_seqs = [
        MAGSequences.from_assembly_dir(x, assembly_dir=assembly) for x in mags
    ]
    logger.info("Extracting contigs...")
    logger.info(
        "Note: Progress bar updated when extraction begins for each assembly"
    )
    files = MAGSequences.extract_contigs(
        mags=mag_seqs,
        output_dir=pth_out,
        binner=binner,
        processes=processes,
        progress_bar=True,
    )

    logger.info("Completed extracting contigs to %s", pth_out)
