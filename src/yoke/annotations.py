"""Link between gene IDs and functional annotations.

Genes may be annotated against several functional schema, such as KEGG
Orthologs, eggNOGs, GBM modules etc. This module provides methods to handle
the connection between gene identities and these annotation schema.
"""

from __future__ import annotations

import csv
import dataclasses
import gzip
import itertools
import logging
import pathlib
import subprocess
import sys
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import Enum
from typing import (
    Any,
    Literal,
)

import numpy as np
import pandas as pd
import tqdm
from multimethod import multimethod
from scipy.stats import barnard_exact, boschloo_exact, fisher_exact
from statsmodels.stats import multitest

from yoke import create
from yoke.create import MAG, MGS
from yoke.data.utils import BundledModules, ModuleInfo
from yoke.utils import GzipManager

logger: logging.Logger = logging.getLogger(__name__)


@dataclass
class GeneFunction:
    """One function, and the gene IDs which have been annotated with it.

    Attributes:
        id: Function ID.
        desc: Description of function.
        color: Color to use when representing this gene.
    """

    id: str
    desc: str | None = None
    color: str | None = None
    gene_ids: set[int] = dataclasses.field(default_factory=set)

    @property
    def sparse_series(self) -> pd.Series:
        a = pd.Series(
            [1] * len(self.gene_ids),
            dtype=pd.SparseDtype(np.dtype("int_"), fill_value=0),
        )
        return a


class FunctionSchema(Enum):
    """Functional annotation schema supported by default load methods."""

    KEGG = 1


class HierarchicalSchema:
    """Functional annotation schema with one or more levels of hierarchy.

    Hierarchical annotation schema can be represented by a tree, allowing
    for multifurcation. Examples include GBM modules or BRITE Hierarchies.
    """

    def __init__(
        self,
        name: str,
        hierarchy: pd.DataFrame,
        gene_map: pd.DataFrame,
        **kwargs,
    ):
        """Initialise a hierarchical annotation schema.

        Args:
            name: Name of the annotation schema.
            hierarchy: Table where the index is function IDs, and each
                column represents a level of the hierarchy.
            gene_map: Table with the index being gene ID, and the only column
                being function ID. Values should not be NA.
            **kwargs: Currently unused.
        """
        self.__name: str = name
        self.__hierarchy: pd.DataFrame = hierarchy
        # Drop NA
        gm = gene_map.loc[~gene_map.iloc[:, 0].isna(), :]
        self.__gene_map: dict = dict(zip(gm.index, gm.iloc[:, 0], strict=True))
        self.__functions: dict[str, GeneFunction] = self.__make_gene_functions(
            hierarchy, gene_map
        )

    @property
    def name(self) -> str:
        return self.__name

    @property
    def hierarchy(self) -> pd.DataFrame:
        return self.__hierarchy

    @property
    def functions(self) -> dict[str, GeneFunction]:
        return self.__functions

    @property
    def gene_map(self) -> pd.DataFrame:
        return self.__gene_map

    @property
    def gene_matrix(self) -> pd.DataFrame:
        mat: pd.DataFrame = pd.concat(
            [f.sparse_series for f in self.functions.values()], axis=1
        )
        return mat

    def __make_hierarchy(
        self, hierarchy: pd.DataFrame | None, annotation: pd.DataFrame
    ) -> pd.DataFrame:
        """Format a hierarchy table, or make a dummy if there is no table."""
        if hierarchy is not None:
            return hierarchy
        # TODO: Unfinished; make a single level hierarchy with all functions
        # from the annotation
        funcs: pd.Series = annotation.iloc[:, 0].unique()  # noqa: F841

    def get_gene_function(self, gene_id) -> GeneFunction | None:
        """Get the annotation for a gene, or None if it has no annotation."""
        if gene_id not in self.gene_map:
            return None
        function_id = self.gene_map[gene_id]
        return self.functions[function_id]

    @multimethod
    def function_counts(self, mag: MAG):
        """Count genes in this MAG with each function in this schema.

        Args:
            mag: MAG to produce counts for.
        """
        f_counts: Counter = Counter(
            gf.id
            for gf in (self.get_gene_function(gene) for gene in mag.genes)
            if gf is not None
        )
        counts: pd.Series = pd.Series(f_counts, name=mag.name).astype("int")
        return counts

    @multimethod
    def function_counts(self, mgs: MGS):  # noqa: F811 (overload)
        """Count genes in each MAG of an MGS with each function in this schema.

        Args:
            mgs: MGS to produce counts for.
        """
        # Make a multiindex with MAG and MGS - MGS is redundant, but it
        # makes columnwise concating easier
        cols: pd.MultiIndex = pd.MultiIndex.from_tuples(
            zip([mag.name for mag in mgs.mags], itertools.cycle([mgs.name])),
            names=["mag", "mgs"],
        )
        df = pd.concat(
            [
                self.function_counts(mag)
                for mag in tqdm.tqdm(mgs.mags, desc=mgs.name)
            ],
            axis=1,
        )
        df.columns = cols
        return df

    @multimethod
    def function_counts(  # noqa: F811 (overload)
        self, mgs: list[MGS]
    ) -> pd.DataFrame:
        """Count genes in each MAG of these MGS with each function in schema.

        Args:
            mgs: Collection of MGS to produce counts for.
        """
        return pd.concat([self.function_counts(m) for m in mgs], axis=1)

    @staticmethod
    def __make_gene_functions(
        hierarchy: pd.DataFrame, gene_map: pd.DataFrame
    ) -> dict[str, GeneFunction]:
        """Create GeneFunction objects for this hierarchy."""
        # Refactor for speed
        # Try group by apply
        ret_dict: dict[str, GeneFunction] = dict()
        for function_id, ann_group in gene_map.groupby(by=1):
            ret_dict[function_id] = GeneFunction(
                id=function_id,
                desc=(
                    None
                    if function_id not in hierarchy.index
                    else hierarchy.loc[function_id].iloc[-1]
                ),
                gene_ids=set(ann_group.index),
            )
        return ret_dict

    @staticmethod
    def from_file(
        name: str,
        hierarchy: str | pathlib.Path,
        annotation: str | pathlib.Path,
        delim: str = "\t",
    ) -> HierarchicalSchema:
        """Read schema from files.

        Create an annotation schema from a definition of the hierarchy,
        and an annotation file which connects genes and function IDS.
        Annotation files should be only two columns, with no headers.
        """
        h_df: pd.DataFrame = pd.read_csv(hierarchy, sep=delim, index_col=0)
        a_df: pd.DataFrame = pd.read_csv(annotation, sep=delim, index_col=0)

        return HierarchicalSchema(name=name, hierarchy=h_df, gene_map=a_df)


class RTKModules:
    """Functional modules for which module presence can be estimated using RTK.

    Groups of genes working together can be defined as functional modules, and
    we can assess whether the whole module is present and it's coverage etc in
    a set of genes. This class handles modules defined in Omixer format,
    and uses external binary `rtk` to estimate their completion.

    As estimation is handled by `rtk`, the modules definitions are never
    loaded or handled in python, only the paths to the components are needed.
    The package also bundles several module definitions, which can be
    loaded from the `RTKModules.bundled()` factory method.
    """

    RTK_DEFAULTS = dict(
        redundancy="5",
        writeExtraModEstimates=True,
        moduleCompl="0.5",
        enzymeCompl="0.5",
        collapseDblModules=True,
    )

    def __init__(
        self,
        definition: str | pathlib.Path,
        description: str | pathlib.Path,
        hierarchy: str | pathlib.Path | None,
        name: str,
    ) -> None:
        self.__definition: pathlib.Path = pathlib.Path(definition)
        self.__description: pathlib.Path = pathlib.Path(description)
        self.__hierarchy: pathlib.Path | None = (
            pathlib.Path(hierarchy) if hierarchy is not None else None
        )
        self.__name: str = name

    @property
    def definition(self) -> pathlib.Path:
        """Path to module definitions (which functions make up a module)."""
        return self.__definition

    @property
    def description(self) -> pathlib.Path:
        """Path to description of modules."""
        return self.__description

    @property
    def hierarchy(self) -> pathlib.Path:
        """Path to hierarchy of modules."""
        return self.__hierarchy

    @property
    def name(self) -> str:
        """Name of the modules, used in output formatting."""
        return self.__name

    @staticmethod
    def from_bundled(name: str) -> RTKModules:
        """Create a module scorer for one of the bundled modules."""
        if name not in BundledModules._member_names_:
            raise ValueError(
                f"Invalid module name: {name}. Valid names are: "
                f"{BundledModules._member_names_}"
            )
        mod: ModuleInfo = BundledModules[name].value
        return RTKModules(
            definition=mod.definition,
            description=mod.description,
            hierarchy=mod.hierarchy,
            name=mod.name,
        )

    @staticmethod
    def _format_rtk_args(args: dict[str, Any]) -> str:
        """Convert dictionary into shell argument list."""
        argline: str = ""
        for key, val in args.items():
            if isinstance(val, bool) and val:
                argline += f" -{key}"
                continue
            argline += f" -{key} {str(val)}"
        return argline

    @staticmethod
    def load_result(prefix: pathlib.Path) -> dict[str, pd.DataFrame]:
        """Load results from an RTK module call as tables.

        Load all resulting tables from a call to `rtk module`. The returned
        dictionary will have a key for each suffix (i.e. result['mat'],
        result['KOused'])

        Args:
            prefix: The prefix for the files (i.e. /home/a/MGS.0.GBM)

        Returns:
            Dictionary will have a key for each suffix (i.e. result['mat'],
            result['KOused'])
        """
        logger.debug("Load RTK results for %s", prefix)
        files: list[pathlib.Path] = list(
            prefix.parent.glob(f"{prefix.name}.*")
        )
        logger.debug("Found %s files: %s", len(files), files)
        res = {
            str(pth).replace(str(prefix) + ".", ""): pd.read_csv(pth, sep="\t")
            for pth in files
        }
        return res

    def score_table(
        self,
        pth_mgs_tbl: pathlib.Path,
        output_prefix: pathlib.Path | None = None,
        load_result: bool = False,
        **kwargs,
    ) -> dict[str, pd.DataFrame] | None:
        """Score the presence of modules for each item in a table.

        The input table should have functions on rows, and items (samples,
        MAGs) on columns, with entries being either coverage or counts.
        `rtk module` is used to perform scoring.

        Args:
            pth_mgs_tbl: Path to table with count of functions to score
            module presence in.
            output_prefix: Path to use as prefix for each file. Will produce
            multiple output files (i.e. prefix /home/a/GBM.MGS.0 will make
            /home/a/GBM.MGS.0.mat, /home/a/GBM.0.KOused, ...). If none,
            pth_mgs_tbl will be used.
            load_result: Load each resulting table and return in a dictionary.
            **kwargs: Extra arguments to be passed to rtk. Any boolean True
            values will be parsed to a flag. No validation is done on these
            arguments.

        Returns:
            List of paths for output files
        """
        output_prefix = pth_mgs_tbl if output_prefix is None else output_prefix
        rtk_params: dict[str, Any] = self.RTK_DEFAULTS | kwargs
        rtk_cmd: str = (
            "rtk module "
            f"-i {pth_mgs_tbl} "
            f"-o {output_prefix} "
            f"-refMods {self.definition} "
            f"-description {self.description} "
        )
        if self.hierarchy is not None:
            rtk_cmd += f"-hierarchy {self.hierarchy} "

        # Convert rtk params to command options
        logger.debug("Format extra RTK arguments: %s", rtk_params)
        extra_args: str = self._format_rtk_args(rtk_params)
        rtk_cmd = rtk_cmd + extra_args

        expect_output: pathlib.Path = output_prefix.with_name(
            f"{output_prefix.name}.mat"
        )

        logger.debug("Run RTK command: %s", rtk_cmd)
        subp: subprocess.CompletedProcess = subprocess.run(
            rtk_cmd, shell=True, capture_output=True
        )

        # Determine success
        if subp.returncode != 0 or not expect_output.exists():
            logger.error("RTK module scoring for %s failed", pth_mgs_tbl)
            logger.error("Stderr")
            logger.error(subp.stderr.decode("utf-8"))
            logger.error("Stdout")
            logger.error(subp.stdout.decode("utf-8"))
            raise Exception(subp)

        if load_result:
            return self.load_result(output_prefix)
        return None


def _combine_mgs(mag_counts: pathlib.Path, mgs: list[str]) -> pd.DataFrame:
    """Combine function counts of multiple MGS into one table.

    The returned table will have one column for each MAG in any of the listed
    MGS. This can be used to get a table of all MGS with a certain phenotype,
    or from a given genus etc.
    """
    # Check whether files appear to be gzipped
    suffix: str = (
        "tsv" if (mag_counts / f"{mgs[0]}.tsv").exists() else "tsv.gz"
    )
    mgs_dfs: list[pd.DataFrame] = [
        pd.read_csv(mag_counts / f"{x}.{suffix}", sep="\t", index_col=0)
        for x in mgs
    ]
    comb_df: pd.DataFrame = pd.concat(mgs_dfs, axis=1)
    return comb_df


def _csv_field_adjust():
    """Increase max csv field size to parse larger entries in MAGvsGC."""
    logger.debug("Increasing max field size to %s", str(sys.maxsize))
    logger.debug(
        "This has been needed in the past to parse large entries in MAGvsGC"
    )
    csv.field_size_limit(sys.maxsize)


_ENRICHMENT_TESTS = {
    "fisher": fisher_exact,
    "barnard": barnard_exact,
    "boschloo": boschloo_exact,
}


def _check_unique(index: pd.Index, what: str) -> None:
    """Raise a ValueError naming duplicated entries in an index."""
    dups = index[index.duplicated()].unique()
    if len(dups) > 0:
        raise ValueError(
            f"Duplicate {what}: {', '.join(str(x) for x in dups[:5])}"
            + (f" and {len(dups) - 5} more" if len(dups) > 5 else "")
        )


def function_enrichment(
    function_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    tests: Iterable[Literal["fisher", "barnard", "boschloo"]] = (
        "fisher",
        "barnard",
    ),
    adj_method: str = "fdr_bh",
) -> pd.DataFrame:
    """Test functions for association to two-class categorical metadata.

    Args:
        function_counts: Table with rows being a functional annoation, columns
        an MGS or MAG, and entries a numerical value. Values greater than 0
        are treated as the function being present.
        metadata: MGS or MAG on rows, and each column should be a metadata
        category to test against, having only 2 levels. NAs will be removed
        from metadata.
        tests: List of the types of test to carry out. Available are
        * fisher - fishers exact test
        * barnard - barnards exact test
        * boschloo - boschloos exact test
        all implemented by scipy
        adj_method: Multiple test correction method to apply, accepting any
        method from statsmodels.stats.multitest

    Returns:
        Table indexed by metadata column and function, with columns
        {test}_stat, {test}_p and {test}_p_adj for each requested test.
        Correction is applied separately for each metadata column and test.
        Metadata columns without exactly 2 levels are skipped with a warning.

    Raises:
        ValueError: Unknown test or correction method, duplicated
        identifiers, non-numeric counts, or no MGS/MAGs shared between
        function_counts and metadata.
    """
    tests = list(tests)
    unknown = set(tests) - set(_ENRICHMENT_TESTS)
    if unknown:
        raise ValueError(f"Unknown enrichment tests: {sorted(unknown)}")
    # Check correction method before running tests, which can be slow
    try:
        multitest.multipletests([0.5], method=adj_method)
    except ValueError:
        raise ValueError(
            f"Unknown multiple test correction method: {adj_method}"
        ) from None

    _check_unique(function_counts.index, "function IDs in function_counts")
    _check_unique(function_counts.columns, "MGS/MAGs in function_counts")
    _check_unique(metadata.index, "MGS/MAGs in metadata")
    non_numeric = [
        str(c)
        for c, d in function_counts.dtypes.items()
        if not pd.api.types.is_numeric_dtype(d)
    ]
    if non_numeric:
        raise ValueError(
            f"function_counts has non-numeric columns: {non_numeric[:5]}"
        )

    # Inform how many MGS/MAGs are not matched in metadata
    dropped_tax = set(function_counts.columns) - set(metadata.index)
    logger.info("%s MGS/MAGs not matched in metadata", len(dropped_tax))
    logger.debug("%s", str(dropped_tax))
    if len(dropped_tax) == len(function_counts.columns):
        raise ValueError(
            "No MGS/MAGs in function_counts match the metadata index. Check "
            "the identifiers are the same, including their type (e.g. int "
            "vs str)."
        )

    if len(function_counts) == 0:
        logger.warning("function_counts has no functions, nothing to test")
        return pd.DataFrame()

    present = function_counts > 0

    res = []
    for col in metadata.columns:
        md_col = metadata[col].dropna()
        shared = present.columns.intersection(md_col.index)
        if len(shared) == 0:
            logger.warning(
                "Metadata column %s has no non-NA values for MGS/MAGs in "
                "function_counts, skipping",
                col,
            )
            continue
        md_col = md_col.loc[shared]
        try:
            levels = sorted(md_col.unique())
        except TypeError:
            # Mixed types can't be compared, so sort on string form
            levels = sorted(md_col.unique(), key=str)
        # Put True first so odds ratios are True vs False
        if set(levels) == {True, False}:
            levels = [True, False]
        if len(levels) != 2:
            logger.warning(
                "Metadata column %s has %s levels (%s) for MGS/MAGs in "
                "function_counts after removing NAs, needs 2, skipping",
                col,
                len(levels),
                ", ".join(str(x) for x in levels[:5]),
            )
            continue
        logger.info("Testing %s with %s MGS/MAGs", col, len(shared))

        # Each 2x2 table is built from counts of presence within each
        # metadata level, so all levels are always represented
        grp = md_col == levels[0]
        fn = present[shared]
        in_a = fn.loc[:, grp.values].sum(axis=1)
        in_b = fn.loc[:, (~grp).values].sum(axis=1)
        n_a, n_b = grp.sum(), (~grp).sum()

        col_res = pd.DataFrame(index=fn.index)
        for test in tests:
            test_fn = _ENRICHMENT_TESTS[test]
            stats = [
                test_fn(np.array([[a, n_a - a], [b, n_b - b]]))
                for a, b in zip(in_a, in_b, strict=True)
            ]
            col_res[f"{test}_stat"] = [r.statistic for r in stats]
            col_res[f"{test}_p"] = [r.pvalue for r in stats]
            col_res[f"{test}_p_adj"] = multitest.multipletests(
                col_res[f"{test}_p"].fillna(1.0), method=adj_method
            )[1]
        col_res.index.name = "function"
        res.append(col_res.assign(metadata=col))

    if len(res) == 0:
        return pd.DataFrame()
    return pd.concat(res).reset_index().set_index(["metadata", "function"])


def _read_id_table(pth: pathlib.Path) -> pd.DataFrame:
    """Read a TSV indexed by the first column, keeping IDs as strings.

    Column headers are always read as strings, so identifiers on the index
    must be too for them to match, even when they look numeric (including
    leading zeros). Passing converters with index_col=0 is not used, as it
    is ignored when the first column has no header.
    """
    tbl: pd.DataFrame = pd.read_csv(pth, sep="\t", converters={0: str})
    tbl = tbl.set_index(tbl.columns[0])
    if str(tbl.index.name).startswith("Unnamed:"):
        tbl.index.name = None
    return tbl


def cli_function_enrichment(
    counts: list[pathlib.Path],
    metadata: pathlib.Path,
    output: pathlib.Path,
    tests: Iterable[str] = ("fisher", "barnard"),
    adj_method: str = "fdr_bh",
    columns: list[str] | None = None,
) -> None:
    """Test functions for association with two-level metadata categories.

    Args:
        counts: Paths to tab separated tables of function counts, with
        functions on rows and MGS or MAGs on columns. Multiple tables are
        combined column-wise, with missing functions treated as absent.
        metadata: Path to tab separated metadata table, with MGS or MAG on
        rows matching the count table columns.
        output: Path to write the results table to.
        tests: Tests to carry out, see function_enrichment.
        adj_method: Multiple test correction method, see function_enrichment.
        columns: Metadata columns to test. All columns tested if not given.
    """
    _csv_field_adjust()
    logger.setLevel(logging.INFO)

    logger.info("function_enrichment")
    logger.info("|- Parameters -|")
    logger.info("Function counts: %s", ", ".join(str(x) for x in counts))
    logger.info("Metadata: %s", metadata)
    logger.info("Output: %s", output)
    logger.info("Tests: %s", ", ".join(tests))
    logger.info("Adjustment method: %s", adj_method)
    logger.info("Columns: %s", columns)

    logger.info("Reading function counts")
    tables: list[pd.DataFrame] = []
    for pth in counts:
        tbl: pd.DataFrame = _read_id_table(pth)
        _check_unique(tbl.index, f"function IDs in {pth}")
        tables.append(tbl)
    _check_unique(
        pd.Index([c for t in tables for c in t.columns]),
        "MGS/MAGs across function count tables",
    )
    function_counts: pd.DataFrame = pd.concat(tables, axis=1).fillna(0)

    logger.info("Reading metadata")
    md: pd.DataFrame = _read_id_table(metadata)
    if columns is not None:
        missing: set[str] = set(columns) - set(md.columns)
        if missing:
            raise ValueError(f"Metadata columns not found: {sorted(missing)}")
        md = md[columns]

    res: pd.DataFrame = function_enrichment(
        function_counts, md, tests=tests, adj_method=adj_method
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Writing results to %s", output)
    res.to_csv(output, sep="\t")


def cli_function_count(
    magvsgc: pathlib.Path,
    annotations: pathlib.Path,
    descriptions: pathlib.Path | None,
    output: pathlib.Path,
    name: str,
    completeness: float = 0,
    contamination: float = 100,
    compress: bool = True,
) -> None:
    """Count how many genes have each functional annotation in all MAGS."""
    _csv_field_adjust()
    logger.setLevel(logging.INFO)

    logger.info("function_count")
    logger.info("|- Parameters -|")
    logger.info("MAGvsGC file: %s", magvsgc)
    logger.info("Gene annotations: %s", annotations)
    logger.info("Annotation descriptions: %s", descriptions)
    logger.info("Output: %s", output)
    logger.info("Name: %s", name)
    logger.info("Completeness: %s", completeness)
    logger.info("Contamination: %s", contamination)
    logger.info("Compress: %s", compress)

    pth_out: pathlib.Path = output / name
    pth_out.mkdir(parents=False, exist_ok=True)

    logger.info("Reading annotation descriptions")
    description: pd.DataFrame = pd.read_csv(
        descriptions, sep="\t", header=None, index_col=0
    )
    logger.info("Reading gene annotations")
    annotation: pd.DataFrame = pd.read_csv(
        annotations, sep="\t", header=None, index_col=0
    )
    logger.info("Making annotation object")
    s: HierarchicalSchema = HierarchicalSchema(
        name=name, hierarchy=description, gene_map=annotation
    )

    logger.info("Creating function counts for each MGS")
    mgs_gen = create.generate_mgs_from_magvsgc(
        magvsgc,
        min_completeness=completeness,
        max_contamination=contamination,
        compressed=True,
    )

    # Collect means to make a summary DF also
    means: list[pd.Series] = []
    comp_args: dict = dict(compression="gzip") if compress else dict()
    for mgs in mgs_gen:
        pth_mgs: pathlib.Path = (
            pth_out / f"{mgs.name}.tsv.gz"
            if compress
            else pth_out / f"{mgs.name}.tsv"
        )
        pth_mgs_pa: pathlib.Path = (
            pth_out / f"presence.{mgs.name}.tsv.gz"
            if compress
            else pth_out / f"presence.{mgs.name}.tsv"
        )
        tbl: pd.DataFrame = s.function_counts(mgs).fillna(0.0)
        tbl.columns = tbl.columns.droplevel("mgs")
        logger.debug("Write to %s", pth_mgs)
        tbl.to_csv(pth_mgs, sep="\t", **comp_args)
        tbl_presence: pd.DataFrame = tbl > 0.0
        logger.debug("Write to %s", pth_mgs_pa)
        tbl_presence.to_csv(pth_mgs_pa, sep="\t", **comp_args)
        mean: pd.Series = tbl.mean(axis=1)
        # The series name becomes the MGS column in the combined table
        mean.name = mgs.name
        means.append(mean)
    pth_all: pathlib.Path = output / f"{name}.tsv"
    logger.debug("Producing per MGS mean table")
    df_means: pd.DataFrame = pd.concat(means, axis=1)
    logger.info("MGS means saved to %s", pth_all)
    df_means.to_csv(pth_all, sep="\t")


def cli_function_hoist(
    mag_counts: pathlib.Path,
    magvsgc_file: pathlib.Path,
    output: pathlib.Path,
    ranks: Sequence[str] = (
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
    ),
    compressed: bool = True,
    cols: dict[str, str] | None = None,
) -> None:
    """Summarise function counts per MAG to different taxonomic levels."""
    _csv_field_adjust()
    logger.setLevel(logging.INFO)

    logger.info("function_summarise")
    logger.info("|- Parameters -|")
    logger.info("MAG function counts: %s", mag_counts)
    logger.info("MAGvsGC file: %s", magvsgc_file)
    logger.info("Ranks: %s", ",".join(ranks))
    logger.info("Output: %s", output)
    logger.info("Compress: %s", compressed)

    # Take per MGS taxonomy from MAGvsGC file
    if cols is None:
        cols = create.DEF_MAGCVSGC_COLS

    logger.info("Reading MAGvsGC file %s", magvsgc_file)
    openfn = (lambda x: gzip.open(x, "rt")) if compressed else open
    with openfn(magvsgc_file) as magvsgc:
        reader = csv.DictReader(magvsgc, delimiter="\t")
        mgs_tax: dict[str, pd.Series] = dict()
        for row in reader:
            mgs_id: str = row[cols["mgs"]]
            if (mgs_id not in mgs_tax) and (len(mgs_id) > 0):
                tdict = {r: row[r] for r in ranks}
                tdict["mgs"] = mgs_id
                mgs_tax[mgs_id] = tdict
    tax_df: pd.DataFrame = pd.DataFrame(
        [pd.Series(x) for x in mgs_tax.values()]
    )

    # Group MGS by each level of each requested rank
    for rank in ranks:
        logger.info("Summarise to rank %s", rank)
        pth_rank = output / rank
        pth_rank.mkdir(parents=True, exist_ok=True)
        means: list[pd.Series] = []
        tax: list[str] = list(tax_df[rank].unique())
        for val in (pbar := tqdm.tqdm(tax)):
            pbar.set_description(f"{rank}: {val}")
            pth_full = pth_rank / f"{val}.tsv"
            df = _combine_mgs(
                mag_counts, tax_df.loc[tax_df[rank] == val, "mgs"].tolist()
            )
            logger.debug("Write to %s", pth_full)
            df.to_csv(pth_full, sep="\t")
            logger.debug("Calculate mean for taxon")
            mean_series = df.mean(axis=1)
            mean_series.name = val
            means.append(mean_series)
        logger.info("Writing taxon mean counts")
        df_mean = pd.DataFrame(means)
        df_mean.T.to_csv(output / f"{rank}.tsv", sep="\t")


def cli_function_modules(
    tables: Iterable[pathlib.Path],
    name: str,
    output_dir: pathlib.Path,
    modules: str | None = None,
    definition: pathlib.Path | None = None,
    description: pathlib.Path | None = None,
    hierarchy: pathlib.Path | None = None,
) -> None:
    """Score module presence based on a table of functional annotations.

    After counting the presence of functional genes in MAGs, this function
    can estimate the presence of pathways (defined roughly as sequence of
    genes) in the MAGs. This is done using external binary RTK, and modules
    defined in Omixer format. Definitions of some modules are bundled with
    this package.

    Args:
        tables: Paths for tables of function counts to be processed. Tables
        should be tab separated, with functions on rows, observations on
        columns
        name: Name of the modules to be scored. Only used for formatting output
        output_dir: Path to directory to write results to
        modules: Name of bundled modules to use. If this is provided, the
        paths for definition, description, and hierarchy are ignored, and
        the bundled modules loaded instead.
        definition: Path to module definitions, in Omixer format
        description: Path to module description
        hierarchy: Pth to module hierarchy

    Returns:
        Produce several output files from RTK call: x.mat, x.desc, x.MODscore,
        x.KOused
    """
    logger.info("module_score")
    logger.info("|- Parameters -|")
    logger.info("Tables: %s tables", len(tables))
    logger.info("Name: %s", name)
    logger.info("Output: %s", output_dir)
    logger.info("Modules: %s", modules)
    logger.info("Definition: %s", definition)
    logger.info("Description: %s", description)
    logger.info("Hierarchy: %s", hierarchy)

    mods: RTKModules
    if modules is not None:
        logger.debug("Loading bundled modules: %s", modules)
        mods = RTKModules.from_bundled(modules)
        name = modules
    else:
        logger.debug("Loading custom modules")
        # Assume argparse types have done the required existence checks
        mods = RTKModules(
            definition=definition,
            description=description,
            hierarchy=hierarchy,
            name=name,
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    for tbl in tables:
        logger.debug("Processing table %s", tbl)
        with GzipManager(str(tbl)) as tbl_dc:
            tbl_name: str = ".".join(tbl_dc.name.split(".")[:-1])
            prefix_out: pathlib.Path = output_dir / f"{tbl_name}.{name}"
            logger.debug("Output prefix: %s", prefix_out)
            mods.score_table(pth_mgs_tbl=tbl_dc, output_prefix=prefix_out)


# Lazy testing - for some reason multimethod dispatch won't work properly
# from pytest
if __name__ == "__main__":
    # print(os.getcwd())
    # mgs: Dict[str, MGS] = create.read_mag_vs_gc(
    #     pathlib.Path("../../tests/data/DIME/MAGvsGC.txt.gz"),
    #     min_completeness=0,
    #     max_contamination=100
    # )
    #
    # annotation: pd.DataFrame = pd.read_csv(
    #     pathlib.Path("../../tests/data/DIME/DIAass_NOG.srt.gzgeneAss.gz"),
    #     sep="\t",
    #     header=None,
    #     index_col=0
    # )
    # description: pd.DataFrame = pd.read_csv(
    #     "../../tests/data/DIME/NOG.descr",
    #     sep="\t",
    #     header=None,
    #     index_col=0
    # )
    #
    # s: HierarchicalSchema = HierarchicalSchema(
    #     name="KEGG",
    #     hierarchy=description,
    #     gene_map=annotation
    # )
    # mgs_list: List[MGS] = list(mgs.values())[0:10]
    # m: MGS = mgs_list[0]
    #
    # # f = s.function_counts(mgs_list[0])
    # tbl: pd.DataFrame = s.function_counts(mgs_list)
    # foo = "bar"
    # cli_function_count(
    #     magvsgc=pathlib.Path("../../tests/data/DIME/MAGvsGC.txt.gz"),
    #     annotations=pathlib.Path(
    #         "../../tests/data/DIME/DIAass_CZy.srt.gzgeneAss.gz"
    #     ),
    #     descriptions=pathlib.Path("../../tests/data/DIME/Cazy.families.fix"),
    #     output=pathlib.Path("../../tests/output/"),
    #     compress=True,
    #     name="CAZY"
    # )
    logger.setLevel(logging.DEBUG)
    gbm_mods: RTKModules = RTKModules(
        definition="../../tests/data/module_defs/GBM/modGBM.list",
        description="../../tests/data/module_defs/GBM/modGBM.descr",
        hierarchy="../../tests/data/module_defs/GBM/modGBM_hiera.txt",
        name="GBM",
    )
    gbm_mods.score_table(
        pth_mgs_tbl=pathlib.Path("../../tests/output/KEGG/MGS.0.tsv"),
        output_prefix=pathlib.Path("../../tests/output/GBM_testo"),
    )
