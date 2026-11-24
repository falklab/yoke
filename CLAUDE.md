# CLAUDE.md

Guidance for AI agents (and people) working on yoke. Workflow details for humans are in
[CONTRIBUTING.md](CONTRIBUTING.md).

## Project
yoke is a Python package and CLI (`yoke <subcommand>`) for pangenomic and functional analysis of
metagenomes. It consumes outputs of the MG-TK / MATAFILER4 pipeline.

Domain terms used throughout the code:
- **MAG**: metagenome-assembled genome (one bin). Names look like `Sample1__SemiBin_177.fa.gz`.
- **MGS**: metagenomic species, a cluster of MAGs from the same species.
- **MAGvsGC**: MG-TK table with one row per MAG: MGS membership, completeness, contamination, and the
  MAG's genes as a comma-separated list (`,,` marks a contig boundary). `create.generate_mgs_from_magvsgc`
  streams it (always gzipped) and starts a new MGS whenever the MGS column changes, so the rows of an MGS
  must be adjacent. The gene list column is `other_Genes` (`DEF_MAGCVSGC_COLS`); MG-TK's Perl path writes
  `other_genes`, and `?` marks a MAG without MGS.
- **Function / schema**: functional annotation of genes, e.g. KEGG orthologs (KOs), eggNOG, CAZy. Gene
  annotation files have one function ID per row; `HierarchicalSchema` counts one function per gene. MG-TK
  `*geneAss` files can hold several functions in one field (`K00001,K00002`), which must be split first.
- **Modules**: sets of functions scored with the external `rtk` tool (`rtk module`; its options are listed
  in rtk's `options.cpp`, e.g. the hierarchy option is `-hiera`). Bundled sets in `src/yoke/data/`: BSB, the
  gut metabolic modules (`MF0001`-`MF0133`, Vieira-Silva et al. 2016), and GBM, the gut-brain modules
  (`MGB001`-`MGB056`, Valles-Colomer et al. 2019, with a hierarchy).

## Commands
```
uv sync                                   # install env from uv.lock
uv run pytest                             # tests
uv run pytest --cov                       # tests with coverage
uv run ruff check --fix . && uv run ruff format .
uv add <pkg> / uv add --group dev <pkg>   # dependencies (never edit uv.lock by hand)
```
Before declaring a change done, run `uv run ruff check .`, `uv run ruff format --check .` and
`uv run pytest`. All must pass.

## Layout and architecture
- `src/yoke/__init__.py`: the argparse CLI (`cli()`). Each subcommand has a subparser and dispatches via
  `set_defaults(func=lambda x: module.cli_<name>(arg=x.arg, ...))`.
- `src/yoke/annotations.py`: function counts, summarising to taxonomic ranks, module scoring (rtk),
  function enrichment.
- `src/yoke/create.py`: `MGS`, `MAG` and `Contig` classes, built from MAGvsGC.
- `src/yoke/extract.py`: writing MAG contigs to FASTA.
- `src/yoke/graphs.py`: pangenome graphs (networkx).
- `src/yoke/data/`: bundled module definitions, accessed with `importlib.resources` through
  `yoke.data.utils.BundledModules`.

Each subcommand follows the same pattern:
1. A `cli_<name>()` function in the feature module, which calls `_csv_field_adjust()`, sets the logger to
   INFO, logs a `|- Parameters -|` block, reads inputs, calls the core function and writes outputs.
2. A pure core function that takes and returns DataFrames (or objects), with no file I/O, so it can be
   tested directly.

The `function_enrichment` code in `annotations.py` (`function_enrichment`, `cli_function_enrichment`,
`_check_unique`, `_read_id_table`) and its tests are the best example of current conventions. Copy that
style for new code.

When adding or changing a subcommand, also update its section in `README.md`.

## Code style
ruff enforces formatting and linting (`pyproject.toml`). The formatter output is authoritative. Beyond what
ruff checks:

- **Layout**: 79 character lines, double quotes, 4-space indentation. Two blank lines between top-level
  definitions.
- **Imports**: absolute only (`from yoke.create import MAG`), never relative. Standard aliases: `numpy as
  np`, `pandas as pd`, `networkx as nx`, `plotnine as pn`. Use `pathlib` for paths (`import pathlib` and
  `pathlib.Path`, or `from pathlib import Path` in `extract.py` and tests, matching the file).
- **Type hints** on every function argument and return value, including `-> None`. Use built-in generics
  and `|` unions (`list[str]`, `dict[str, MGS]`, `pathlib.Path | None`), not `typing.List` or `Optional`.
  Use `from __future__ import annotations` in modules whose classes refer to themselves. Module-level
  variables are often annotated too (`logger: logging.Logger = ...`).
- **Docstrings**: Google style. The summary is one line ending with a period, on the same line as the opening
  `"""`. Then come `Args:`, `Returns:` and `Raises:` sections. Don't repeat types in `Args:`, because they
  are in the hints. Indent continuation lines by 4 spaces. Document constants and NamedTuple/class
  attributes with a string literal on the line after them:
  ```python
  ASSEMBLY_PATH: Path = Path("assemblies") / "metag"
  """Location of the assembly within a sample directory."""
  ```
- **Naming**: `snake_case` functions and variables, `CapWords` classes, keeping domain acronyms upper case
  (`MGS`, `MAG`, `RTKModules`). `UPPER_SNAKE` constants. Path variables are prefixed `pth_`
  (`pth_out`, `pth_mgs`). Alternative constructors are `from_*` static methods. Private module helpers
  use one leading underscore. Existing classes keep private state in `self.__attr` behind a
  read-only `@property`; follow that within those classes.
- **CLI**: subcommands and long options are `snake_case` (`--output_dir`). Use `type=pathlib.Path` for
  path arguments. Help and description text goes in parenthesised implicit string concatenation.
- **Logging**: use the module logger, `logger: logging.Logger = logging.getLogger(__name__)`, with lazy
  `%s` arguments (`logger.info("Write to %s", pth)`), not f-strings. Never use `print`.
- **Errors**: raise `ValueError(f"...")` (or another specific built-in) with a clear message for bad input.
  Use `raise ... from err` / `from None` when re-raising. Don't raise bare `Exception`. Validate cheap
  things before slow work.
- **pandas**: method chains go in parentheses, one call per line, each starting with `.`. Tabular files are
  tab-separated (`sep="\t"`), optionally gzip-compressed. Stream large MAGvsGC files with `csv.DictReader`.
- **Subprocesses**: new code passes an argument list to `subprocess.run` (no `shell=True`) and checks
  the return code.
- **Comments**: sparse, full-line, explaining why. Don't leave commented-out code.

## Tests
- pytest functions (no classes) in `tests/test_<module>.py`, named `test_<unit>_<behaviour>`, with a
  one-line docstring stating the expected behaviour.
- Fixtures use `@pytest.fixture` (never import from `_pytest`). Prefer small in-memory DataFrames.
  Write outputs to `tmp_path`.
- Never hard-code absolute or machine-specific paths. Use the `data_dir`, `expected_dir` and `local_data_dir`
  fixtures from `tests/conftest.py`. Small committed inputs go in `tests/data/`, golden outputs in
  `tests/expected/<test_name>/`.
- Mark tests needing large local data `@pytest.mark.local_data`, needing rtk `@pytest.mark.requires_rtk`,
  and taking more than a couple of seconds `@pytest.mark.slow`. Unmarked tests must run anywhere,
  quickly.
- Use `pytest.raises(..., match=...)`, `pytest.approx`, `pd.testing.assert_frame_equal` and `caplog` as
  appropriate. Every test must assert something.

## Git and releases
- Commit messages and PR titles use Conventional Commits (`feat(scope): ...`, `fix: ...`, `test: ...`,
  `docs: ...`, `refactor: ...`, `build: ...`, `ci: ...`, `chore: ...`). The scope is the subcommand or
  module. Use `!` for breaking CLI or output changes.
- Work on a branch, never directly on `main`.
- Don't change the version in `pyproject.toml`, `CHANGELOG.md` or `.release-please-manifest.json`;
  release-please manages them.
- Conda channels: anaconda.org is blocked on the team's network, so use the prefix.dev mirrors
  (`https://prefix.dev/conda-forge`, `https://prefix.dev/bioconda`).
