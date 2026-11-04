# Yoke
A python package with CLI tools for doing assorted bits of pangenomic and functional analysis with MG-TK outputs.

## Installation
Yoke requires Python 3.12 or later.
Releases are tagged `vX.Y.Z`; the available versions and their changes are listed on the
[releases page](https://github.com/falklab/yoke/releases) and in `CHANGELOG.md`.

Install a specific release with pip
```
pip install git+https://github.com/falklab/yoke.git@vX.Y.Z
```
or as a standalone command line tool with [uv](https://docs.astral.sh/uv/)
```
uv tool install git+https://github.com/falklab/yoke.git@vX.Y.Z
```
Omitting `@vX.Y.Z` installs the latest development version from `main`.

### Install in a conda/mamba env
To install in a conda or mamba environment, make sure you have pip installed in the environment rather than using
the system pip, then install as above.
```
micromamba create --name yoke --override-channels \
    -c https://prefix.dev/conda-forge -c https://prefix.dev/bioconda \
    python=3.12 pip rtk
micromamba activate yoke
pip install git+https://github.com/falklab/yoke.git@vX.Y.Z
```

### External tools
`yoke module_score` uses [rtk](https://github.com/hildebra/rtk), which must be available on your `PATH`.
It can be installed from bioconda, as in the example above.
The examples use the [prefix.dev](https://prefix.dev) mirrors of conda-forge and bioconda, as anaconda.org is
blocked on some institute networks.

## Update
To move to a different release, install again with the new tag
```
pip install --upgrade git+https://github.com/falklab/yoke.git@vX.Y.Z
```
`yoke --version` shows the installed version.

## Development
See [CONTRIBUTING.md](CONTRIBUTING.md) for setting up a development environment, code style, testing, and
the branching and release process.

## CLI Tools
All of the CLI tools are accessible via `yoke`.
Run `yoke` to see a list of available tools.

Each individual tool will have its own help page, accessible via 
`yoke <tool> --help`.

Some of the tools are described below as well.

### Function Counts
You can produce tables with counts of how many times each function appears in each MAG, grouped by MGS.
A single table per MGS will be produced.
For help on usage

```yoke function_counts --help```

Example usage for KEGG ortholog counts

```
yoke function_counts \
-m Bin_SB/MAGvsGC.txt.gz \
-a Anno/Func/DIAass_KGM.srt.gzgeneAss.gz \
-d Anno/Func/KEGG.descr \
-o Bin_SB/function_counts \
-n KEGG \
--compress
```

You can then summarise these results to different taxonomic levels using `yoke function_summarise`.
See `yoke function_summarise --help` for usage.

### Module Scores
You can produce a table of module scores for the MAGs in a MGS.
There are a number of bundled modules which can be selected, or you can provide your own.
The input is the table of functional annotations per MAG produced by `function_counts`.

Example usage for scoring BSB gut-brain modules, which are based on KEGG orthologs.
```
yoke module_score \
-i KEGG/MGS.* \
-o module_scores/ \
-m BSB \
-n BSB
```

For each MGS, four output files will be produced:
- `MGS.1.{mod}.mat`: Matrix of counts of genes in each module for each MAG in the MGS.
- `MGS.1.{mod}.MODscore`: Module scores for each MAG in the MGS.
- `MGS.1.{mod}.KOused`: Which functions participated in each module for each MAG in the MGS.
- `MGS.1.{mod}.descr`: Decription of each module which appears in the MGS.

The bundled modules are:
- BSB: Human gut-brain modules. Based on KEGG orthologs.
- GBM: Human gut metabolism modules. Based on KEGG orthologs.
  ([Vieira-Silva et al., 2016](https://doi.org/10.1038/nmicrobiol.2016.88))

### Function Enrichment
Test whether the presence of each function is associated with a metadata category which has two levels
(e.g. case/control, or a True/False trait).
For each function, a 2x2 table of presence/absence in each level is tested using Fisher's, Barnard's or
Boschloo's exact test (from `scipy`), and p-values are corrected for multiple testing within each metadata
column and test (any method from `statsmodels.stats.multitest`, default Benjamini-Hochberg).

Inputs are:
- One or more function count tables, with functions on rows and MGS or MAGs on columns, such as those
  produced by `function_counts` or `function_summarise`. Values greater than 0 count as present. Multiple
  tables are combined, with functions missing from a table treated as absent.
- A tab separated metadata table, with MGS or MAG identifiers in the first column and one column per
  category to test. NAs are removed per column, and columns without exactly two levels are skipped.

Example usage
```
yoke function_enrichment \
-i KEGG/MGS.*.tsv.gz \
-m mag_metadata.tsv \
-o enrichment/KEGG.tsv \
-t fisher barnard \
-c disease
```

The output is indexed by metadata column and function, with `{test}_stat`, `{test}_p` and
`{test}_p_adj` columns for each test. For Fisher's test the statistic is the odds ratio of presence in the
first level (alphabetically, or `True` for boolean columns) relative to the second.

### Extract Contigs
Produce FASTA files with the contigs for each MAG.
This can be useful if you want to use MAGs for some analysis that needs a whole(ish)
genome context.
I have tried to make sensible decision to minimise runtime: all MAGs from one assembly
get written during the same linear scan through the assembly, and multiprocessing is
supported to handle multiple assemblies simultaneously (though I'm not sure it actually 
improves speed).
Output is one gzipped fasta file per MAG, with name `{mag_name}.contigs.fna.gz`.

When you have run multiple binners on the same assemblies, please ensure you specify
`--binner`, to ensure the correct binner results are used.

Example usage
```
yoke extract_contigs \
-m magvsgc/MAGvsGC.txt.gz \
-a assemblies/ \
-o output/ \
--completeness 50 \
--contamination 10 \
--binner SB \
--processes 4
```