# Yoke
A python package with CLI tools for doing assorted bits of pangenomic and functional analysis with MG-TK outputs.

## Installation
There are currently some problems building recent `scipy` with the versions given in the current `toml`.
I will resolve these eventually, but for now, what works is to make an environment with the following
version
```micromamba create --name yoke pip python==3.12 cmake gfortran libopenblas scipy==1.15.1```
then
```pip install git+https://github.com/falklab/yoke.git```

## Update
I'm not doing proper versioning currently.
To force update, use:
```pip install --upgrade --force-reinstall git+https://github.com/falklab/yoke.git```

### Install in a conda/mamba env
To install in a conda or mamba environment, make sure you have pip installed in the environment rather than using
the system pip, then install as above.

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