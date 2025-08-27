# Yoke
A python package with CLI tools for doing assorted bits of pangenomic and functional analysis with MG-TK outputs.

## Installation
```pip install git+https://github.com/apduncan/yoke.git```

## Update
I'm not doing proper versioning currently.
To force update, use:
```pip install --upgrade --force-reinstall git+https://github.com/apduncan/yoke.git```

### Install in a conda/mamba env
To install in a conda or mamba environment, make sure you have pip installed in the environment rather than using
the system pip, then install as above.

## CLI Tools
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

