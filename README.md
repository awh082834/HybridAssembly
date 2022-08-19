HybridAssembly
==============

### Overview
HybridAssembly is written in Nextflow DSL2 with the purpose to perform hybrid assemblies using UniCycler as well as return metrics. This workflow also has the ability
to isolate contigs that are circularized by UniCycler and that are extra chromosomal in order to perform downstream analysis on those isolated contigs. 

### Inputs
|Flag     |Option     |Description             |
|:-------:|:---------:|:----------------------:|
|--out    |String of output directory|Name of an output directory that does not already exist|
|--sample |String used for the sample name|Will set the name of the sample used thoughout the workflow and in UniCycler|
|--isolate|true, false|Sets whether or not the workflow will isolate contigs from the UniCycler assembly that are extra chromosomal and circularized|

### Tools Used
1. Trimmomatic
2. Porechop
3. NanoPlot
4. FiltLong
5. UniCycler
6. QUAST

## Caveats 
Workflow needs to be run from a directory containing both the short and long reads. 