# Characterizing the informativeness of pathogen genome sequence datasets about transmission between population groups

## Abstract

Pathogen genome analysis helps characterizing transmission between population groups. The information carried by pathogen sequences comes from the accumulation of mutations within their genomes. This means that the pace at which mutations accumulate should determine the granularity of transmission processes that pathogen sequences can characterize. Here, we investigate how the complex interplay between mutation, transmission, mixing and sampling impacts the power of phylogeographic studies. First, we develop a conceptual probabilistic framework to quantify the ability of pairs of sequences of capturing migration history. This allows us to comprehensively explore the space of possible phylogeographic analyses by explicitly considering the pace at which mutations accumulate and the pace at which migration events occur. Using this framework, we identify a pathogen-intrinsic limit in the mixing scale at which their sequence data remains informative, with faster mutating pathogen enabling finer spatial characterization. Secondly, we perform a simulation study exploring a range of assumptions regarding sequencing intensity. We find that sample size further imposes a limit on the characterization of mixing processes. This work highlights inherent horizons of observability for population mixing processes that depend on the interaction between evolution, transmission, mixing and sampling. Such considerations are important for the design of phylogeographic studies. 

## Repository organization
This repository is organized in sub-folders as follows:
- ```input/``` contains input parameters values and data used to generate analyses. Further information is available on the folder-level [README](https://github.com/blab/phylogeo-signal/blob/main/input/README.md) file. 
- ```figures/``` contains the figures (both from the main text and the supplementary information) associated with the manuscript.
- ```manuscript/``` contains the manuscript.
- ```scripts/``` contains the code used to analyse the data and reproduce the figures. Further information is available on the folder-level [README](https://github.com/blab/phylogeo-signal/blob/main/scripts/README.md) file. 
- ```remaster/``` contains the code used to simulate the phylogenies depicted in Figure 1. 
- ```remaster-sample-size/``` contains the code used to perform the simulation study used for Figure 7 and 8.

## Install
The code is written in R and relies on some packages, which can be installed using:

```bash
Rscript ./scripts/install_requirements.R "scripts/requirements.txt"
```
