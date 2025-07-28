# Simulation study used to depict the phylogenetic trees used in Figure 1

We provide the `Snakefile` we used to run workflow.
Running the simulations requires the installation of BEAST2 (version 2.7 or later) along with the `ReMASTER` and `feast` packages.
For the workflow to work, make sure the BEAST2 file path in the `simulate_from_xml` rule refers to the local BEAST2 application path. 

The simulations additionally rely on three additional R packages which can be installed as part of a conda environment using

```
conda env create -f envs/remaster.yaml
```

The environment can then be activated using 
```
conda activate remaster
```
and deactivated using 
```
conda deactivate
```

It is then possible to run the worflow using (example using 1 thread):
```
snakemake --use-conda --cores 1
```
