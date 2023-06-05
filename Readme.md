# MicroMiner Utils

MicroMiner searches similar local residue-centered 3D micro-environments in protein structure databases.

This repository provides auxiliary scripts and utility to annotate mutation datasets (like ProTherm)
with structures for the mutants.

MicroMiner can be used interactively at the [proteins.plus](https://proteins.plus) web server or
downloaded as part of the [NAOMI ChemBio Suite](https://uhh.de/naomi).

Single mutations extracted from the PDB can be downloaded [here](TODO). These data sets contain amino acid
pairs in protein structures exemplifying single mutations' local structural effects with multiple
similarity measurements, like RMSD: 
* single mutations in monomers (255 million pairs)
* single mutations in monomers (filtered) (~5 million pairs)
* single mutation at PPIs (~45 million pairs)
* single mutation at PPIs (filtered) (~800,000 pairs)

Mutations from and to non-standard residues can be also downloaded
in a separate file from the link above.

## Installation

Install Python dependencies for example with pip from the requirements.txt.

 ```bash
pip install -r requirements.txt
 ```

Download the [NAOMI suite](https://uhh.de/naomi) to get the MicroMiner executable and package. 
Follow instructions in the MicroMiner package readme file.

Set paths and parameters in [config.ini](config.ini) accordingly. 

Verify that everything works by running the unittests:

```bash
python -m unittest discover helper
```

## Data sets of mutation effect measurements

If you want to use mutation effect datasets you need to download
them from their respective sites and set the paths in config.ini: 

* protein mutations: [ProTherm](https://doi.org/10.1093/nar/gkj103),
[ThermoMutDB](https://doi.org/10.1093/nar/gkaa925),
[ProThermDB](https://doi.org/10.1093/nar/gkaa1035),
[FireProtDB](https://doi.org/10.1093/nar/gkaa981)
* protein-protein complex mutations:
[SKEMPI2](https://doi.org/10.1093/bioinformatics/bty635)
* protein-ligand complex mutations:
[Platinum](https://doi.org/10.1093/nar/gku966)
* structural effects of mutations: [Shanthirabalan et al.](https://doi.org/10.1002/prot.25499) 

Note that some of these datasets like SKEMPI2 and Platinum come with custom prepared PDB files.

## Reproduce experiments from the paper

To reproduce the results from the paper run the following commands.

### single mutation benchmark
```bash
sh run_mutation_benchmark.sh
```
Then run the notebooks `plot_single_mutation_benchmark.ipynb` and `protein_flexibility.ipynb` to 
obtain the plots from the paper.

### annotate mutation effect measurements with structures for the mutant

```bash
sh run_mutation_annotation.sh
python annotate_mutation_datasets.py -d protherm thermomutdb platinum prothermdb skempi2 fireprotdb -m results/mutation_annotation/search -o .
python annotation_statistics.py -d protherm thermomutdb platinum prothermdb skempi2 fireprotdb -m results/mutation_annotation/search -o .
```
Then run the notebook `annotation_analysis.ipynb` to generate the plot.

### single mutations in the PDB

Set `run_on_hpc=false` in `run_pdb_experiments.sh` to run on the local machine. Note that this 
experiment is computationally expensive, and we optimized the code for our in-house HPC (SGE).
It is not expected that the HPC code works out of the box in other environments.

```bash
sh run_pdb_experiments.sh

# collect results
head -n1 results/pdb_experiments/search/pdb/results/100D/resultStatistic.csv > pdb_all_ppi.tsv && find results/pdb_experiments/search/pdb/results -name resultStatistic.csv -exec tail -n+2 -q {} >> pdb_all_ppi.tsv +
```

Then run the `pdb_stats_database.ipynb` notebook.
