This project provides scripts to evaluate and apply MicroMiner.

To run all experiments regarding single mutations
```bash
sh run_mutation_experiments.sh
```






===================================

This project contains the scripts/libs for evaluating and applying my mutation
search tool using 3D micro-environments. It is used for the following purposes:
  1. Validate the 3D mutation search on mutation data for
     which both wild type and mutant have known protein structures in the PDB.
     It is evaluated whether the 3D mutation search can find these mutations
     retrospectively.
  2. Test how much more experimental structures representing mutants can be
     retrieved with the tool from the PDB.
  3. Annotate mutation entries in thermodynamic mutation data set without
     a known mutant structure with an experimental PDB structure using my
     mutation search tool.

My mutation search tool is based on the technologie of Bietz and Rarey:

<cite>ASCONA: Rapid Detection and Alignment of Protein Binding Site Conformations
Stefan Bietz and Matthias Rarey
Journal of Chemical Information and Modeling 2015 55 (8), 1747-1756
DOI: 10.1021/acs.jcim.5b00210</cite>

<cite>SIENA: Efficient Compilation of Selective Protein Binding Site Ensembles
Stefan Bietz and Matthias Rarey
Journal of Chemical Information and Modeling 2016 56 (1), 248-259
DOI: 10.1021/acs.jcim.5b00588</cite>

## Methods

### Mutation Datasets

Mutation data sets contain different data on experimentally performed mutations. Usually for the wild
type protein a PDB structure is available. In contrast, for the mutant a protein structure
is rarely available.

#### General Protein Mutations
<cite>M. D. Shaji Kumar, K. Abdulla Bava, M. Michael Gromiha, Ponraj Prabakaran, Koji Kitajima, Hatsuho Uedaira, 
Akinori Sarai, **ProTherm and ProNIT: thermodynamic databases for proteins and protein–nucleic acid interactions**, 
Nucleic Acids Research, Volume 34, Issue suppl_1, 1 January **2006**, Pages D204–D206,
https://doi.org/10.1093/nar/gkj103</cite>

<cite>Joicymara S Xavier, Thanh-Binh Nguyen, Malancha Karmarkar, Stephanie Portelli, Pâmela M Rezende, João P L Velloso,
David B Ascher, Douglas E V Pires,
**ThermoMutDB: a thermodynamic database for missense mutations**,
Nucleic Acids Research, Volume 49, Issue D1, 8 January **2021**, Pages D475–D479,
https://doi.org/10.1093/nar/gkaa925</cite>

<cite>Rahul Nikam, A Kulandaisamy, K Harini, Divya Sharma, M Michael Gromiha,
**ProThermDB: thermodynamic database for proteins and mutants revisited after 15 years**,
Nucleic Acids Research, Volume 49, Issue D1, 8 January **2021**, Pages D420–D424,
https://doi.org/10.1093/nar/gkaa1035</cite>

This is a TODO:
<cite>Jan Stourac, Juraj Dubrava, Milos Musil, Jana Horackova, Jiri Damborsky, Stanislav Mazurenko, David Bednar,
**FireProtDB: database of manually curated protein stability data**,
Nucleic Acids Research, Volume 49, Issue D1, 8 January **2021**, Pages D319–D324,
https://doi.org/10.1093/nar/gkaa981</cite>

#### Protein-Protein Interface mutations

<cite>Justina Jankauskaitė, Brian Jiménez-García, Justas Dapkūnas, Juan Fernández-Recio, Iain H Moal,
**SKEMPI 2.0: an updated benchmark of changes in protein–protein binding energy, kinetics and thermodynamics
upon mutation**, Bioinformatics, Volume 35, Issue 3, 01 February **2019**, Pages 462–469,
https://doi.org/10.1093/bioinformatics/bty635</cite>

#### Protein-Ligand Complex mutations

<cite>Douglas E.V. Pires, Tom L. Blundell, David B. Ascher,
**Platinum: a database of experimentally measured effects of mutations on structurally defined protein–ligand complexes**,
Nucleic Acids Research, Volume 43, Issue D1, 28 January **2015**, Pages D387–D391, https://doi.org/10.1093/nar/gku966</cite>


### Automatic protein similarity analysis

For all PDB structure pairs of wild type and mutant protein similarity can be
calculated using the TM-Align software. TM-Align is protein structure alignment tool. 
Usually we expect a wild type protein, and a mutant structure to be highly similar.
In the best case the structures should differ only at one position which
is the position of the mutation. Unfortunately, such ideal data representing
point mutations is rarely available. Therefore, in the general case, we expect 
the protein structure pairs to differ in multiple position.

Use `python add_tmalign.py` to add columns of protein similarity based on TM-Align
and it's TM-score. Usually protein structure pairs with a TM-score < 0.5 are 
likely unrelated in their overall folds.

### Manual evaluated PDB mutations

For the subset of mutations from the ProTherm and ThermoMutDB our tool
could not correctly identify we evaluated each mutation, and it's protein structure
pair manually (except for the ones with a TM-score < 0.5). Look at the `problematic_mutations.csv`
for a list of problematic mutations in the data sets with reason why they are
erroneous based on manual, visual examination of superposed structures.
