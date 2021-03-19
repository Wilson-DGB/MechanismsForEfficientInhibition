# MechanismsForEfficientInhibition
This repository contains primarily julia code for a model of receptor inhibition.

There are currently (as of 19/03/2021) six main scripts:

1 CaseBoth_WellMixedLigandDimerisation.jl
2 CaseNeither_WellMixedLigandDimerisation.jl
3 CaseCD28_WellMixedLigandDimerisation.jl
4 CasePD1_WellMixedLigandDimerisation.jl
5 CaseNeither_SpatialLigandDimerisation.jl
6 VisualisationCode.jl

Codes 1 through 4 implement a model of inhibition where the ligands are well-mixed and associate to the various receptors through first order rates. The two ligands can dimerise and the product can either activate both costimulatory and coinhibitory receptors (1), activate neither pathways (2), activate just the costimulatory pathway (3), or just the coinhibitory pathway (4). 

Code 5 realises the ligand molecules spatially. Currently only the case where the dimerised product activates neither pathway is available.

Code 6 presents code that visualises the position of receptors and ligands in real time as the simulation progresses.

There is additionally a folder storing all pretabulated reaction rates for the CRDME and a second folder with some matlab scipts useful for plotting data from .csv files

