# MELTRON 

_a statistical framework to detect differences in chromatin contact density at genomic regions of interest_

[![DOI:10.1101/2020.04.02.020990](http://img.shields.io/badge/DOI-10.1101/2020.04.02.020990-B31B1B.svg)](https://www.biorxiv.org/content/10.1101/2020.04.02.020990v1)

<img src="data/IS_gif.gif" width="400">

<img src="./data/meltron_pipeline.png" width="900">

Deposited scripts and files allow to calculate the melting score of long genes in oligodendroglia (OLG) from the somatosensory cortex, dopaminergic neurons (DNs) from the midbrain VTA and pyramidal glutamatergic neurons (PGNs) from the hippocampus CA1 relative to embryonic stem cells (ESCs).

<img src="./data/melting_examples.png" width="900">

### Required packages
```r
library(tidyverse)
library(lemon)
library(ggpubr)
```

### Available scripts
- domain_melting/MELTRON.R:   
   Compares insulation score (IS) distributions over long genes in OLGs, DNs, PGNs to ESCs and calculates melting score per gene.   
- domain_melting/plot_ecdf.R:  
   Plots empirical cumulative density functions (ECDF) for IS values of individual genes.   
- domain_melting/plot_rna_atac_over_meltingScore.R:  
   Plots expression and chromatin accessibility as a function of the melting score per cell-type. Density plots with median lines indicate population trends.   
- trans_cis_ratio/calculate_trans_cis_ratio.R:  
   Calculates trans-cis NPMI ratios for GAM data genome-wide for all cell-types and replicates.   
- trans_cis_ratio/plot_trans_cis_raincloud.R:  
   Plots trans-cis ratio of long genes in the different cell-types stratified by melting.   


Developed and tested with R version 3.6.0 Planting of a Tree

Command line tool that calulates melting scores for user-specified input is under development


Please check our [preprint](https://www.biorxiv.org/content/10.1101/2020.04.02.020990v1):

__Winick-Ng, W., Kukalev, A., Harabula, I., Zea Redondo, L., Szabo, D.:  
Cell-type specialization in the brain is encoded by specific long-range chromatin topologies__



