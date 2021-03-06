# MELTRON 

_a statistical framework to detect differences in chromatin contact density at genomic regions of interest_

[![DOI:10.1038/s41586-021-04081-2](https://zenodo.org/badge/DOI/110.1038/s41586-021-04081-2.svg)](https://www.nature.com/articles/s41586-021-04081-2)

<img src="data/IS_gif.gif" width="350">

<img src="./data/meltron_pipeline.png" width="900">

Deposited scripts and files for calculation of [melting scores](https://www.nature.com/articles/s41586-021-04081-2) of long genes in oligodendroglia (OLG) from the somatosensory cortex, dopaminergic neurons (DNs) from the midbrain VTA and pyramidal glutamatergic neurons (PGNs) from the hippocampus CA1 relative to embryonic stem cells (ESCs).

### Required packages
```r
library(tidyverse)
library(ggpubr)
```

### Available scripts
- code/MELTRON.R:   
   Compares insulation score (IS) distributions over long genes in OLGs, DNs, PGNs to ESCs and calculates melting score per gene.   
- code/plot_ecdf.R:  
   Plots empirical cumulative density functions (ECDF) for IS values of individual genes.   
- code/plot_rna_atac_over_meltingScore.R:  
   Plots expression and chromatin accessibility as a function of the melting score per cell-type. Density plots with median lines indicate population trends. 
- code/plot_domain_melting_gene_characteristics.R 
   Plots heatmaps and summary violin plots for each of the clusters.


### Available command line scripts:
- command_line_apps/matrix_wide_to_long.R: 
    Converts a square matrix into a long matrix for IS calculation. Accepts wildcards for processing of multiple chromosomes.
- command_line_apps/long_matrix_to_IS.R:
    Calculates insulation scores at multiple distances (default 100kb - 1Mb, steps of 100kb).
- command_line_apps/Meltron.R:
    Calculation of melting scores over genomic regions of interest
    
type
```bash
Rscript command_line_apps/matrix_wide_to_long.R --help 
```
for explanations
   
   
Developed and tested with R version 3.6.0 Planting of a Tree.  
Developed and maintained by Dominik Szabó [<img src="https://cloud.githubusercontent.com/assets/1810515/4228292/6b03dc88-3958-11e4-9094-d3c1771ccfea.png" width="15">](https://orcid.org/0000-0001-8109-5088) with intellectual input from Christoph Thieme [<img src="https://cloud.githubusercontent.com/assets/1810515/4228292/6b03dc88-3958-11e4-9094-d3c1771ccfea.png" width="15">](https://orcid.org/0000-0002-1566-0971).  
Please get in touch for questions and issues: dominik.szabo at mdc-berlin.de  


Please check our manuscript:  
__Winick-Ng, W., Kukalev, A., Harabula, I., Zea Redondo, L., Szabó, D. et al.:  
[Cell-type specialization is encoded by specific chromatin topologies](https://www.nature.com/articles/s41586-021-04081-2)__



