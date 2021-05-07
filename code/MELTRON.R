### calculate the melting score of long genes based on Genome Architecture Mapping (GAM) insulation scores.
### developed by Dominik Szabo (dominik.szabo@mdc-berlin.de) and Christoph J. Thieme (christoph.thieme@mdc-berlin.de)

library(tidyverse)

#set working directory to the directory code is in, to make sure that files are found
#for Rstudio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#for command line R:
#setwd(getSrcDirectory()[1])

#read in insulation scores (IS) at 10 different length scales (100 kb to 1000 kb in increments of 100 kb) in 4 cell types:
## mouse embryonic stem cells (ESC), oligodendroglia (OLG), dopaminergic neurons (DN) in two replicates, pyramidal glutamatergic neurons (PGN) in two replicates
tbl_is_scores <- read_tsv('../data/IS_scores_genomewide_100kbTO1mb.tsv.gz')

#read in list of long genes (GCF_000001635.26_GRCm38.p6_genomic.gtf)
gene_list <- read_tsv('../data/long_genes.tsv.gz')

#prepare the list and tibble in which the loop results will be written
scores <- list()
scores_red <- tibble(gene_id = NA, 
                     OLG_ks_pval = NA, DN_R1_ks_pval = NA, DN_R2_ks_pval = NA, PGN_R1_ks_pval = NA, PGN_R2_ks_pval = NA)


#loop through all long genes
for(i in 1:dim(gene_list)[1]){
  #store the IS values of each gene in the list
  scores[[gene_list[i,]$gene_id]]$ESC <- tbl_is_scores %>% dplyr::filter(chrom == as.character(gene_list[i,]$chrom),
                                                                                start >= gene_list[i,]$start_bin,
                                                                                stop <= gene_list[i,]$end_bin) %>% dplyr::select(matches("ESC"))
  scores[[gene_list[i,]$gene_id]]$OLG <- tbl_is_scores %>% dplyr::filter(chrom == as.character(gene_list[i,]$chrom),
                                                                                 start >= gene_list[i,]$start_bin,
                                                                                 stop <= gene_list[i,]$end_bin) %>% dplyr::select(matches("OLG"))
  scores[[gene_list[i,]$gene_id]]$DN_R1 <- tbl_is_scores %>% dplyr::filter(chrom == as.character(gene_list[i,]$chrom),
                                                                                 start >= gene_list[i,]$start_bin,
                                                                                 stop <= gene_list[i,]$end_bin) %>% dplyr::select(matches("DN_R1"))
  scores[[gene_list[i,]$gene_id]]$DN_R2 <- tbl_is_scores %>% dplyr::filter(chrom == as.character(gene_list[i,]$chrom),
                                                                                 start >= gene_list[i,]$start_bin,
                                                                                 stop <= gene_list[i,]$end_bin) %>% dplyr::select(matches("DN_R2"))
  scores[[gene_list[i,]$gene_id]]$PGN_R1 <- tbl_is_scores %>% dplyr::filter(chrom == as.character(gene_list[i,]$chrom),
                                                                                 start >= gene_list[i,]$start_bin,
                                                                                 stop <= gene_list[i,]$end_bin) %>% dplyr::select(matches("PGN_R1"))
  scores[[gene_list[i,]$gene_id]]$PGN_R2 <- tbl_is_scores %>% dplyr::filter(chrom == as.character(gene_list[i,]$chrom),
                                                                                 start >= gene_list[i,]$start_bin,
                                                                                 stop <= gene_list[i,]$end_bin) %>% dplyr::select(matches("PGN_R2"))
  #calculate p-values with a one-sided Kolmogorov-Smirnov test
    summary_statistics <- tibble(
      OLG_ks_pval=ks.test(unlist(scores[[gene_list[i,]$gene_id]]$ESC),
                            unlist(scores[[gene_list[i,]$gene_id]]$OLG),
                            alternative = "less")$p.value,
      DN_R1_ks_pval=ks.test(unlist(scores[[gene_list[i,]$gene_id]]$ESC),
                            unlist(scores[[gene_list[i,]$gene_id]]$DN_R1),
                            alternative = "less")$p.value,
      DN_R2_ks_pval=ks.test(unlist(scores[[gene_list[i,]$gene_id]]$ESC),
                            unlist(scores[[gene_list[i,]$gene_id]]$DN_R2),
                            alternative = "less")$p.value,
      PGN_R1_ks_pval=ks.test(unlist(scores[[gene_list[i,]$gene_id]]$ESC),
                            unlist(scores[[gene_list[i,]$gene_id]]$PGN_R1),
                            alternative = "less")$p.value,
      PGN_R2_ks_pval=ks.test(unlist(scores[[gene_list[i,]$gene_id]]$ESC),
                            unlist(scores[[gene_list[i,]$gene_id]]$PGN_R2),
                            alternative = "less")$p.value)
    
    scores_red <- add_row(scores_red, gene_id = gene_list[i,]$gene_id, summary_statistics)
}

save(scores, file='../data/IS_values_perGene.RData')

#drop the first empty row
scores_red_mod <- drop_na(scores_red) 

#run bonferroni multiple testing correction
scores_red_corrected <- scores_red_mod %>% 
  dplyr::mutate_if(is.double, ~ p.adjust(., method = "bonferroni", n=dim(scores_red_mod)[1] * (dim(scores_red_mod)[2] -1))) %>% 
  dplyr::mutate(meltingScore_OLG = -log10(OLG_ks_pval),
                meltingScore_DN_R1 = -log10(DN_R1_ks_pval),
                meltingScore_DN_R2 = -log10(DN_R2_ks_pval),
                meltingScore_PGN_R1 = -log10(PGN_R1_ks_pval),
                meltingScore_PGN_R2 = -log10(PGN_R2_ks_pval))

#plot denisty of melting scores to help decide on a melting threshold
scores_red_corrected %>% 
  pivot_longer(cols=c(meltingScore_OLG, meltingScore_DN_R1, meltingScore_DN_R2, meltingScore_PGN_R1, meltingScore_PGN_R2), names_to='cell_type', names_prefix='meltingScore_', values_to='melting_score') %>% 
  ggplot() + 
  geom_density(aes(x=melting_score, color=cell_type))+ 
  geom_vline(xintercept=5) # corresponds to a p-value of 1e-05 

#set melting score threshold to 5
scores_red_corrected_cat <- scores_red_corrected %>% 
  dplyr::mutate(melting_in_OLG = case_when(meltingScore_OLG > 5 ~ TRUE, TRUE~FALSE),
                melting_in_DN_R1 = case_when(meltingScore_DN_R1 > 5 ~ TRUE, TRUE~FALSE),
                melting_in_DN_R2 = case_when(meltingScore_DN_R2 > 5 ~ TRUE, TRUE~FALSE),
                melting_in_PGN_R1 = case_when(meltingScore_PGN_R1 > 5 ~ TRUE, TRUE~FALSE),
                melting_in_PGN_R2 = case_when(meltingScore_PGN_R2 > 5 ~ TRUE, TRUE~FALSE),
                DN_replicate_confirm = case_when(melting_in_DN_R1 == TRUE & melting_in_DN_R2 == TRUE ~ TRUE, TRUE~FALSE),
                PGN_replicate_confirm = case_when(melting_in_PGN_R1 == TRUE & melting_in_PGN_R2 == TRUE ~ TRUE, TRUE~FALSE))

write_tsv(scores_red_corrected_cat, '../data/melting_scores.tsv.gz')
