##plotting of heatmaps, summary violin and bar plots displayed in ED figures 5f-i based on supplementary Table 8 

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
`%notin%` <- Negate(`%in%`)
#set working directory to the directory code is in, to make sure that files are found
#for Rstudio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#for command line R:
#setwd(getSrcDirectory()[1])

melting_scores <- read_tsv('../data/master_table_long_genes_withClustIdent.tsv.gz') %>% 
  dplyr::mutate(DN_R1_comp_transition = paste0(ESC_compartment, ' -> ', DN_R1_compartment),
                DN_R2_comp_transition = paste0(ESC_compartment, ' -> ', DN_R2_compartment),
                PGN_R1_comp_transition = paste0(ESC_compartment, ' -> ', PGN_R1_compartment),
                PGN_R2_comp_transition = paste0(ESC_compartment, ' -> ', PGN_R2_compartment),
                OLG_comp_transition = paste0(ESC_compartment, ' -> ', OLG_compartment),
                rnaRPMM_l2FC_DNs = log2(rnaRPMM_DN) - log2(rnaRPMM_ESC),
                rnaRPMM_l2FC_PGNs = log2(rnaRPMM_PGN) - log2(rnaRPMM_ESC),
                rnaRPMM_l2FC_OLGs = log2(rnaRPMM_OLGs) - log2(rnaRPMM_ESC),
                atacRPMM_l2FC_DNs = log2(atacRPMM_DN) - log2(atacRPMM_ESC),
                atacRPMM_l2FC_PGNs = log2(atacRPMM_PGN) - log2(atacRPMM_ESC),
                atacRPMM_l2FC_OLGs = log2(atacRPMM_OLGs) - log2(atacRPMM_ESC)) %>% 
  dplyr::mutate_at(c('DN_R1_comp_transition', 'DN_R2_comp_transition', 'PGN_R1_comp_transition', 'PGN_R2_comp_transition', 'OLG_comp_transition'), factor, levels=c('B -> B', 'A -> B', 'B -> A', 'A -> A'))


#ED figure 5f: topoisomerase inhibition sensitivity
scores_for_plotting <- melting_scores %>% 
  dplyr::mutate(vtaR1_1000_top15 = case_when(gene_symbol %in% (top_n(melting_scores, 15, meltingScore_OLGs))$gene_symbol ~ T, T~F),
                vtaR2_1000_top15 = case_when(gene_symbol %in% (top_n(melting_scores, 15, meltingScore_DN_R1))$gene_symbol ~ T, T~F),
                ca1R1_1000_top15 = case_when(gene_symbol %in% (top_n(melting_scores, 15, meltingScore_DN_R2))$gene_symbol ~ T, T~F),
                ca1R2_1000_top15 = case_when(gene_symbol %in% (top_n(melting_scores, 15, meltingScore_PGN_R1))$gene_symbol ~ T, T~F),
                oligo_1000_top15 = case_when(gene_symbol %in% (top_n(melting_scores, 15, meltingScore_PGN_R2))$gene_symbol ~ T, T~F),
                #vta_1000_replicate_confirm = case_when(vtaR1_1000_top15 == T & vtaR2_1000_top15 == T ~ T, T~F),
                #ca1_1000_replicate_confirm = case_when(ca1R1_1000_top15 == T  & ca1R2_1000_top15 == T ~ T, T~F),
                top15_inany = case_when(vtaR1_1000_top15 == T | vtaR2_1000_top15 == T | ca1R1_1000_top15 == T | ca1R2_1000_top15 == T | oligo_1000_top15 == T ~ T, T~F),
                melting_below_top = case_when(top15_inany == F & (meltingScore_OLGs > 5 | meltingScore_DN_R1 > 5 | meltingScore_DN_R2 > 5 | meltingScore_PGN_R1 > 5 | meltingScore_PGN_R2 > 5) ~ T, T~F))

# how many of the melting genes are topo sensitive:
threebar <- scores_for_plotting %>% group_by(top15_inany, melting_below_top, topoisomerase_inhibit_sens) %>% summarise(count=n()) %>% 
  dplyr::mutate(percent = count / sum(count),
                melting_status = case_when(top15_inany==FALSE & melting_below_top == FALSE ~ 'not_melting',
                                           top15_inany==TRUE  ~ 'top3',
                                           melting_below_top == TRUE ~ 'intermediate'),
                melting_status = factor(melting_status, levels=c('top3', 'intermediate', 'not_melting')))
threebar %>% 
  ggplot() + 
  geom_col(aes(x=melting_status, y=percent, fill=topoisomerase_inhibit_sens)) + 
  scale_fill_manual(values = c('grey90', 'grey30')) + 
  #  scale_x_discrete(limits=c('TRUE', 'FALSE')) + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12))

#stats:
chisqtable <- scores_for_plotting %>% 
  dplyr::mutate(melting_status = case_when(top15_inany==FALSE & melting_below_top == FALSE ~ 'not_melting',
                                           top15_inany==TRUE  ~ 'top3',
                                           melting_below_top == TRUE ~ 'intermediate'),
                melting_status = factor(melting_status, levels=c('top3', 'intermediate', 'not_melting')),
                topoisomerase_inhibit_sens = case_when(topoisomerase_inhibit_sens ==TRUE ~ 'sensitive', 
                                                       topoisomerase_inhibit_sens ==FALSE ~ 'insensitive')) %>% 
  dplyr::select(topoisomerase_inhibit_sens, melting_status)
table <- table(chisqtable$melting_status, chisqtable$topoisomerase_inhibit_sens)
chisq.test(table)


#ED5g heatmaps per cell type
#universal heatmap colors:
col_runif_rna = colorRamp2(c(0.4, 2, 3.4), c('blue', "white", "red"))
col_runif_atac = colorRamp2(c(2,2.7,3.2), c('blue', "white", "red"))
col_black_white = colorRamp2(c(0,100), c('white', 'grey30'))
col_transcis <- colorRamp2(c(-100,0,100), c('blue', "white", "red"))
col_discr = structure(c('#009E73', '#F0E442'), names = c("interior", "periphery")) # black, red, green, blue

#OLG
esc_olig <- melting_scores %>% 
  dplyr::filter(meltingScore_OLGs > 5) %>% 
  dplyr::arrange( ESC_compartment, OLG_compartment)
esc_olig <- rowid_to_column(esc_olig)

col_runif_melting_olig = colorRamp2(c(0, 107), c("white", "#800080"))
plot_data1 <- as.matrix(esc_olig[,c('ESC_compartment', 'OLG_compartment') ])
plot_data2 <- as.matrix(esc_olig[,c('meltingScore_OLGs') ])
plot_data3 <- log10(as.matrix(esc_olig[,c('rnaRPMM_ESC', 'rnaRPMM_OLGs')]))
plot_data4 <-  as.matrix(esc_olig[,c('percent_genebody_lad', 'percent_genebody_nad')])
plot_data4[is.na(plot_data4)] <- 0
plot_data5 <- log10(as.matrix(esc_olig[,c('atacRPMM_ESC', 'atacRPMM_OLGs')]))

rownames(plot_data2) <- esc_olig$gene_symbol
compartments_h <- Heatmap(plot_data1, col=c('#9E31DA', '#929292',  '#F49D1F'), heatmap_legend_param=list(title='Compartment'), cluster_rows = FALSE)
meltingscore <- Heatmap(plot_data2, col=col_runif_melting_olig, heatmap_legend_param=list(title='Melting score'),cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontsize = 9))
rpkm_h <- Heatmap(plot_data3, col=col_runif_rna, heatmap_legend_param=list(title='RPMM'), name='rpmm', cluster_columns = FALSE, row_km = 4,  row_km_repeats = 1000, cluster_row_slices = FALSE, cluster_rows = FALSE)
loc_pos <-  Heatmap(plot_data4, col=col_black_white, cluster_columns=FALSE, cluster_rows=FALSE, heatmap_legend_param=list(title='periphery assoc [%]'))
atac_h <- Heatmap(plot_data5, col=col_runif_atac, heatmap_legend_param=list(title='atac RPMM'), cluster_columns = FALSE, cluster_rows = FALSE)
ht_list <- rpkm_h + atac_h + compartments_h + loc_pos + meltingscore 
draw(ht_list, ht_gap=unit(0, 'cm'), merge_legends = TRUE, heatmap_legend_side = "bottom")
dev.off()

olig_toplot <- melting_scores %>% 
  drop_na(OLG_cluster_ident) %>% 
  dplyr::mutate(OLG_cluster_ident = as.character(OLG_cluster_ident))

#log2FC of RNA per cluster
olig_toplot %>% 
  ggplot()+
  geom_violin(aes(x=OLG_cluster_ident, y=rnaRPMM_l2FC_OLGs))+
  geom_boxplot(aes(x=OLG_cluster_ident, y=rnaRPMM_l2FC_OLGs), outlier.shape = NA, width=0.1)+
  geom_jitter(aes(x=OLG_cluster_ident, y=rnaRPMM_l2FC_OLGs), size=0.5, width=0.2)+ 
  geom_hline(yintercept = 0, linetype='dotted')
#log2FC of ATAC per cluster
olig_toplot %>% 
  ggplot()+
  geom_violin(aes(x=OLG_cluster_ident, y=atacRPMM_l2FC_OLGs))+
  geom_boxplot(aes(x=OLG_cluster_ident, y=atacRPMM_l2FC_OLGs), outlier.shape = NA, width=0.2)+
  geom_jitter(aes(x=OLG_cluster_ident, y=atacRPMM_l2FC_OLGs), size=0.5, width=0.2)+
  geom_hline(yintercept = 0, linetype='dotted')
#melting score per cluster
olig_toplot %>% 
  ggplot()+
  geom_violin(aes(x=OLG_cluster_ident, y=meltingScore_OLGs))+
  geom_boxplot(aes(x=OLG_cluster_ident, y=meltingScore_OLGs), outlier.shape = NA, width=0.1)+
  geom_jitter(aes(x=OLG_cluster_ident, y=meltingScore_OLGs), size=0.5, width=0.2)
#compartment transition
olig_toplot %>% 
  group_by(OLG_cluster_ident, OLG_comp_transition) %>% 
  summarise(n=n()) %>%
  dplyr::filter(!is.na(OLG_comp_transition)) %>% 
  dplyr::mutate(percent= n/ sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=OLG_cluster_ident, y=percent, fill=OLG_comp_transition), position = 'dodge')+
  labs(fill = "transition")+
  scale_fill_manual(values=c('A -> A' = '#9E31DA', 'A -> B' ='#825515','B -> A' =   '#c88ee8', 'B -> B' = '#F49D1F'))


#PGNs: 
esc_ca1 <- melting_scores %>% 
  dplyr::filter(meltingScore_PGN_R1 > 5|  meltingScore_PGN_R2 >5) %>% # | meltingScore_PGN_R2 >5
  dplyr::mutate(DN_cluster_ident = as.character(DN_cluster_ident)) %>% 
  dplyr::arrange(ESC_compartment, PGN_R1_compartment, PGN_R2_compartment)
esc_ca1 <- rowid_to_column(esc_ca1)

col_runif_melting_ca1 = colorRamp2(c(0, 90), c("white", "#6367DC"))
col_runif_rna_ca1 = colorRamp2(c(0.4, 2, 3.6), c('blue', "white", "red"))
col_runif_atac_ca1 = colorRamp2(c(2,2.7,3.4), c('blue', "white", "red"))
plot_data1 <- as.matrix(esc_ca1[,c('ESC_compartment', 'PGN_R1_compartment', 'PGN_R2_compartment') ])
plot_data2 <- as.matrix(esc_ca1[,c('meltingScore_PGN_R1', 'meltingScore_PGN_R2') ])
plot_data3 <- log10(as.matrix(esc_ca1[,c('rnaRPMM_ESC', 'rnaRPMM_PGN')]))
plot_data4 <- as.matrix(esc_ca1[,c('percent_genebody_lad', 'percent_genebody_nad')])
plot_data4[is.na(plot_data4)] <- 0
plot_data5 <- log10(as.matrix(esc_ca1[,c('atacRPMM_ESC', 'atacRPMM_PGN')]))

rownames(plot_data2) <- esc_ca1$gene_symbol
compartments_h <- Heatmap(plot_data1, col=c('#9E31DA', '#929292',  '#F49D1F'), heatmap_legend_param=list(title='Compartment'), cluster_rows=FALSE)
#compartments_ca1 <- Heatmap(plot_data15, col=c('#9E31DA', '#929292',  '#F49D1F'), heatmap_legend_param=list(title='Compartment')) #, width = unit(2, "cm")
meltingscore <- Heatmap(plot_data2, col=col_runif_melting_ca1, heatmap_legend_param=list(title='Melting score'),cluster_columns = FALSE, cluster_rows=FALSE, row_names_gp = gpar(fontsize = 7))
rpkm_h <- Heatmap(plot_data3, col=col_runif_rna_ca1, heatmap_legend_param=list(title='RPMM'), name='rpmm', cluster_columns = FALSE, row_km = 4,  row_km_repeats = 1000, cluster_row_slices = FALSE, cluster_rows = FALSE)
loc_pos <- Heatmap(plot_data4, col=col_black_white, cluster_rows=FALSE, cluster_columns = FALSE,  heatmap_legend_param=list(title='periphery assoc [%]'))
atac_h <- Heatmap(plot_data5, col=col_runif_atac_ca1, heatmap_legend_param=list(title='atac RPKM'), cluster_columns = FALSE, cluster_rows=FALSE)
ht_list <- rpkm_h + atac_h + compartments_h + loc_pos + meltingscore 
draw(ht_list, ht_gap=unit(0, 'cm'), merge_legends = TRUE, heatmap_legend_side = "bottom")
dev.off()

ca1_toplot <- melting_scores %>% 
  drop_na(PGN_cluster_ident) %>% 
  dplyr::mutate(PGN_cluster_ident = as.character(PGN_cluster_ident))

#log2FC of RNA per cluster
ca1_toplot %>% 
  ggplot()+
  geom_violin(aes(x=PGN_cluster_ident, y=rnaRPMM_l2FC_PGNs))+
  geom_boxplot(aes(x=PGN_cluster_ident, y=rnaRPMM_l2FC_PGNs), outlier.shape = NA, width=0.1)+
  geom_jitter(aes(x=PGN_cluster_ident, y=rnaRPMM_l2FC_PGNs), size=0.5, width=0.2)+ 
  geom_hline(yintercept = 0, linetype='dotted')
#log2FC of ATAc per cluster
ca1_toplot %>% 
  ggplot()+
  geom_violin(aes(x=PGN_cluster_ident, y=atacRPMM_l2FC_PGNs))+
  geom_boxplot(aes(x=PGN_cluster_ident, y=atacRPMM_l2FC_PGNs), outlier.shape = NA, width=0.2)+
  geom_jitter(aes(x=PGN_cluster_ident, y=atacRPMM_l2FC_PGNs), size=0.5, width=0.2)+
  geom_hline(yintercept = 0, linetype='dotted')

#melting score per cluster
ca1_toplot %>% 
  ggplot()+
  geom_violin(aes(x=PGN_cluster_ident, y=meltingScore_PGN_R1))+
  geom_boxplot(aes(x=PGN_cluster_ident, y=meltingScore_PGN_R1), outlier.shape = NA, width=0.1)+
  geom_jitter(aes(x=PGN_cluster_ident, y=meltingScore_PGN_R1), size=0.5, width=0.2)
ca1_toplot %>% 
  ggplot()+
  geom_violin(aes(x=PGN_cluster_ident, y=meltingScore_PGN_R2))+
  geom_boxplot(aes(x=PGN_cluster_ident, y=meltingScore_PGN_R2), outlier.shape = NA, width=0.1)+
  geom_jitter(aes(x=PGN_cluster_ident, y=meltingScore_PGN_R2), size=0.5, width=0.2)

#compartment transition
ca1_toplot %>% 
  group_by(PGN_cluster_ident, PGN_R1_comp_transition) %>% 
  summarise(n=n()) %>%
  dplyr::filter(!is.na(PGN_R1_comp_transition)) %>% 
  dplyr::mutate(percent= n/ sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=PGN_cluster_ident, y=percent, fill=PGN_R1_comp_transition), position = 'dodge')+
  labs(fill = "transition")+
  scale_fill_manual(values=c('A -> A' = '#9E31DA', 'A -> B' ='#825515','B -> A' =   '#c88ee8', 'B -> B' = '#F49D1F'))
ca1_toplot %>% 
  group_by(PGN_cluster_ident, PGN_R2_comp_transition, .drop=FALSE) %>% 
  summarise(n=n()) %>%
  dplyr::filter(!is.na(PGN_R2_comp_transition)) %>% 
  dplyr::mutate(percent= n/ sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=PGN_cluster_ident, y=percent, fill=PGN_R2_comp_transition), position = 'dodge')+
  labs(fill = "transition")+
  scale_fill_manual(values=c('A -> A' = '#9E31DA', 'A -> B' ='#825515','B -> A' =   '#c88ee8', 'B -> B' = '#F49D1F'))

#DN combined
esc_dn <- melting_scores %>% 
  dplyr::filter(meltingScore_DN_R1 > 5 | meltingScore_DN_R2 > 5) %>% # meltingScore_DN_R2 > 5
  dplyr::arrange(ESC_compartment, DN_R1_compartment, DN_R2_compartment)
esc_dn <- rowid_to_column(esc_dn)

col_runif_melting_vta = colorRamp2(c(0, 70), c("white", "#259A37"))
plot_data1 <- as.matrix(esc_dn[,c('ESC_compartment', 'DN_R1_compartment', 'DN_R2_compartment') ])
plot_data2 <- as.matrix(esc_dn[,c('meltingScore_DN_R1', 'meltingScore_DN_R2')])
plot_data3 <- log10(as.matrix(esc_dn[,c('rnaRPMM_ESC','rnaRPMM_DN') ]))
rownames(plot_data3) <- esc_dn$gene_symbol
plot_data4 <-  as.matrix(esc_dn[,c('percent_genebody_lad', 'percent_genebody_nad')]) 
plot_data4[is.na(plot_data4)] <- 0
plot_data5 <- as.matrix(esc_dn[,c('atacRPMM_ESC', 'atacRPMM_DN') ])

rownames(plot_data2) <- esc_dn$gene_symbol
compartments_h <- Heatmap(plot_data1, col=c('#9E31DA', '#929292',  '#F49D1F'), heatmap_legend_param=list(title='Compartment')) #, width = unit(2, "cm")
meltingscore <- Heatmap(plot_data2, col=col_runif_melting_vta, heatmap_legend_param=list(title='Melting score'),cluster_columns = FALSE, cluster_rows=FALSE, row_names_gp = gpar(fontsize = 7))
rpkm_h <- Heatmap(plot_data3, col=col_runif_rna, heatmap_legend_param=list(title='log10(lsRRPM)'), name='rpmm', cluster_columns = FALSE, row_km = 4, row_km_repeats = 1000, cluster_row_slices = FALSE, cluster_rows = FALSE) #cluster_rows = TRUE, row_km = 4, row_km_repeats = 100
loc_pos <- Heatmap(plot_data4, col=col_black_white, cluster_columns=FALSE, cluster_rows=FALSE, heatmap_legend_param=list(title='periphery assoc [%]'))
atac_h <- Heatmap(log10(plot_data5), col=col_runif_atac, heatmap_legend_param=list(title='log10(lsARPM)'), cluster_columns = FALSE, cluster_rows = FALSE)
ht_list <-  rpkm_h + atac_h + compartments_h   + loc_pos + meltingscore #+ rowAnnotation(log10_length = anno_points(log10(esc_dn$width), pch = 16, size = unit(1, "mm"), 
                                                                        #                                                     axis_param = list(at = c(5, 5.5, 6, 6.4)),
                                                                        #                                                     width = unit(2, "cm")))
draw(ht_list, ht_gap=unit(0, 'cm'),merge_legends = TRUE, heatmap_legend_side = "bottom")
dev.off()

vta_toplot <- melting_scores %>% 
  drop_na(DN_cluster_ident) %>% 
  dplyr::mutate(DN_cluster_ident = as.character(DN_cluster_ident))

#log2FC of RNA per cluster
vta_toplot  %>% 
  ggplot()+
  geom_violin(aes(x=DN_cluster_ident, y=rnaRPMM_l2FC_DNs))+
  geom_boxplot(aes(x=DN_cluster_ident, y=rnaRPMM_l2FC_DNs), outlier.shape = NA, width=0.2)+
  geom_jitter(aes(x=DN_cluster_ident, y=rnaRPMM_l2FC_DNs), size=0.5, width=0.2)+ 
  geom_hline(yintercept = 0, linetype='dotted')
#log2FC of ATAC per cluster
vta_toplot %>% 
  ggplot()+
  geom_violin(aes(x=DN_cluster_ident, y=atacRPMM_l2FC_DNs))+
  geom_boxplot(aes(x=DN_cluster_ident, y=atacRPMM_l2FC_DNs), outlier.shape = NA, width=0.2)+
  geom_jitter(aes(x=DN_cluster_ident, y=atacRPMM_l2FC_DNs), size=0.5, width=0.2)+
  geom_hline(yintercept = 0, linetype='dotted')

#melting_score in R1 and R2 per cluster
vta_toplot %>% 
  ggplot()+
  geom_violin(aes(x=DN_cluster_ident, y=meltingScore_DN_R1))+
  geom_boxplot(aes(x=DN_cluster_ident, y=meltingScore_DN_R1), outlier.shape = NA, width=0.2)+
  geom_jitter(aes(x=DN_cluster_ident, y=meltingScore_DN_R1), size=0.5, width=0.2)
vta_toplot %>% 
  ggplot()+
  geom_violin(aes(x=DN_cluster_ident, y=meltingScore_DN_R2))+
  geom_boxplot(aes(x=DN_cluster_ident, y=meltingScore_DN_R2), outlier.shape = NA, width=0.2)+
  geom_jitter(aes(x=DN_cluster_ident, y=meltingScore_DN_R2), size=0.5, width=0.2)

#compartment transitions ESC -> DN
vta_toplot %>% 
  group_by(DN_cluster_ident, DN_R1_comp_transition, .drop=FALSE) %>% 
  summarise(n=n()) %>%
  dplyr::filter(!is.na(DN_R1_comp_transition)) %>% 
  dplyr::mutate(percent= n/ sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=DN_cluster_ident, y=percent, fill=DN_R1_comp_transition), position = 'dodge')+
  labs(fill = "transition")+
  scale_fill_manual(values=c('A -> A' = '#9E31DA', 'A -> B' ='#825515','B -> A' =   '#c88ee8', 'B -> B' = '#F49D1F'))
vta_toplot %>% 
  group_by(DN_cluster_ident, DN_R2_comp_transition) %>% 
  summarise(n=n()) %>%
  dplyr::filter(!is.na(DN_R2_comp_transition)) %>% 
  dplyr::mutate(percent= n/ sum(n)) %>% 
  ggplot()+
  geom_col(aes(x=DN_cluster_ident, y=percent, fill=DN_R2_comp_transition), position = 'dodge')+
  labs(fill = "transition")+
  scale_fill_manual(values=c('A -> A' = '#9E31DA', 'A -> B' ='#825515','B -> A' =   '#c88ee8', 'B -> B' = '#F49D1F'))

