## plot expression and chromatin accessability over the melting score per gene

library(tidyverse)
library(ggpubr)

#read in files
scores_red <- read_tsv('../data/melting_scores.tsv.gz')
rna_atac <- read_tsv('../data/rna_atac.tsv.gz')
scores_expr_mod <- left_join(scores_red, rna_atac)

#RNA over melting score
plot_rna <- function(celltype=c('OLG', 'DN_R1', 'DN_R2', 'PGN_R1', 'PGN_R2')){
  if(celltype == 'OLG'){
    x='meltingScore_OLG'
    y='rnaRPMM_OLGs'
    top='melting_in_OLG'
    color="#800080" #bd06bd
  }  
  if(celltype == 'DN_R1'){
    x='meltingScore_DN_R1'
    y='rnaRPMM_DN'
    top='melting_in_DN_R1'
    color="#259A37"
  }
  if(celltype == 'DN_R2'){
    x='meltingScore_DN_R2'
    y='rnaRPMM_DN'
    top='melting_in_DN_R2'
    color="#259A37" #
  }
  if(celltype == 'PGN_R1'){
    x='meltingScore_PGN_R1'
    y='rnaRPMM_PGN'
    top='melting_in_PGN_R1'
    color="#6367DC"
  }
  if(celltype == 'PGN_R2'){
    x='meltingScore_PGN_R2'
    y='rnaRPMM_PGN'
    top='melting_in_PGN_R2'
    color="#6367DC"
  }
  medi <- scores_expr_mod %>% group_by(filt =!!sym(top) == TRUE) %>% 
    summarise(median = median(log10(!!sym(y)), na.rm = T), n=n()) %>% 
    dplyr::mutate(col = case_when(filt == T ~ color, T ~ 'grey50'))
  test <- wilcox.test(log10(pull((scores_expr_mod %>% dplyr::filter(!!sym(top)==TRUE))[sym(y)])), 
                      log10(pull((scores_expr_mod %>% dplyr::filter(!!sym(top)==FALSE))[sym(y)])))
  pval <- format(test$p.value,scientific = TRUE, digits = 3)
  atac_col <- colorRampPalette(c('#000000', color))
  
  p1 <- ggplot() +     
    geom_point(data=scores_expr_mod %>% dplyr::filter(!!sym(top) == FALSE), aes(x=!!sym(x), y=log10(!!sym(y))), color='grey40', size=3.5, alpha=0.3)  +
    geom_point(data=scores_expr_mod %>% dplyr::filter(!!sym(top) == TRUE), aes(x=!!sym(x), y=log10(!!sym(y))), alpha=0.8, color=color, size=3.5)  +
    ylab('atac RPMM')+
    scale_y_continuous(breaks=c(1,2,3))+
    theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.title.y = element_text(size=14),axis.text = element_text(size=14))

    p2 <- ggplot() + 
    geom_density(data=scores_expr_mod, aes(x=log10(!!sym(y)), fill=!!sym(top)), alpha=0.5) +
    geom_vline(data=medi, aes(xintercept=median, color=col), linetype='twodash', size=2)+
    scale_fill_manual(values=c('grey40',color)) + 
    scale_color_manual(values=c(color, 'grey40')) +
    scale_y_continuous(breaks=c(0,1,2))+ #,limits = c(0,1) for oligo 1.0 tickmark
    theme(legend.position = 'none', axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text = element_text(size=14)) + 
    coord_flip()

    p3 <- ggplot()+
    geom_text(aes(y= 0.1, x =1, label=pval), size=8, angle=-90) +
    theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_blank())  
  
  together <- ggarrange(p1,p2, ncol=2, widths = c(1,0.3), align = 'hv') 
  print(pval)
  return(together)
}
plot_rna(celltype='PGN_R2') # options: OLG, DN_R1, DN_R2, PGN_R1, PGN_R2


#ATAC over melting score
plot_atac <- function(celltype=c('OLG', 'DN_R1', 'DN_R2', 'PGN_R1', 'PGN_R2')){
  if(celltype == 'OLG'){
    x='meltingScore_OLG'
    y='atacRPMM_OLGs'
    top='melting_in_OLG'
    color="#800080" #bd06bd
  }  
  if(celltype == 'DN_R1'){
    x='meltingScore_DN_R1'
    y='atacRPMM_DN'
    top='melting_in_DN_R1'
    color="#259A37"
  }
  if(celltype == 'DN_R2'){
    x='meltingScore_DN_R2'
    y='atacRPMM_DN'
    top='melting_in_DN_R2'
    color="#259A37" #
  }
  if(celltype == 'PGN_R1'){
    x='meltingScore_PGN_R1'
    y='atacRPMM_PGN'
    top='melting_in_PGN_R1'
    color="#6367DC"
  }
  if(celltype == 'PGN_R2'){
    x='meltingScore_PGN_R2'
    y='atacRPMM_PGN'
    top='melting_in_PGN_R2'
    color="#6367DC"
  }
  medi <- scores_expr_mod %>% group_by(filt =!!sym(top) == TRUE) %>% 
    summarise(median = median(log10(!!sym(y)), na.rm = T), n=n()) %>% 
    dplyr::mutate(col = case_when(filt == T ~ color, T ~ 'grey50'))
  test <- wilcox.test(log10(pull((scores_expr_mod %>% dplyr::filter(!!sym(top)==TRUE))[sym(y)])), 
                      log10(pull((scores_expr_mod %>% dplyr::filter(!!sym(top)==FALSE))[sym(y)])))
  pval <- format(test$p.value,scientific = TRUE, digits = 3)
  atac_col <- colorRampPalette(c('#000000', color))
  
  p1 <- ggplot() +     
    geom_point(data=scores_expr_mod %>% dplyr::filter(!!sym(top) == FALSE), aes(x=!!sym(x), y=log10(!!sym(y))), color='grey40', size=3.5, alpha=0.3)  +
    geom_point(data=scores_expr_mod %>% dplyr::filter(!!sym(top) == TRUE), aes(x=!!sym(x), y=log10(!!sym(y))), alpha=0.8, color=color, size=3.5)  +
    ylab('atac RPMM')+
    scale_y_continuous(breaks=c(1,2,3))+
    theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.title.y = element_text(size=14),axis.text = element_text(size=14))
  
  p2 <- ggplot() + 
    geom_density(data=scores_expr_mod, aes(x=log10(!!sym(y)), fill=!!sym(top)), alpha=0.5) +
    geom_vline(data=medi, aes(xintercept=median, color=col), linetype='twodash', size=2)+
    scale_fill_manual(values=c('grey40',color)) + 
    scale_color_manual(values=c(color, 'grey40')) +
    scale_y_continuous(breaks=c(0,1,2))+ #,limits = c(0,1) for oligo 1.0 tickmark
    theme(legend.position = 'none', axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text = element_text(size=14)) + 
    coord_flip()

  p3 <- ggplot()+
    geom_text(aes(y= 0.1, x =1, label=pval), size=8, angle=-90) +
    theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), axis.line = element_blank())  
  
  together <- ggarrange(p1,p2, ncol=2, widths = c(1,0.3), align = 'hv') 
  print(pval)
  return(together)
}
plot_atac(celltype='OLG') # options: OLG, DN_R1, DN_R2, PGN_R1, PGN_R2
