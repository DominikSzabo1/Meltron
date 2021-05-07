## plot the empirical cumulative density functions (ECDF) of IS values of long genes.

library(tidyverse)

##load data generated with smelte.R script
#list that contains 
load('../data/IS_values_perGene.RData')

#####PLOT empirical cumulative density function (ECDF) plots of IS value distributions per gene:
####ECDF plot
plot_smooth_ecdf <- function(gene, celltype1=c('OLG', 'DN_R1', 'PGN_R1', 'all'), celltype2=F, distance_label=T){
  if(celltype1=='OLG'){
    scores_gene <- scores[[gene]]
    scores_gene_stem <- unlist(scores_gene[['ESC']]) # raw insulation values
    scores_gene_stem_ecdf <-ecdf(scores_gene_stem) # ecdf values
    scores_gene_comp <- unlist(scores_gene[[celltype1]]) 
    scores_gene_comp_ecdf <-ecdf(scores_gene_comp) # ecdf values
    stat <- ks.test(scores_gene_stem, scores_gene_comp, alternative = 'less')
    stat_D <- formatC(round(stat$statistic, 2), format='f', digits=2) # maximum distance between curves
    minMaxR1 <- seq(min(scores_gene_stem, scores_gene_comp), max(scores_gene_stem, scores_gene_comp), length.out=length(scores_gene_stem)) #
    x0R1 <- minMaxR1[which( abs(scores_gene_stem_ecdf(minMaxR1) - scores_gene_comp_ecdf(minMaxR1)) == max(abs(scores_gene_stem_ecdf(minMaxR1) - scores_gene_comp_ecdf(minMaxR1))) )] 
    y0R1 <- scores_gene_stem_ecdf(x0R1)
    y1R1 <- scores_gene_comp_ecdf(x0R1)
    ##smoothing for plotting puropses:
    dens_stem = density(scores_gene_stem, adjust=0.3)
    dens_compR1 = density(scores_gene_comp, adjust=0.3)
    dens = data.frame(stem_x=dens_stem$x, stem_y=dens_stem$y,
                      compR1_x=dens_compR1$x, compR1_y=dens_compR1$y)
    
    #if distance label needs to be added, execute this part
    if(distance_label==TRUE){
      ecdf_smooth <- ggplot() + 
        geom_line(data=dens, aes(x=stem_x, y=cumsum(stem_y)/sum(stem_y)), size=1.7, colour='#DE9132') +
        geom_line(data=dens, aes(x=compR1_x, y=cumsum(compR1_y)/sum(compR1_y)), size=1.7, colour='#800080') +
        geom_segment(aes(x = x0R1[1], y = y0R1[1], xend = x0R1[1], yend = y1R1[1]),
                     linetype = "dotted", color = "grey30", size=1.5) +
        annotate(geom='text', x =  x0R1[1] + 0.08, y = (y1R1[1] - y0R1[1])/2 + y0R1[1] , label=stat_D, size=6)
      return(ecdf_smooth)
    }
    #if distance label not needed, execute this part
    else{
      ecdf_smooth <- ggplot() + 
        geom_line(data=dens, aes(x=stem_x, y=cumsum(stem_y)/sum(stem_y)), size=1.7, colour='#DE9132') +
        geom_line(data=dens, aes(x=compR1_x, y=cumsum(compR1_y)/sum(compR1_y)), size=1.7, colour='#800080') +
        geom_segment(aes(x = x0R1[1], y = y0R1[1], xend = x0R1[1], yend = y1R1[1]),
                     linetype = "dotted", color = "grey30", size=1.5) 
      
      return(ecdf_smooth)
    } 
  }
  
  #If you want to plot DN or PGN, you need to specify celltype2
  if(celltype1 %in% c('DN_R1', 'PGN_R1')){ # 2 replicates available
    scores_gene <- scores[[gene]]
    scores_gene_stem <- unlist(scores_gene[['ESC']]) # raw insulation values
    scores_gene_stem_ecdf <-ecdf(scores_gene_stem) # ecdf values
    scores_gene_compR1 <- unlist(scores_gene[[celltype1]]) 
    scores_gene_compR1_ecdf <-ecdf(scores_gene_compR1) # ecdf values
    statR1 <- ks.test(scores_gene_stem, scores_gene_compR1, alternative = 'less')
    stat_DR1 <- formatC(round(statR1$statistic, 2), format='f', digits=2) # maximum distance between curves
    minMaxR1 <- seq(min(scores_gene_stem, scores_gene_compR1), max(scores_gene_stem, scores_gene_compR1), length.out=length(scores_gene_stem)) 
    x0R1 <- minMaxR1[which( abs(scores_gene_stem_ecdf(minMaxR1) - scores_gene_compR1_ecdf(minMaxR1)) == max(abs(scores_gene_stem_ecdf(minMaxR1) - scores_gene_compR1_ecdf(minMaxR1))) )] 
    y0R1 <- scores_gene_stem_ecdf(x0R1)
    y1R1 <- scores_gene_compR1_ecdf(x0R1)
    
    scores_gene_compR2 <- unlist(scores_gene[[celltype2]]) 
    scores_gene_compR2_ecdf <-ecdf(scores_gene_compR2) # ecdf values
    statR2 <- ks.test(scores_gene_stem, scores_gene_compR2, alternative = 'less')
    stat_DR2 <- formatC(round(statR2$statistic, 2), format='f', digits=2) # maximum distance between curves
    minMaxR2 <- seq(min(scores_gene_stem, scores_gene_compR2), max(scores_gene_stem, scores_gene_compR2), length.out=length(scores_gene_stem)) 
    x0R2 <- minMaxR2[which( abs(scores_gene_stem_ecdf(minMaxR2) - scores_gene_compR2_ecdf(minMaxR2)) == max(abs(scores_gene_stem_ecdf(minMaxR2) - scores_gene_compR2_ecdf(minMaxR2))) )] 
    y0R2 <- scores_gene_stem_ecdf(x0R2)
    y1R2 <- scores_gene_compR2_ecdf(x0R2)
    ##smoothing:
    dens_stem = density(scores_gene_stem, adjust=0.3)
    dens_compR1 = density(scores_gene_compR1, adjust=0.3)
    dens_compR2 = density(scores_gene_compR2, adjust=0.3)
    dens = data.frame(stem_x=dens_stem$x, stem_y=dens_stem$y,
                      compR1_x=dens_compR1$x, compR1_y=dens_compR1$y,
                      compR2_x=dens_compR2$x, compR2_y=dens_compR2$y)
    if(distance_label==TRUE){
      if(celltype1 == 'DN_R1'){
        color1='#259A37'
        color2='#1d822c'
      }
      if(celltype1=='PGN_R1'){
        color1='#6367DC'
        color2='#4c50ad'
      }
      ecdf_smooth <- ggplot() + 
        geom_line(data=dens, aes(x=stem_x, y=cumsum(stem_y)/sum(stem_y)), size=1.7, colour='#DE9132') +
        geom_line(data=dens, aes(x=compR1_x, y=cumsum(compR1_y)/sum(compR1_y)), size=1.7, colour=color1) +
        geom_line(data=dens, aes(x=compR2_x, y=cumsum(compR2_y)/sum(compR2_y)), size=1.7, colour=color2) +
        geom_segment(aes(x = x0R1[1], y = y0R1[1], xend = x0R1[1], yend = y1R1[1]),
                     linetype = "dotted", color = "grey30", size=1.5) +
        annotate(geom='text', x =  x0R1[1] + 0.08, y = (y1R1[1] - y0R1[1])/2 + y0R1[1] , label=stat_DR1, size=6) +
        geom_segment(aes(x = x0R2[1], y = y0R2[1], xend = x0R2[1], yend = y1R2[1]),
                     linetype = "dotted", color = "grey30", size=1.5) +
        annotate(geom='text', x =  x0R2[1] + 0.08, y = (y1R2[1] - y0R2[1])/2 + y0R2[1] + 0.15, label=stat_DR2, size=6)
      return(ecdf_smooth)
    }
    else{
      if(celltype1 == 'DN_R1'){
        color1='#259A37'
        color2='#1d822c'
      }
      if(celltype1=='PGN_R1'){
        color1='#6367DC'
        color2='#4c50ad'
      }
      ecdf_smooth <- ggplot() + 
        geom_line(data=dens, aes(x=stem_x, y=cumsum(stem_y)/sum(stem_y)), size=1.7, colour='#DE9132') +
        geom_line(data=dens, aes(x=compR1_x, y=cumsum(compR1_y)/sum(compR1_y)), size=1.7, colour=color1) +
        geom_line(data=dens, aes(x=compR2_x, y=cumsum(compR2_y)/sum(compR2_y)), size=1.7, colour=color2) +
        geom_segment(aes(x = x0R1[1], y = y0R1[1], xend = x0R1[1], yend = y1R1[1]),
                     linetype = "dotted", color = "grey30", size=1.5) +
        geom_segment(aes(x = x0R2[1], y = y0R2[1], xend = x0R2[1], yend = y1R2[1]),
                     linetype = "dotted", color = "grey30", size=1.5) 
      
      return(ecdf_smooth)
    }
  }
  if(celltype1=='all'){
    scores_gene <- scores[[gene]]
    scores_gene_stem <- unlist(scores_gene[['ESC']]) # raw insulation values
    scores_gene_oligo <- unlist(scores_gene[['OLG']])
    scores_gene_vtaR1 <- unlist(scores_gene[['DN_R1']])
    scores_gene_vtaR2 <- unlist(scores_gene[['DN_R2']])
    scores_gene_ca1R1 <- unlist(scores_gene[['PGN_R1']])
    scores_gene_ca1R2 <- unlist(scores_gene[['PGN_R2']])
    ##smoothing:
    dens_stem = density(scores_gene_stem, adjust=0.3)
    dens_oligo = density(scores_gene_oligo, adjust=0.3)
    dens_vtaR1 = density(scores_gene_vtaR1, adjust=0.3)
    dens_vtaR2 = density(scores_gene_vtaR2, adjust=0.3)
    dens_ca1R1 = density(scores_gene_ca1R1, adjust=0.3)
    dens_ca1R2 = density(scores_gene_ca1R2, adjust=0.3)
    dens = data.frame(stem_x=dens_stem$x, stem_y=dens_stem$y,
                      oligo_x=dens_oligo$x, oligo_y=dens_oligo$y,
                      vtaR1_x=dens_vtaR1$x, vtaR1_y=dens_vtaR1$y,
                      vtaR2_x=dens_vtaR2$x, vtaR2_y=dens_vtaR2$y,
                      ca1R1_x=dens_ca1R1$x, ca1R1_y=dens_ca1R1$y,
                      ca1R2_x=dens_ca1R2$x, ca1R2_y=dens_ca1R2$y)
    #plotting:
    ecdf_smooth <- ggplot() + 
      geom_line(data=dens, aes(x=stem_x, y=cumsum(stem_y)/sum(stem_y)), size=1.7, colour='#DE9132') +
      geom_line(data=dens, aes(x=oligo_x, y=cumsum(oligo_y)/sum(oligo_y)), size=1.7, colour='#800080') +
      geom_line(data=dens, aes(x=vtaR1_x, y=cumsum(vtaR1_y)/sum(vtaR1_y)), size=1.7, colour='#259A37') +
      geom_line(data=dens, aes(x=vtaR2_x, y=cumsum(vtaR2_y)/sum(vtaR2_y)), size=1.7, colour='#1d822c') +
      geom_line(data=dens, aes(x=ca1R1_x, y=cumsum(ca1R1_y)/sum(ca1R1_y)), size=1.7, colour='#6367DC') +
      geom_line(data=dens, aes(x=ca1R2_x, y=cumsum(ca1R2_y)/sum(ca1R2_y)), size=1.7, colour='#4c50ad') 
    return(ecdf_smooth)
  }
}  

###if you want to plot DN or PGN: R1 needs to be celltype1 and you need to specify corresponding R2 as celltype2
#KS test is one-sided: negative distances between curves will not be calculated
plot_smooth_ecdf(gene = 'Rbfox1', celltype1 = 'PGN_R1', celltype2 = 'PGN_R2', distance_label=TRUE)
plot_smooth_ecdf(gene = 'Rbfox1', celltype1 = 'DN_R1', celltype2 = 'DN_R2', distance_label=TRUE)
plot_smooth_ecdf(gene = 'Rbfox1', celltype1 = 'OLG', distance_label=TRUE)
plot_smooth_ecdf(gene = 'Rbfox1', celltype1 = 'all', distance_label=TRUE)
