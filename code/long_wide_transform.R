#library("reshape2")
library("dplyr")
library("readr")
library("tidyr")

options(scipen=999999)


#library(doParallel)
#library(foreach)
#registerDoParallel(cores=30)

res=50000


pth="CN"
rstudioapi::jobRunScript("reshape.R", importEnv = TRUE)
pth="NPC"
rstudioapi::jobRunScript("reshape.R", importEnv = TRUE)
pth="ES"
rstudioapi::jobRunScript("reshape.R", importEnv = TRUE)



      #args = commandArgs(trailingOnly = TRUE)
      #infile=args[1]
      #chrom=args[2]
      #res=as.integer(args[3])


for (norm in c("raw", "KR")){
  for (norm2 in c("NONE", "KR")){
    for (chm in paste0("chr", c(1:19))) {
      
      infile=sprintf("/data/pombo/christoph/datasets/bonev20201002/%s_%s/observed_%s_%s_%s.txt",  pth, norm, norm2, chm, chm)
      if(file.exists(infile)==F){next()}
      
      chrom=chm

      print(paste(infile, chrom, res))
      
      longdata = read_delim(infile, col_names = c("start.x", "start.y", "value"), delim = "\t")
      
      # is symmetric?
      tmp = longdata[17,]
      tmpstart.x=tmp$start.y
      tmpstart.y=tmp$start.x
      
      if (nrow(longdata %>% filter(start.x==tmpstart.x, start.y==tmpstart.y))==0){
        #missing bins
        p=seq(0,max(longdata$start.x,longdata$start.y),res)
        
        longdata=longdata %>% 
          #make data symmetrically
          bind_rows(data.frame("start.x"=longdata$start.y, "start.y"=longdata$start.x, "value"=longdata$value)) %>%
          #add missing bins at diagonal
          bind_rows(data.frame("start.x"=p, "start.y"=p, "value"=NA))
      }else{
        print("please check input")}
      
      
      #x=c(NA,NA,NA)
      aggfun = function(x){
        if(sum(!is.na(x))==0){return(NA)
        }else{return(max(x, na.rm = T))} 
      }
      
      
      widedata = longdata %>% 
        group_by(start.x, start.y) %>% summarize(value = aggfun(value)) %>% 
        mutate(value = ifelse(is.infinite(value), NA, value)) %>%
        ungroup() %>%
        arrange(start.y, start.x) %>% 
        mutate(name.x=sprintf("%s:%s-%s", chrom, start.x, start.x+res),
               name.y=sprintf("%s:%s-%s", chrom, start.y, start.y+res)) %>%
        pivot_wider(id_cols = c(name.x, start.x), names_from=name.y, values_from=value, values_fn=max) %>%
        #widedata2 = widedata %>% 
        arrange(start.x) 
      
      
      colnames(widedata)[1]="HiC"
      
      
      write.table(widedata %>% dplyr::select(-start.x), paste0(infile, ".wideR"),col.names = T, row.names = F, quote=F, sep="\t")
      
    }
  }
}


