rm(list = ls())


library(ggplot2)
library(tidyverse)
library(hypogen)
library(dplyr)
library(plyr)

setwd('/Users/fco/Desktop/BREMEN_OP/chapter1_2/GxP')

#df_gxp_rand <- read.table(file = 'H.randallorum.lm.50k.5k.txt.gz', sep = '\t', header = TRUE)
# df_gxp_bars <- read.table(file = 'bars_body.lm.50k.5k.txt.gz', sep = '\t', header = TRUE)
# 

files <- list.files(pattern = "k.txt.gz")
traits <- list("bars_body", "bars_head", "peduncle", "snout", "H.unicolor", "H.aberrans", "H.nigricans", "H.maya", "H.gemma", "H.puella", "H.indigo","H.gummiguta", "Tan.hamlet", "H.guttavarius", "H.floridae", "H.chlorurus", "H.randallorum")

for(i in 1:length(traits)) {
  
  pdf(paste0(traits[i],"_plots"))
  
  par(mfrow=c(2,2))

  for(j in 1:length(files)){
    
    
    if(substring(files[j], 1, 1)=="H" | substring(files[j], 1, 1)=="T"){
      print("starts with H")
      traitfile <- strsplit(files[j], ".lm")[[1]][1]
      print(traitfile)
    } else{
      print("Dont start with H")
      traitfile <- regmatches(x=files[j],gregexpr('[.]',files[j]),invert=TRUE)[[1]][1]
      print(traitfile)
      
    }
    
    if(traitfile==traits[i]) {
      print("trait and file match")
      data <- read_tsv(files[j]) %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
      name <- tools::file_path_sans_ext(files[j])
      name <- tools::file_path_sans_ext(name)
      print(name)
      
      
      p <- ggplot(data = data, aes(x = GPOS, y = AVG_p_wald))+
        geom_hypo_LG()+
        geom_point(size = .2)+
        scale_fill_hypo_LG_bg()+
        scale_x_hypo_LG()+
        theme_hypo() + ggtitle(name)
      
      print(p)
       
    }
  
  }
  print("Close the plot")
  dev.off()
  
  
}


