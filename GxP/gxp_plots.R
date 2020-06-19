rm(list = ls())


library(ggplot2)
library(tidyverse)
library(hypogen)
library(dplyr)
library(plyr)
library(ggpubr)


setwd('/Users/fco/Desktop/BREMEN_OP/chapter1_2/GxP')

#df_gxp_rand <- read.table(file = 'H.randallorum.lm.50k.5k.txt.gz', sep = '\t', header = TRUE)
df_gxp_bars <- read.table(file = 'bbody.lm.50k.5k.txt.gz', sep = '\t', header = TRUE)
# 

df_gxp_bars <- read_tsv('bbody.lm.50k.5k.txt.gz') %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)




files <- list.files(pattern = "k.txt.gz")
traits <- list("bhead", "bbody", "pue", "nig", "combo_spec", "ind", "abe", "may",  "ran", "gem", "gut", "flo", "chl", "tan", "snout", "ped", "gum", "uni")

for(i in 1:length(traits)) {
  
  png(paste0(traits[i],"_plots.png"))
  #par(mfrow=c(2,2))
  h <- 0
  pltList <- list()
  

  for(j in 1:length(files)){

      traitfile <- regmatches(x=files[j],gregexpr('[.]',files[j]),invert=TRUE)[[1]][1]
      print(traitfile)

    
    if(traitfile==traits[i]) {
      h <- h + 1
      print(h)
      print("trait and file match")
      data <- read_tsv(files[j]) %>% left_join(hypo_chrom_start) %>% mutate(GPOS = MID_POS + GSTART)
      name <- tools::file_path_sans_ext(files[j])
      name <- tools::file_path_sans_ext(name)
      print(name)
      
      
      pltName <- paste( 'p', h, sep = '_' )
      print(pltName)
      pltList[[ pltName ]] <- ggplot(data = data, aes(x = GPOS, y = AVG_p_wald))+
        geom_hypo_LG()+
        geom_point(size = .2)+
        scale_fill_hypo_LG_bg()+
        scale_x_hypo_LG()+
        theme_hypo() + ggtitle(name)
      

    }
  
  }
  figure <- ggarrange(pltList[[1]],pltList[[2]],pltList[[3]] ,pltList[[4]] ,ncol = 2, nrow = 2)
  print(figure)
  print("Close the plot")
  dev.off()
  
  
}


