rm(list = ls())


library(ggspatial)
library(ggplot2)
theme_set(theme_bw())
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(tidyverse)
library(reshape2)
library(scatterpie)
library(GenomicOriginsScripts)
library(ggtext)
library(hypoimg)


meta <- read.csv("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/metadata_gxp_ben_floridae_complete", sep = ";")

meta <- meta %>% mutate(spec = ifelse(label == "PL17_23nigpue", "tan", spec), 
                        label = gsub("PL17_23nigpue", "PL17_23tanpue", label),
                        spec = ifelse(label == "PL17_35puepue", "ind", spec),
                        label = gsub("PL17_35puepue", "PL17_35indpue", label)) %>%
                 filter(., !label %in% c("PL17_101maybel", "AG9RX_47pueboc", "PL17_98indbel", "PL17_79abepue")) %>%
                 summarise(spec, geo) %>%
                 group_by(geo, spec) %>%
                 summarize(count = n())

p <- ggplot(meta, aes(x=geo, y=spec, label=count, color=spec)) +
     geom_point(aes(size = count)) +
  geom_text(color = "darkgrey", nudge_x = 0.25) +
  scale_x_discrete(labels = c('Belize','Panama','USA (Florida)', 'Puerto Rico'),
                   position = "top",
                   expand = c(0.2,0)) +
  coord_fixed(ratio = 0.45) +
  scale_y_discrete(expand = c(0.05,0.01),
                   labels = c("nig" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_nigricans.l.cairo.png' width='40' /><br>*H. nigricans*",
                              "chl" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_chlorurus.l.cairo.png' width='40' /><br>*H. chlorurus*",
                              "abe" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_aberrans.l.cairo.png' width='40' /><br>*H. aberrans*",
                              "gut" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_guttavarius.l.cairo.png' width='40' /><br>*H. guttavarius*",
                              "gum" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_gumigutta.l.cairo.png' width='40' /><br>*H. gummigutta*",
                              "ran" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_randallorum.l.cairo.png' width='40' /><br>*H. randallorum*",
                              "gem" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_gemma.l.cairo.png' width='40' /><br>*H. gemma*",
                              "may" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_maya.l.cairo.png' width='40' /><br>*H. maya*",
                              "ind" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_indigo.l.cairo.png' width='40' /><br>*H. indigo*",
                              "pue" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_puella.l.cairo.png' width='40' /><br>*H. puella*",
                              "flo" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_floridae.l.cairo.png' width='40' /><br>*H. floridae*",
                              "tan" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_tan.l.cairo.png' width='40' /><br>*Tan hamlet*",
                              "uni" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_unicolor.l.cairo.png' width='40' /><br>*H. unicolor*")) +
  scale_color_manual(values=c("nig" = '#FF0033', "chl" = '#9900CC', "abe" = '#996600', "gut" = '#0000FF',
                              "gum" = '#FF00FF', "ran" = '#666699', "gem" = '#CC0000', "may" = '#FF9933',
                              "ind" = '#66CCFF', "pue" = '#FFCC00', "flo" = '#33FFCC', "tan" = '#333333',
                              "uni" = '#66CC00')) +
  theme(legend.position="none",
        legend.title=element_blank(),
        legend.box = "vertical",
        legend.text =  element_markdown(size = 10),
        panel.background = element_blank(),
        text = element_text(size=10),
        legend.key=element_blank(),
        axis.text.x = element_text(size = 9, face = "bold"),
        axis.text.y = element_markdown(size = 9, face = "italic"),
        panel.grid.major = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank())

p

hypo_save(filename = "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/T1.png",
          type = "cairo",
          plot = p,
          width = 5.4,
          height = 6.5)

hypo_save(filename = "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/T1.pdf",
          plot = p,
          width = 5.4,
          height = 6.5)
