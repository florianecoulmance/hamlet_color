# by: Floriane Coulmance: 20/06/2022
# usage:
# Rscript sample_map.R <gxp_metadata>
# -------------------------------------------------------------------------------------------------------------------
# gxp_metadata in : "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/metadata_gxp_ben_floridae_complete"
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
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


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6]
print(args)

gxp_metadata <- as.character(args[1]) # Path to phenotype PCA files folder

# gxp_metadata <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/metadata_gxp_ben_floridae_complete"


# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


custom_geom_scatterpie_legend <- function (radius, x, y, n = 5, labeller, textsize=1.9) {
  
  # This is a custom function to change the size of piechart labels specifically
  # Taken from :
  # https://community.rstudio.com/t/pie-chart-world-map-for-genetics/97192/2
  
  if (length(radius) > n) {
    radius <- unique(sapply(seq(min(radius), max(radius), length.out = n), scatterpie:::round_digit))
    }
  label <- FALSE
  
    
  if (!missing(labeller)) {
    if (!inherits(labeller, "function")) {
      stop("labeller should be a function for converting radius")
      }
    label <- TRUE
    }
    
  dd <- data.frame(r = radius, start = 0, end = 2 * pi, x = x, y = y + radius - max(radius), maxr = max(radius))
  
  if (label) {
    dd$label <- labeller(dd$r)
    }
  else {
    dd$label <- dd$r
    }
  
  
  list(ggforce:::geom_arc_bar(aes_(x0 = ~x, y0 = ~y, r0 = ~r, r = ~r, start = ~start, end = ~end),
                              data = dd,
                              inherit.aes = FALSE), 
                 geom_segment(aes_(x = ~x, xend = ~x + maxr * 1.5, y = ~y + r, yend = ~y + r),
                              data = dd,
                              inherit.aes = FALSE), 
                 geom_text(aes_(x = ~x + maxr * 1.6, y = ~y + r, label = ~label),
                           data = dd,
                           hjust = "left",
                           inherit.aes = FALSE,
                           size=textsize))
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Read metadata table for the GxP dataset
sample_list <- read.csv(gxp_metadata, sep=";") %>%
               mutate(spec = ifelse(label == "PL17_23nigpue", "tan", spec),
                      label = gsub("PL17_23nigpue", "PL17_23tanpue", label)) %>%
               filter(., !label %in% c("PL17_101maybel", "AG9RX_47pueboc", "PL17_98indbel", "PL17_79abepue"))

# Retrieve GPS coordinations and select columns in metadata
sample_list <- sample_list %>%
               mutate(coord_N.x = case_when(endsWith(geo, "bel") ~ 16.765278,
                                            endsWith(geo, "boc") ~ 9.332778,
                                            endsWith(geo, "hon") ~ 16.030000,
                                            endsWith(geo, "pue") ~ 17.952820,
                                            endsWith(geo, "gun") ~ 9.290722,
                                            endsWith(geo, "qui") ~ 20.978806,
                                            endsWith(geo, "san") ~ 12.501617,
                                            endsWith(geo, "arc") ~ 20.203940,
                                            endsWith(geo, "liz") ~ 19.155547,
                                            endsWith(geo, "tam") ~ 21.475878,
                                            endsWith(geo, "are") ~ 22.115278,
                                            endsWith(geo, "ala") ~ 22.399490,
                                            endsWith(geo, "flo") ~ 24.752580,
                                            endsWith(geo, "bar") ~ 13.223694,
                                            endsWith(geo, "hai") ~ 19.677000),
                      coord_W.x = case_when(endsWith(geo, "bel") ~ -88.14417,
                                            endsWith(geo, "boc") ~ -82.25472,
                                            endsWith(geo, "hon") ~ -83.32861,
                                            endsWith(geo, "pue") ~ -67.05643,
                                            endsWith(geo, "gun") ~ -78.14039,
                                            endsWith(geo, "qui") ~ -86.800111,
                                            endsWith(geo, "san") ~ -81.731367,
                                            endsWith(geo, "arc") ~ -91.970770,
                                            endsWith(geo, "liz") ~ -95.863903,
                                            endsWith(geo, "tam") ~ -97.227000,
                                            endsWith(geo, "are") ~ -91.398333,
                                            endsWith(geo, "ala") ~ -89.656060,
                                            endsWith(geo, "flo") ~ -80.76065,
                                            endsWith(geo, "bar") ~ -59.64506,
                                            endsWith(geo, "hai") ~ -71.843000)) %>%
               summarise(label, spec, geo, coord_N.x, coord_W.x) %>%
               setNames(., nm = c("SampleID", "spec", "geo", "coord_N.x", "coord_W.x"))

# Create a species count for each location and add a radius column for each location to use in the plot
pie_dat <- dcast(sample_list,geo~spec, fill = 0, value.var="spec", fun.aggregate = length) %>%
           mutate(coord_N.x = case_when(endsWith(geo, "bel") ~ 16.765278,
                                        endsWith(geo, "boc") ~ 9.332778,
                                        endsWith(geo, "hon") ~ 16.030000,
                                        endsWith(geo, "pue") ~ 17.952820,
                                        endsWith(geo, "gun") ~ 9.290722,
                                        endsWith(geo, "qui") ~ 20.978806,
                                        endsWith(geo, "san") ~ 12.501617,
                                        endsWith(geo, "arc") ~ 20.203940,
                                        endsWith(geo, "liz") ~ 19.155547,
                                        endsWith(geo, "tam") ~ 21.475878,
                                        endsWith(geo, "are") ~ 22.115278,
                                        endsWith(geo, "ala") ~ 22.399490,
                                        endsWith(geo, "flo") ~ 24.752580,
                                        endsWith(geo, "bar") ~ 13.223694,
                                        endsWith(geo, "hai") ~ 19.677000),
                  coord_W.x = case_when(endsWith(geo, "bel") ~ -88.14417,
                                        endsWith(geo, "boc") ~ -82.25472,
                                        endsWith(geo, "hon") ~ -83.32861,
                                        endsWith(geo, "pue") ~ -67.05643,
                                        endsWith(geo, "gun") ~ -78.14039,
                                        endsWith(geo, "qui") ~ -86.800111,
                                        endsWith(geo, "san") ~ -81.731367,
                                        endsWith(geo, "arc") ~ -91.970770,
                                        endsWith(geo, "liz") ~ -95.863903,
                                        endsWith(geo, "tam") ~ -97.227000,
                                        endsWith(geo, "are") ~ -91.398333,
                                        endsWith(geo, "ala") ~ -89.656060,
                                        endsWith(geo, "flo") ~ -80.76065,
                                        endsWith(geo, "bar") ~ -59.64506,
                                        endsWith(geo, "hai") ~ -71.843000),
                  sum = as.numeric(rowSums(.[,2:14])),
                  radius = sum/15)

# Create the world map data to use in plot
world <- map_data('world')

# Create the plot
p <- ggplot(data = world, aes(long, lat)) +
     geom_map(map=world, aes(map_id=region), fill="peachpuff1", color = "antiquewhite") +
     coord_quickmap() +
     geom_scatterpie(aes(x=coord_W.x, y=coord_N.x, group=geo, r=radius),
                     data=pie_dat, cols=colnames(pie_dat)[2:14], color=NA, alpha=.8) +
     custom_geom_scatterpie_legend(pie_dat$radius, x=-64, y=27, n = 4, labeller = function(x) round(x*15, digits = 0)) +
     annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", fontface = "italic", color = "grey22", size = 4) +
     annotate(geom = "text", x = -77, y = 15, label = "Carribbean Sea", fontface = "italic", color = "grey22", size = 4) +
     annotate(geom = "text", x = -73, y = 30, label = "Atlantic", fontface = "italic", color = "grey22", size = 4) +
     annotate(geom = "text", x = -64, y = 30, label = "Sample size", fontface = "bold", color = "black", size = 3) +
     # annotation_scale(location = "bl", width_hint = 0.5) +
     # annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
     coord_sf(xlim = c(-100.15, -55.12), ylim = c(7.00, 31.97), expand = FALSE) +
     # xlab("Longitude") + ylab("Latitude") +
     # ggtitle("Map of the Gulf of Mexico and the Caribbean Sea") +
     scale_fill_manual(values=c("abe" = '#996600',
                                "chl" = '#9900CC',
                                "flo" = '#33FFCC',
                                "gem" = '#CC0000',
                                "gum" = '#FF00FF',
                                "gut" = '#0000FF',
                                "ind" = '#66CCFF',
                                "may" = '#FF9933',
                                "nig" = '#FF0033',
                                "pue" = '#FFCC00',
                                "ran" = '#666699',
                                "uni" = '#66CC00',
                                "tan" = '#333333'),
                       labels = c("abe" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_aberrans.l.cairo.png' width='40' /><br>*H. aberrans*",
                                  "chl" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_chlorurus.l.cairo.png' width='40' /><br>*H. chlorurus*",
                                  "flo" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_floridae.l.cairo.png' width='40' /><br>*H. floridae*",
                                  "gem" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_gemma.l.cairo.png' width='40' /><br>*H. gemma*",
                                  "gum" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_gumigutta.l.cairo.png' width='40' /><br>*H. gummigutta*",
                                  "gut" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_guttavarius.l.cairo.png' width='40' /><br>*H. guttavarius*",
                                  "ind" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_indigo.l.cairo.png' width='40' /><br>*H. indigo*",
                                  "may" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_maya.l.cairo.png' width='40' /><br>*H. maya*",
                                  "nig" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_nigricans.l.cairo.png' width='40' /><br>*H. nigricans*",
                                  "pue" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_puella.l.cairo.png' width='40' /><br>*H. puella*",
                                  "ran" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_randallorum.l.cairo.png' width='40' /><br>*H. randallorum*",
                                  "uni" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_unicolor.l.cairo.png' width='40' /><br>*H. unicolor*",
                                  "tan" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_tan.l.cairo.png' width='40' /><br>*Tan hamlet*")) +
     guides(fill = guide_legend(nrow = 7)) +
     theme(legend.title=element_blank(),
           panel.background = element_rect(fill = "lightblue2"),
           panel.grid.major = element_blank(),
           legend.text =  element_markdown(size = 10),
           legend.key=element_blank(),
           legend.box = "vertical",
           axis.title = element_blank(),
           axis.ticks = element_blank(),
           axis.text = element_blank(),
           panel.border = element_blank(),
           text = element_text(size = 3))

p

# Save the plot
hypo_save(filename = "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/Hypoplectrus_gxp_map.png",
          type = "cairo",
          plot = p,
          width = 11,
          height = 5)

hypo_save(filename = "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/Hypoplectrus_gxp_map.pdf",
          plot = p,
          width = 11,
          height = 5)










# world <- ne_countries(scale = "medium", returnclass = "sf")
# ggplot(data = world) + geom_sf(fill= "antiquewhite") + 
#   geom_point(data = test, mapping = aes(x = coord_W, y = coord_N), color = "red") +
#   annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", fontface = "italic", color = "grey22", size = 3) +
#   annotation_scale(location = "bl", width_hint = 0.5) +
#   annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
#   coord_sf(xlim = c(-102.15, -60.12), ylim = c(5.00, 33.97), expand = FALSE) +
#   xlab("Longitude") + ylab("Latitude") +
#   ggtitle("Map of the Gulf of Mexico and the Caribbean Sea") +
#   theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))
