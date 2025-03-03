# by: Floriane Coulmance: 20/06/2022
# usage:
# Rscript sample_map.R <gxp_metadata> <coverage_table>
# -------------------------------------------------------------------------------------------------------------------
# gxp_metadata in : "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/metadata_gxp_ben_floridae_complete"
# coverage_table in : $BASE_DIR/outputs/coverage/coverage_table
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(ggspatial)
# theme_set(theme_bw())
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(tidyverse)
library(reshape2)
library(scatterpie)
library(GenomicOriginsScripts)
library(ggplot2)
library(ggtext)
library(hypoimg)
library(gridExtra)
library(dplyr)



# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:7]
print(args)

gxp_metadata <- as.character(args[1]) # Path to phenotype PCA files folder
# gxp_metadata <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/hamlet_color/metadata/metadata_gxp_ben_floridae_complete"
coverage_table <- as.character(args[2]) # Path to coverage table file
# coverage_table <- "/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/hamlet_color/coverage_table"


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
                      label = gsub("PL17_23nigpue", "PL17_23tanpue", label),
                      spec = ifelse(label == "PL17_35puepue", "ind", spec),
                      label = gsub("PL17_35puepue", "PL17_35indpue", label)) %>%
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
               summarise(label, spec, geo, date, coord_N.x, coord_W.x) %>%
               setNames(., nm = c("SampleID", "spec", "geo", "Date", "coord_N.x", "coord_W.x"))

cov_list <- read.csv(coverage_table, sep=" ", header = FALSE, col.names = c("SampleID", "coverage")) %>%
            mutate(SampleID = gsub("PL17_23nigpue", "PL17_23tanpue", SampleID),
                   SampleID = gsub("PL17_35puepue", "PL17_35indpue", SampleID),
                   coverage = round(coverage, digits = 1)) %>%
            filter(., !SampleID %in% c("PL17_101maybel", "AG9RX_47pueboc", "PL17_98indbel", "PL17_79abepue"))

cov_sample <- left_join(sample_list, cov_list, by = "SampleID")

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

geo_label <- c("abe" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_aberrans.l.cairo.png' width='60' /><br> *H. aberrans*",
               "chl" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_chlorurus.l.cairo.png' width='60' /><br>*H. chlorurus*",
               "flo" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_floridae.l.cairo.png' width='60' /><br>*H. floridae*",
               "gem" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_gemma.l.cairo.png' width='60' /><br>*H. gemma*",
               "gum" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_gumigutta.l.cairo.png' width='60' /><br>*H. gummigutta*",
               "gut" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_guttavarius.l.cairo.png' width='60' /><br>*H. guttavarius*",
               "ind" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_indigo.l.cairo.png' width='60' /><br>*H. indigo*",
               "may" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_maya.l.cairo.png' width='60' /><br>*H. maya*",
               "nig" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_nigricans.l.cairo.png' width='60' /><br>*H. nigricans*",
               "pue" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_puella.l.cairo.png' width='60' /><br>*H. puella*",
               "ran" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_randallorum.l.cairo.png' width='60' /><br>*H. randallorum*",
               "uni" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_unicolor.l.cairo.png' width='60' /><br>*H. unicolor*",
               "tan" = "<img src='/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/logos/H_sp.l.cairo.png' width='60' /><br>*Hypoplectrus* sp.")


# Create the plot
p <- ggplot(data = world, aes(x=long, y=lat)) +
     geom_map(map=world, aes(map_id=region), fill="peachpuff1", color = "antiquewhite") +
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
                       labels = geo_label) +
     theme(legend.title=element_blank(),
           # text = element_text(size = 3),
           panel.background = element_rect(fill = "lightblue2"),
           rect = element_rect(fill = "transparent"),
           panel.grid.major = element_blank(),
           legend.text = element_markdown(size = 10),
           legend.key=element_blank(),
           legend.box = "vertical",
           axis.title = element_blank(),
           axis.ticks = element_blank(),
           axis.text = element_blank(),
           panel.border = element_blank(),
           plot.background=element_rect(fill="transparent", colour=NA)) +
  guides(fill = guide_legend(nrow = 7))

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






sample_list2 <- cov_sample %>%
                mutate(ID = substr(SampleID, 1, nchar(SampleID)-6),
                       Nr = seq.int(nrow(sample_list))) %>%
                summarise(Nr, ID, spec, geo, Date, coord_N.x, coord_W.x, coverage) %>%
                setNames(., nm = c("Nr", "ID", "Species", "Location", "Date", "Latitude", "Longitude", "Coverage"))


ben_ena <- data.frame("ID" = c("28366", "28377", "28384", "28386", "28387", "28388", "28390", "28391", "28392", "28394", "28399", "AG9RX_46", "AG9RX_48", "AG9RX_49", "AG9RX_50", "AG9RX_51", "AG9RX_53", "PL17_01", "PL17_02", "PL17_04", "PL17_05", "PL17_103", "PL17_100", "PL17_104", "PL17_105", "PL17_106", "PL17_107", "PL17_108", "PL17_109", "PL17_110", "PL17_111", "PL17_112", "PL17_125", "PL17_117", "PL17_127", "PL17_128", "PL17_132", "PL17_134", "PL17_135", "PL17_136", "PL17_137", "PL17_138", "PL17_139", "PL17_140", "PL17_141", "PL17_143", "PL17_149", "PL17_155", "PL17_157", "PL17_159", "PL17_23", "PL17_35", "PL17_37", "PL17_38", "PL17_39", "PL17_40", "PL17_41", "PL17_42", "PL17_43", "PL17_44", "PL17_50", "PL17_53", "PL17_55", "PL17_54", "PL17_56", "PL17_57", "PL17_60", "PL17_62", "PL17_63", "PL17_64", "PL17_65", "PL17_67", "PL17_68", "PL17_69", "PL17_70", "PL17_71", "PL17_72", "PL17_73", "PL17_74", "PL17_93", "PL17_91", "PL17_90", "PL17_88", "PL17_86", "PL17_87", "PL17_75", "PL17_76", "PL17_77", "PL17_82", "PL17_85", "PL17_97", "PL17_94", "PL17_96", "PL17_99", "28383", "28385", "28389","PL17_66","PL17_89", "PL17_95", "PL17_119", "PL17_120", "PL17_121",
                               "PL17_122", "PL17_123", "PL17_124", "PL17_126", "PL17_142",
                               "PL17_144", "PL17_145", "PL17_148", "PL17_153", "PL17_160"),
                      "Accession Number" = c("ERS8632035", "ERS14948427", "ERS14948442", "ERS14948430", "ERS14948426", "ERS14948436", "ERS14948425", "ERS14948444", "ERS14948440", "ERS14948432", "ERS14948433", "ERS14948429", "ERS14948431", "ERS14948424", "ERS14948437", "ERS14948434", "ERS14948435", "ERS14948439", "ERS14948428", "ERS14948438", "ERS14948441", "ERS14948415", "ERS14948407", "ERS14948408", "ERS14948419", "ERS14948401", "ERS14948420", "ERS14948411", "ERS14948417", "ERS14948422", "ERS14948399", "ERS14948405", "ERS14948421", "ERS14948404", "ERS14948402", "ERS14948406", "ERS14948416", "ERS14948448", "ERS14948455", "ERS14948452", "ERS14948453", "ERS14948451", "ERS14948450", "ERS14948456", "ERS14948457", "ERS14948454", "ERS14948447", "ERS14948449", "ERS14948446", "ERS14948445", "ERS14948480", "ERS14948482", "ERS14948469", "ERS14948464", "ERS14948460", "ERS14948466", "ERS14948462", "ERS14948467", "ERS14948461", "ERS14948463", "ERS14948475", "ERS14948476", "ERS14948485", "ERS14948483", "ERS14948484", "ERS14948479", "ERS14948481", "ERS14948465", "ERS14948487", "ERS14948471", "ERS14948490", "ERS14948493", "ERS14948470", "ERS14948486", "ERS14948489", "ERS14948478", "ERS14948472", "ERS14948491", "ERS14948492", "ERS14948412", "ERS14948400", "ERS14948418", "ERS14948398", "ERS14948468", "ERS14948413", "ERS14948458", "ERS14948477", "ERS14948488", "ERS14948473", "ERS14948474", "ERS14948409", "ERS14948410", "ERS14948414", "ERS14948403", "ERS14948443", "ERS14948423", "ERS8632036", "ERS14948459", "ERS2899590", "ERS2899591", "ERS2899593",
                                             "ERS2899594", "ERS2899595", "ERS2899596",
                                             "ERS2899597", "ERS2899598", "ERS2899599",
                                             "ERS2899137", "ERS2899138", "ERS2899139",
                                             "ERS2899140", "ERS2899141", "ERS4141277"))

tableS1 <- merge(sample_list2, ben_ena, all = TRUE) %>%
           arrange(Nr) %>%
           summarise(Nr, ID, Species, Location, Date, Latitude, Longitude, Coverage, Accession.Number) %>%
           setNames(., nm = c("Nr", "ID", "Species", "Location", "Date", "Latitude", "Longitude", "Coverage", "Accesion Number"))

tableS1 <- tableS1 %>% mutate(Species = sub("tan", "sp.", Species),
                              Coverage = ifelse(ID == "PL17_111", "26.4", Coverage),
                              Coverage = ifelse(ID == "PL17_108", "24.6", Coverage))

tableS1[is.na(tableS1)] <- "_"
table_h1 <- tableS1[1:57,]
table_h1
table_h2 <- tableS1[58:113,]
table_h2

# pdf("/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/figures/Tableh1.pdf", height = 17, width = 10)
# # p<-tableGrob(cmp_glob)
# grid.table(table_h1, rows = NULL)
# # grid.table(table_h2, rows = NULL)
# dev.off()

# pdf("/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/figures/Tableh2.pdf", height = 17, width = 10)
# # p<-tableGrob(cmp_glob)
# # grid.table(table_h1, rows = NULL)
# grid.table(table_h2, rows = NULL)
# dev.off()


# h <- knitr::include_graphics(c("/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/figures/Tableh1.pdf", "/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/chapter1/figures/Tableh2.pdf"))
# h <- knitr::kable(list(table_h1, table_h2))
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 2.0)),
  colhead = list(fg_params=list(cex = 2.0)),
  rowhead = list(fg_params=list(cex = 2.0)))


Table1 <- tableGrob(table_h1, rows = NULL, theme = mytheme)
Table2 <- tableGrob(table_h2, rows = NULL, theme = mytheme)
pdf("/Users/fco/Desktop/PHD/1_CHAPTER1/1_GENETICS/hamlet_color/figures/TableS1.pdf", height = 24, width = 32)
# p<-tableGrob(cmp_glob)
# grid.table(table_h1, rows = NULL)

grid.arrange(Table1, Table2, ncol = 2)
# grid.table(h, rows = NULL)
dev.off()



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
