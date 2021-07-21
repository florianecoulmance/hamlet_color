rm(list = ls())


library("ggspatial")
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("rgeos")
library(tidyverse)

world <- ne_countries(scale = "medium", returnclass = "sf")
table <- read.csv("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/metadata/metadata_gxp_ben_floridae_complete", sep=";")
table2 <- read_tsv("/Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/1_GENOTYPING/samples_meta.txt", col_names = FALSE)
table3 <- table2[!is.na(table2$X7),]
#table <- subset(table, site.country. == "Mexiko")
joined_df <- merge(table, table3, by.x = "id", 
                   by.y = "X1")

joined_df["coord_N"] <- joined_df["X7"]
joined_df["coord_W"] <- joined_df["X8"]
drops <- c("X2","X3","X4","X5","X6","X7","X8","X9","X10","X11")
test <- joined_df[ , !(names(joined_df) %in% drops)]

ggplot(data = world) +
  geom_sf() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-102.15, -74.12), ylim = c(7.65, 33.97))


ggplot(data = world) + geom_sf(fill= "antiquewhite") + 
  geom_point(data = test, mapping = aes(x = coord_W, y = coord_N), color = "red") +
  annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", fontface = "italic", color = "grey22", size = 3) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-102.15, -60.12), ylim = c(5.00, 33.97), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("Map of the Gulf of Mexico and the Caribbean Sea") +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))








