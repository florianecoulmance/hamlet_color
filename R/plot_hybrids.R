#!/usr/bin/env Rscript
# by: Floriane Coulmance: 02/02/2022
# usage:
# Rscript plot_hybrids.R <newhybrid_directory> <result_directory>
# -------------------------------------------------------------------------------------------------------------------
# newhybrid_directory : $BASE_DIR/outputs/9_newhyb/
# result_directory : $BASE_DIR/figures/9_newhyb/
# -------------------------------------------------------------------------------------------------------------------


# Clear the work space
rm(list = ls())

# Load needed library
library(paletteer)
library(ggthemes)
library(prismatic)
library(patchwork)
library(ggtext)
library(hypoimg)
library(stringr)


# -------------------------------------------------------------------------------------------------------------------
# ARGUMENTS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# Get the arguments in variables
args = commandArgs(trailingOnly=FALSE)
args = args[6:7]
print(args)

base_dir <- as.character(args[1]) # Path to the general NEWHYBRID result folder
print(base_dir)
fig_path <- as.character(args[2]) # Path to figure folder
print(fig_path)

# -------------------------------------------------------------------------------------------------------------------
# FUNCTIONS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


get_subfolder <- function(path, dir_list) {
  
  # Takes a list of directory and find path to each of their subsirectories
  
  res_folders <- list() # Create empty list to store path to folders
  i <- 0
  print(i)

  # Loop through list of folder in base directory to retrieve important sub folders
  for (folder in dir_list) {
    
    print(folder)
    label <- str_split(folder, fixed("_")) # Split folder name to get each element separated by "_"
    print(label)
    print(label[[1]][1])
    print(label[[1]][2])
    print(label[[1]][3])
    pop1 <- label[[1]][1] # Get 1st population from the pairwise comparison
    print(pop1)
    pop2 <- label[[1]][2] # Get 2nd population from the pairwise comparison
    print(pop2)
    
    runname <- paste0(pop1,"_",pop2) # Get the pairwise comparison label
    print(runname)
    
    result_folder <- paste0(path,folder,"/NH.Results/newHyb.",runname,".80SNPs.txt_Results/")
    print(result_folder)

    i <- i+1
    print(i)
    res_folders[[i]] <- result_folder
    print(res_folders)

    
  }
  
  print(res_folders)
  return(res_folders)
  
}


getPofZ <- function(res_paths) {
  
  # Takes a list of subdirectories, for each of them extract files and make a table of data from these extracted files
  
  # Get the pairwise comparison label
  runname <- res_paths %>% 
             #str_remove("newHyb.") %>% 
             str_remove(".80SNPs.txt_Results/") %>%
             sub("^.+newHyb",) %>%
             str_remove(".")
  print(runname)
  
  pops <- c(str_sub(runname,1,6),str_sub(runname,-6,-1)) # Get each population of the pairwise comparison separately
  print(pops)
  
  pofz <- dir(res_paths, pattern = "PofZ.txt") # List the Pofz files for each pairwise comparison of the list 
  print(pofz)
  inds <- dir(res_paths, pattern = "_individuals.txt") # List the individual files for each pairwise comparison of the list
  print(inds) 
  
  colN <- c("P1", "P1_bc", "P2", "P2_bc") # Set future column names
  
  # Open Pofz files and combine them in big table and add list of individual sample names to column from individuals files 
  NHres <- vroom::vroom(str_c(res_paths, pofz),
                        delim = '\t',
                        skip = 1,
                        col_names  = c('indNR', 'IndivName', colN[1], colN[3], 'F1', 'F2', colN[2], colN[4])) %>%
    mutate(IndivName = vroom::vroom(str_c(res_paths, inds), 
                                    delim = ',',
                                    col_names =  c('IndivName'))[,1] %>%
             unname() %>%
             unlist())
  
  # Rename column names for big table
  bin_tib <- tibble(bin_generic = c(colN, "F1", "F2"),
                    bin = c(paste0(c(pops[1],pops[1],pops[2],pops[2]), c('_pure','_BC','_pure','_BC')), "F1", "F2"))
  
  # Create final data table
  data <- NHres %>%
    pivot_longer(names_to = 'bin_generic',
                 values_to= "post_prob",
                 cols = c(-indNR, -IndivName)) %>%
    left_join(bin_tib) %>%
    mutate(run = runname,
           loc = str_sub(run,-3,-1),
           ind_order = str_c(str_sub(IndivName,-6,-1),"_",str_sub(IndivName,1,-7)))
  
  return(data)
  
}


plot_loc <- function(loc){
  
  # Takes a list of files and subset according to locations to create a plot for the considered location
  
  print(fold_paths[str_sub(fold_paths,-23,-21) == loc])
  
  # Create big data table combining all files from different pairwise comparison and keep track of pairwise comparison label 
  data <- map_dfr(.x = fold_paths[str_sub(fold_paths,-23,-21) == loc],
                  .f = getPofZ)
  print(data)
  
  # Find samples that are hybrids
  is_hybr <- data %>%
    filter(!(grepl(pattern = "_pure", bin))) %>%
    filter(post_prob > .99) %>%
    .$ind_order %>%
    unique()
  print(is_hybr)
  
  # List of species names and their abbreviation
  sp_names <- c(abe = "aberrans",
                flo = "floridae",
                gum = "gummigutta",
                ind = "indigo",
                may = "maya",
                nig = "nigricans",
                pue = "puella",
                ran = "randallorum",
                tab = "tabacarius",
                tor = "tortugarum",
                uni = "unicolor",
                chl = "chlorurus")
  
  data <- data %>%
          mutate(ind_label = ifelse(ind_order %in% is_hybr, str_c("**", ind_order, "**"), ind_order),
                 run2 = str_c("*H. ", sp_names[str_sub(string = run,1,3)],"* - *H. ", sp_names[str_sub(string = run,8,10)],"*"))
  print(data)
  
  # Create labels
  data_labs <- data %>%
               filter(!duplicated(ind_order)) %>%
               select(ind_order, ind_label)
  print(data_labs)
  
  # Create color palette for hybrid, pure line and backcrosses
  clr <- paletteer_c("ggthemes::Red-Green-Gold Diverging",3) %>%
         c(.,clr_lighten(.)) %>%
         color() %>% 
         .[c(1,4,2,5,6,3)] %>%
         set_names(nm = c("P1", "P1_bc","F1", "F2", "P2_bc", "P2"))
  print(clr)
  
  # List of location names and labels
  loc_names <- c(bel = "Belize",
                 pue = "Puerto Rico",
                 boc = "Panama",
                 flo = "Florida")
  
  # Create plot panel for a particular location 
  data  %>%
    ggplot(aes(x = ind_order, y = post_prob, fill = bin_generic)) +
    geom_bar(position = 'stack', stat = "identity") +
    scale_fill_manual(values = clr) +
    scale_x_discrete(breaks = data_labs$ind_order, labels = data_labs$ind_label) +
    labs(y = "Posterior probability", title = loc_names[loc]) +
    facet_grid(run2~.) +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text.x  = element_markdown(),
          strip.text.y  = element_markdown(),
          axis.text.x = element_markdown(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_text(vjust = 4))
  
}


theme_hyb <-  function(legend.position = "none",...) {
  
  # Define a theme for plot presentation
  
  list(scale_y_continuous(breaks = c(0,.5,1)),
       theme(legend.position = legend.position,
             legend.background = element_rect(fill = "white",colour = rgb(1,1,1,0)),
             legend.direction = "horizontal",
             legend.justification = c(1,1),
             strip.text.y = element_markdown(angle = 0,hjust = 0), ...))
  
}


final_plot <- function(plots, path) {
  
  # Compose figure from the each location panels and save it as a pdf figure 
  
  p <- (plots[[1]] +  guides(fill = guide_legend(title = "Hybrid Class")) + theme_hyb(legend.position = c(1,1))) +
       (plots[[2]] + theme_hyb()) +
       (plots[[3]] + theme_hyb()) +
       plot_layout(ncol = 1, heights = c(10,15,3) %>% label_spacer()) +
       plot_annotation(tag_levels = 'a')
  
  hypo_save(filename = paste0(path,"hybrids.pdf"),
            plot = p,
            height = 16,
            width = 10)
  
}


# -------------------------------------------------------------------------------------------------------------------
# ANALYSIS

# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# sp_labs <- c(abe = expression(italic(H.~aberrans)),
#              flo = expression(italic(H.~floridae)),
#              gum = expression(italic(H.~gummigutta)),
#              ind = expression(italic(H.~indigo)),
#              may = expression(italic(H.~maya)),
#              nig = expression(italic(H.~nigricans)),
#              pue = expression(italic(H.~puella)),
#              ran = expression(italic(H.~randallorum)),
#              tab = expression(italic(S.~tabacarius)),
#              tor = expression(italic(S.~tortugarum)),
#              uni = expression(italic(H.~unicolor)),
#              chl = expression(italic(H.~chlorurus)))

sub_dir <- list.dirs(path = base_dir, full.names = F, recursive = F)
print(sub_dir)

fold_paths <- get_subfolder(base_dir, sub_dir)
print(fold_paths)

p_loc <- c("boc", "bel", "por") %>% map(plot_loc)
print(p_loc)

# bocdat <- p_loc[[1]]$data
# beldat <- p_loc[[2]]$data
# puedat <- p_loc[[3]]$data

final_plot(p_loc, fig_path)
