#!/usr/bin/python
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

import sys # to get arguments from command line
import seaborn as sns
from PIL import Image
import os, numpy, PIL
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, ListedColormap, LinearSegmentedColormap
from skimage import data, io, color, filters, util
from skimage.color import rgb2hsv, rgb2lab
from numpy import *
from pylab import *
import pandas as pd
import cv2
from sklearn.decomposition import PCA
import numpy as np
from scratch_1 import blur_mask, modify_image, Pca



# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def make_pc_list(table_name) :
    '''
    Create the string vector corresponding to the PCs taken into account in the multivariate association
    input : path_to_file as string
    '''

    nb1 = int(table_snp.split('_')[0].split('PC')[1]) # <-- split string by _ and take 1st string then split by "PC" and take the following string to have the nb
    # print(type(nb1))
    nb2 = int(table_snp.split('_')[1].split('.')[0]) # <-- split string by _ and take 2nd string then split by . and take the 1st part of the string to have the nb
    # print(type(nb2))

    p = [] # <-- empty list of PCs

    # Iterate over nb of PCs to create PC string list 
    for i in range(nb1, nb2 + 1):
        s = "PC" + str(i)
        print(s)
        p.append(s)

    return(p)



# ----------------------------------------------------------------------------------------------

def make_heat(b_m, rgb_m, pca, weights, pcs):
    '''
    Create the phenotype image for a particular SNP based on multivariate PCs involved in association and their relative importance in the association
    inputs : boolean_mask matrix the size of the image (1 layer only),
             rgb_image corresponding to boolean mask in black and white (raster matrix 3 layers),
             output of PCA with 15 components,
             list of weights corresponding to the nb of PCs in the multivariate association and corresponding to a particular snp
             list of PCs as string corresponding to the PCs involved in the multivariate association
    '''

    # get the egeinvectors for each PCs : matrix = 15 PCs x nb of pixels flattened (nb_pixels * 3)
    feat_imp = pca.components_ 

    # Create empty table to store vectors of calculated values for each PCs considered in the association
    weights_norms = numpy.zeros(shape=(len(pcs),273217))

    # count starts at different number depending on the first PC in the pc list
    # important to consider the right row in the feat_imp table for further analysis 
    if (pcs[0] == "PC1"):
        count = 0
    elif (pcs[0] == "PC2"):
        count = 1
    elif (pcs[0] == "PC3"):
        count = 2
    elif (pcs[0] == "PC4"):
        count = 3
    elif (pcs[0] == "PC5"):
        count = 4
    elif (pcs[0] == "PC6"):
        count = 5
    elif (pcs[0] == "PC7"):
        count = 6
    elif (pcs[0] == "PC8"):
        count = 7
    elif (pcs[0] == "PC9"):
        count = 8
    elif (pcs[0] == "PC10"):
        count = 9
    else:
        print("not implemented")

    
    # Loop through list of PC and perform calculations with weights for each PCs
    # count can be different than i
    for i, pc in enumerate(pcs):
        print(count)
        print(i)
        print(pc)

        im_PC = feat_imp[count].reshape((273217, 3), order='C') # <-- use count as index since we want to retrieve the specific PC row in feat_imp (the matrix of eigenvectors for each PC)

        im_PCscores = [v for v in weights[i]*np.linalg.norm(im_PC,axis = 1)] # <-- use i as index to retrieve the specific PC weight in the list of weights
                                                                             # <-- calculation of the norm of 3 layers of raster per pixel and multiply by the weights
                                                                             # <-- return a vector the lenght of number of pixel 

        weights_norms[i] = im_PCscores # <-- fill the table with the vector just produced
        
        count += 1 # <-- increment count to go to the next PC row

    im_pcs = np.sum(weights_norms,axis=0) # <-- sum the PCs column by column

    # find minimum and maximum values
    min_val, max_val = min(im_pcs), max(im_pcs)
    print(min_val, max_val)
    print(-max_val, -min_val)
    
    # Create colormap for futur transformation of im_pcs to raster
    norm = matplotlib.colors.Normalize(-0.005,0.005)
    colors = [[norm(-0.005), "darkblue"],
          [norm(-0.003), "blue"],
          [norm(-0.001), "lightblue"],
          [norm(0), "grey"],
          [norm(0.001), "yellow"],
          [norm(0.003), "orange"],
          [norm( 0.005), "darkred"]]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    mapper = cm.ScalarMappable(cmap=cmap,norm=norm)

    # transform the im_pcs vector to an RGB image where the pixel as a color proportional to the value in the vector
    color_mask = np.array([(r, g, b) for r, g, b, a in mapper.to_rgba(im_pcs)])

    rgb_im = np.array(rgb_m) # <-- transform the black and white input image as a numpy array
    rgb_im[b_m] = color_mask # <-- fill the positions of the black and white numpy image where the boolean mask is TRUE with the color values

    return rgb_im



# ----------------------------------------------------------------------------------------------

def plot_heat(im_table, figure_path, region_id):
    '''
    Creates a panel of 5 heatmaps per regions corresponding to the 5 highest associated SNPs of the region
    '''

    g = sns.FacetGrid(im_table, col="SNP", col_wrap = 2, height = 4) # <-- create facetgrid per SNP
    
    cax = g.fig.add_axes([0.910, 0.100, 0.009, 0.775]) # <-- create ax for the colorbar scale
    
    bounds = [-0.005, 0, 0.005] # <-- create min and max for the colorbar scale

    norm = matplotlib.colors.Normalize(-0.005,0.005) # <-- norm to use further to create own color map

    # create own colormap gradient
    colors = [[norm(-0.005), "darkblue"],
        [norm(-0.003), "blue"],
        [norm(-0.001), "lightblue"],
        [norm(0), "grey"],
        [norm(0.001), "yellow"],
        [norm(0.003), "orange"],
        [norm( 0.005), "darkred"]]
    
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors) # <-- with gradient create the colormap for heatmaps

    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks = bounds, orientation='vertical') # <-- Create a colorbar axes
    
    g.map(lambda x, **kwargs : (plt.imshow(x.values[0]),plt.grid(False)), # <-- feed each subplot with the heatmap image correponding to each SNP
                                'IMAGE', cbar_ax=cmap, vmin=-0.005, vmax=0.005) 

    g.fig.tight_layout() # <-- tights the subplots between them
    g.despine(left=True, bottom=True) # <-- remove x and y lines
    g.set(xlabel=None) # <-- remove x label
    g.set_xticklabels([]) # <-- remove x tick labels
    g.set_yticklabels([]) # <-- remove y tick labels
    g.set(xticks=[]) # <-- remove x ticks
    g.set(yticks=[]) # <-- remove y ticks 
    g.fig.subplots_adjust(right=.9)  # <-- Add space so the colorbar doesn't overlap the plot

    # plt.show()

    g.savefig(figure_path+region_id+"_heatmaps.png") # <-- save in appropriate figure folder with region id as file title



# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def main():
    '''
    Executes all essential commands
    '''
    
    snp_arr = pd.read_csv(path_data+table_snp, sep=' ') # <-- open the table of SNP highly associated
    arr_modified = pd.read_csv(path_images, sep=',') # <-- open the table of modified images from the dataset
    mask = cv2.imread(path_mask, 0) # <-- open the mask file
    print(snp_arr)
    print(arr_modified)
    
    pc_list = make_pc_list(table_snp) # <-- create the PC string list from the input file string
    
    bool_mask, rgb_mask, mask_blur = blur_mask(mask) # <-- create all the mask matrix and images needed for further analysis

    # arr_modified, list_f = modify_image(path_images, bool_mask, mask_blur, "LAB") # <-- create table of sample * modified image (refer to function in scratch_1.py)

    # save_modifiedImage(arr_modified, list_f, "LAB")
    
    pca_im, new_arr = Pca(arr_modified) # <-- perform PCA on the sample * image created at the previous step
    
    # variances(pca_im, "LAB")
    # pc1_im = plot_heatmap(bool_mask, rgb_mask, pca_im, 1, "AB")
    # pc2_im = plot_heatmap(bool_mask, rgb_mask, pca_im, 1, "AB")
    # pc3_im = plot_heatmap(bool_mask, rgb_mask, pca_im, 1, "AB") 

    region_id = snp_arr["ID"].values.tolist() # <-- retrieve the region id column
    region_list =  list(set(region_id)) # <-- make a list of regions 

    # loop over regions list to create plot with 5 SNPs for each region in the list
    for i, region in enumerate(region_list) :
        print(i)
        print(region)
        df = snp_arr[snp_arr["ID"].str.contains(region)] # <-- find rows in table that have the same region id and make a subset of this table for the region

        wei = df["WEIGHTS"].values.tolist() # <-- make a list of weight list from the WEIGHT column of the region table
        print(wei)

        pos = df["POS"].values.tolist() # <-- make a list of 5 positions from each row corresponding to the 5 SNPs

        columns = ["SNP", "IMAGE"] # <-- create column names
        heat_table = pd.DataFrame(index = [1,2,3,4,5], columns = columns) # <-- create empty table that will store the 5 SNP ids and the corresponding images

        # loop over the list of weight lists 
        for i, load in enumerate(wei) :
            print(i)
            print(load)

            lab = region + '_' + str(pos[i]) # <-- create label for each SNP consisting of the chromosome region id and the position in this region
            heat_table["SNP"][i+1] = lab # <-- store the label in the table 

            w = wei[i].split(',') # <-- make the string of weight a list separated by ,
            w = list(map(float, w)) # <-- change the list os string separated by , to a list of floats
            print(w)

            im = make_heat(bool_mask, rgb_mask, pca_im, w, pc_list) # <-- create image heatmap based on weights list for the SNP

            heat_table["IMAGE"][i+1] = im # <-- store the image in the table on the row corresponding to the SNP


        plot_heat(heat_table, path_figure, region) # <-- plot 5 images per chromosome region id  



# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

if __name__ == "__main__":

    print ("The script is called %s" % (sys.argv[0])) # <-- name of the script running
    arguments = len(sys.argv) - 1 # <-- count the arguments
    print ("The script is called with %i arguments" % (arguments))
    
    path_data = sys.argv[1]    
    table_snp = sys.argv[2] # <-- path to the input table of highest SNP associated
    # print(table_snp) 
    path_images = sys.argv[3] # <-- path to the modified image table data
    path_mask = sys.argv[4] # <-- path to the mask file
    path_figure = sys.argv[5] # <-- path to the figure folder
    
    print(path_data)
    print(table_snp)
    print(path_images)
    print(path_mask)
    print(path_figure)

    main() # <-- execute the main






# (testA==testB).all() check if two numpy arrays are equivalent
