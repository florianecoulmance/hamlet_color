#!/usr/bin/python

import sys
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
matplotlib.use('TkAgg')
from scratch_1 import blur_mask, modify_image, Pca



# ----------------------------------------------------------------------------------------------

def make_pc_list(table_name) :
    '''
    Create the string vector corresponding to the PCs taken into account in the multivariate association
    '''

    nb1 = int(table_snp.split('_')[0].split('PC')[1])
    # print(type(nb1))
    nb2 = int(table_snp.split('_')[1].split('.')[0])
    # print(type(nb2))
    p = []
    for i in range(nb1, nb2 + 1):
        s = "PC" + str(i)
        print(s)
        p.append(s)
    # print(p)

    return(p)


# ----------------------------------------------------------------------------------------------

def make_heat(b_m, rgb_m, pca, weights, pcs):
    '''
    Create the fish heatmap corresponding to the importance of each pixels in the PCs variations
    '''

    # get the egeinvectors for each PCs : matrix = PCs x pixels
    feat_imp = pca.components_ 
    # print(feat_imp)
    # print(feat_imp[0])
    # print(feat_imp.shape)

    # print(len(pcs))

    # Create empty list to store norm values 
    weights_norms = numpy.zeros(shape=(len(pcs),273217))
    # print(weights_norms.shape)

    # count starts at different number depending on the first PC in the list
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
    else:
        print("not implemented")

    

    for i, pc in enumerate(pcs):
        print(count)
        print(i)
        print(pc)
        im_PC = feat_imp[count].reshape((273217, 3), order='C') # 347361
        # print(weights[count])
        # print(im_PC)
        im_PCscores = [v for v in weights[count]*np.linalg.norm(im_PC,axis = 1)]
        # print(len(im_PCscores))
        # list_norms.append(im_PCscores)
        weights_norms[i] = im_PCscores
        count += 1

    # print(weights_norms.shape)
    # print(weights_norms)

    im_pcs = np.sum(weights_norms,axis=0)
    # print(im_pcs.shape)  
    # print(len(im_pcs))

    min_val, max_val = min(im_pcs), max(im_pcs)
    print(min_val, max_val)
    print(-max_val, -min_val)
    
    # cmap = cm.icefire
    # cmap = sns.color_palette("Spectral", as_cmap=True)

    # if abs(max_val)<abs(min_val):
    #     bounds = [min_val, -max_val, 0, max_val, -min_val]
    #     norm = cm.colors.Normalize(vmin=min_val, vmax=-min_val)
    # elif abs(max_val)>abs(min_val):
    #     bounds = [-max_val, -min_val, 0, min_val, max_val]
    #     norm = cm.colors.Normalize(vmin=-(max_val), vmax=max_val)
    
    # norm = cm.colors.Normalize(vmin=-0.005, vmax=0.005)
    
    bounds = [-0.005, 0, 0.005]
    norm = matplotlib.colors.Normalize(-0.005,0.005)
    colors = [[norm(-0.005), "darkblue"],
          [norm(-0.003), "blue"],
          [norm(-0.001), "lightblue"],
          [norm(0), "grey"],
          [norm(0.001), "yellow"],
          [norm(0.003), "orange"],
          [norm( 0.005), "darkred"]]

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    # print(bounds)
    mapper = cm.ScalarMappable(cmap=cmap,norm=norm)



    # print(mapper)
    # color_list = cmap(im_PCscores)
    # print(len(color_list))
    # print(color_list[0:10])
    color_mask = np.array([(r, g, b) for r, g, b, a in mapper.to_rgba(im_pcs)])
    # print(len(color_mask))
    # print(color_mask[0:10])
    rgb_m[b_m] = color_mask
    print(rgb_m.shape)
    
    return rgb_m

    # fig, ax = plt.subplots(figsize=(12,8))
    # im = a.imshow(rgb_m)
    # cax = f.add_axes([0.910, 0.100, 0.009, 0.775])
    # cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks = bounds, orientation='vertical')
    # a.set_axis_off()  
    # plt.xticks([])
    # plt.yticks([]) 
    #cb.ax.tick_params(labelsize=5)
    # cb.ax.ticklabel_format(style='scientific',useOffset=True, useMathText=True)
    # ax.xaxis.set_major_formatter(cb.ScalarFormatter(useMathText=True))
    #cb.ax.set_yticklabels(['{:e}'.format(x) for x in bounds])

    # return feat_imp


# ----------------------------------------------------------------------------------------------

def plot_heat(im_table):
    '''
    '''

    g = sns.FacetGrid(im_table, col="SNP")
    # g.map(plt.imshow(im_table["IMAGE"]))
    cax = g.fig.add_axes([0.910, 0.100, 0.009, 0.775])
    
    bounds = [-0.005, 0, 0.005]

    norm = matplotlib.colors.Normalize(-0.005,0.005)

    colors = [[norm(-0.005), "darkblue"],
        [norm(-0.003), "blue"],
        [norm(-0.001), "lightblue"],
        [norm(0), "grey"],
        [norm(0.001), "yellow"],
        [norm(0.003), "orange"],
        [norm( 0.005), "darkred"]]
    
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks = bounds, orientation='vertical') # <-- Create a colorbar axes
    
    g.map(lambda x, **kwargs : (plt.imshow(x.values[0]),plt.grid(False)), 
                                'IMAGE', cbar_ax=cmap, vmin=-0.005, vmax=0.005)
                                
    # g.map(plt.imshow(x.values[0]), 'IMAGE')
    g.set_yticklabels([])
    g.set_xticklabels([])
    g.fig.subplots_adjust(right=.9)  # <-- Add space so the colorbar doesn't overlap the plot

    plt.show()

    # g.savefig("output.png")



# ----------------------------------------------------------------------------------------------

def main():
    
    pc_list = make_pc_list(table_snp)
    # print(pc_list)
    
    bool_mask, rgb_mask, mask_blur = blur_mask(mask)
    # print(bool_mask)
    # print(rgb_mask)
    # print(bool_mask.shape)
    # print(rgb_mask.shape)

    arr_modified, list_f = modify_image(path_images, bool_mask, mask_blur, "LAB")
    # print(arr_modified)
    # print(arr_modified.shape)

    # save_modifiedImage(arr_modified, list_f, "LAB")
    
    pca_im, new_arr = Pca(arr_modified)
    # print(pca_im)
    
    # variances(pca_im, "LAB")
    # pc1_im = plot_heatmap(bool_mask, rgb_mask, pca_im, 1, "AB")
    # pc2_im = plot_heatmap(bool_mask, rgb_mask, pca_im, 1, "AB")
    # pc3_im = plot_heatmap(bool_mask, rgb_mask, pca_im, 1, "AB") 



    # retrieve the chrom column and make a list of chromosomes 
    chrom = snp_arr["CHROM"].values.tolist()
    # print(chrom)
    chrom_list =  list(set(chrom))
    # print(chrom_list)



    for i, ch in enumerate(chrom_list) :
        print(i)
        print(ch)
        df = snp_arr[snp_arr['CHROM'].str.contains(ch)]
        # print(df)

        # retrieve the weights column
        wei = df["WEIGHTS"].values.tolist()
        print(wei)

        pos = df["POS"].values.tolist()


        # fig, ax = plt.subplots(3, 2, figsize=(30,20))
        # print(fig)
        # print(ax)

        columns = ["SNP", "IMAGE"]
        heat_table = pd.DataFrame(index = [1,2,3,4,5], columns = columns)
        # print(heat_table)

        for i, load in enumerate(wei) :
            print(i)
            print(load)
            lab = ch + '_' + str(pos[i])
            heat_table["SNP"][i+1] = lab
            print(heat_table)

            w = wei[i].split(',')
            w = list(map(float, w))
            print(w)

            # if (i == 0) :
            #     b = ax[0][0]
            # elif (i == 1) :
            #     b = ax[0][1]
            # elif (i == 2) :
            #     b = ax[1][0]
            # elif (i == 3) :
            #     b = ax[1][1]
            # else :
            #     b = ax[2][0]


            im = make_heat(bool_mask, rgb_mask, pca_im, w, pc_list)
            heat_table["IMAGE"][i+1] = im
            print(heat_table)


        plot_heat(heat_table) 
        # plt.show()



    # # save_pcpict(new_arr, list_f, "LAB")
    


# ----------------------------------------------------------------------------------------------


if __name__ == "__main__":
    # execute only if run as a script
    print ("The script is called %s" % (sys.argv[0]))
    # Count the arguments
    arguments = len(sys.argv) - 1
    print ("The script is called with %i arguments" % (arguments))
    table_snp = sys.argv[1]
    print(table_snp)
    snp_arr = pd.read_csv(table_snp, sep=' ')
    print(snp_arr)
    path_images = "/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/smallDatasetl1/3-registred/Modalities/RGB/all/"
    mask = cv2.imread("/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/smallDatasetl1/mask1.tif", 0)
    main()