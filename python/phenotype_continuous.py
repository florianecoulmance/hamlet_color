#!/usr/bin/python
# by: Floriane Coulmance: 11/05/2021
# usage : 
# pythonm3 python/phenotype_continuous.py <IMAGE_PATH> <MASK_PATH> <COLOR_SPACE> <OUTPUT_PATH> <FIGURE_PATH> <MASK_LABEL> <DATA_LABEL>
# ----------------------------------------------------------------------------------------------
# IMAGE_PATH : $BASE_DIR/ressources/images/left_54off_59on/3-registred/Modalities/RGB/(all, on, off, barred, unbarred)/
# MASK_PATH : $BASE_DIR/ressources/images/(body_mask.tif, full_mask.tif)
# COLOR_SPACE : LAB, L, AB
# OUTPUT_PATH : $BASE_DIR/images/continuous/LAB/$DATASET/
# FIGURE_PATH : $BASE_DIR/figures/7_gxp/continuous/LAB/$DATASET/
# MASK_LABEL : bodym or fullm
# DATA_LABEL : left_54off, left_59on, left_54off_59on, left_38barred, left_75unbarred
# ----------------------------------------------------------------------------------------------



# ********** Import all necessary library *********
# -------------------------------------------------

from PIL import Image
import seaborn as sns
import sys
import os, numpy, PIL
import matplotlib
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
from matplotlib import gridspec



# FUNCTIONS ------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def blur_mask(m):
    '''
    Transformation of the mask into a boolean mask and blur the mask
    '''

    mask = m > 128

    mask2 = np.full((1000, 1639, 3), np.NaN)
    mask2[mask == 0] = (255, 255, 255)
    mask2[mask == 1] = (0, 0, 0) 
    # cv2.imshow('mask2', mask2)

    mask_blurred = cv2.blur(mask2,(5,5),cv2.BORDER_DEFAULT)
    # cv2.imshow('mask blur', mask_blurred)

    return mask, mask2, mask_blurred



# ----------------------------------------------------------------------------------------------

def modify_image(im_path, bo_mask, bl_mask, effect, path_fig):
    '''
    Perform blur and color space modification on aligned images
    Creates a flatten modified images and store it into a vector
    Store all the flatten images into a table of flatten images
    '''

    file_list = [f for f in os.listdir(im_path) if f.endswith('.png')]
    
    print(bo_mask.sum())
    
    if (effect == "LAB") or (effect == "RGB") or (effect == "HSV"):
        pixel_array = np.full((len(file_list), bo_mask.sum() * 3), np.NaN)
    elif (effect == "L") or (effect == "A") or (effect == "B"):
        pixel_array = np.full((len(file_list), bo_mask.sum()), np.NaN) # 273217 347361
    elif (effect == "AB"):
        pixel_array = np.full((len(file_list), bo_mask.sum() *2), np.NaN) # 546434 694722
    else:
        print("specify correct effect")

    # print(pixel_array.shape)

    for i, f in enumerate(file_list):
        print(i, f)
        # print(im_path,f)
        # print(path_images+f)
        image = cv2.imread(im_path+f)
        # print(image.shape)
        # print(image)
        image_blurred = image.copy()
        # print(image_blurred.shape)

        image_blurred[bo_mask == 0] = (11, 102, 35)
        # print(image_blurred.shape)
        # image_blurred[bo_mask == 1] = (11, 102, 35)

        if (f == "PL17_109puebel-l2-s4-f4-c2-d1.png"):
            # im2 = Image.fromarray(image_blurred)
            # im2.save(path_fig+f, format="png")
            cv2.imwrite(path_fig+f, image_blurred)
             

        image_blurred = cv2.GaussianBlur(image_blurred,(5,5),cv2.BORDER_DEFAULT)
        # print(image_blurred.shape)

        image_blurred[bo_mask == 0] = image_blurred[bo_mask == 0] / bl_mask[bo_mask == 0]
        # print(image_blurred.shape)

        if (effect == "LAB") or (effect == "AB") or (effect == "L") or (effect == "A") or (effect == "B"):
            im = cv2.cvtColor(image_blurred, cv2.COLOR_BGR2LAB)
        elif (effect == "HSV"):
            im = cv2.cvtColor(image_blurred, cv2.COLOR_BGR2HSV)
        elif (effect == "RGB"):
            im = image_blurred
        else:
            print("specify a correct effect")
        
        # print(im.shape)
        # print(im[mask].shape)
        # print(im[mask])

        if (effect == "L"):
            im = im[:, :, 0]
            # print(im.shape)

        if (effect == "A"):
            im = im[:, :, 1]

        if (effect == "B"):
            im = im[:, :, 2]

        if (effect == "AB"):
            im = im[:, :, 1:3]
            # print(ab.shape)
            # a = lab[:, :, 1]
            # b = img[:, :, 2]

        # deroulement de l'image
        pixel_array[i, :] = im[bo_mask].ravel()
    
    # print(pixel_array.shape)
    
    return pixel_array, file_list



# ----------------------------------------------------------------------------------------------

def save_modifiedImage(pix_arr, bo_mask, lof, effect, res_path, data_name):
    '''
    Save the flatten & modified images table with sample names
    '''

    if (effect == "LAB") or (effect == "RGB") or (effect == "HSV"):
        n_col = bo_mask.sum() * 3  # 819651 1042083
    elif (effect == "L") or (effect == "A") or (effect == "B"):
        n_col = bo_mask.sum() # 273217 347361
    elif (effect == "AB"):
        n_col = bo_mask.sum() * 2 # 546434 sq694722
    else:
        print("specify a correct effect")


    columns_name = ['{}'.format(i+1) for i in range(n_col)]
    pixel_table = pd.DataFrame(pix_arr, columns=columns_name)
    pixel_table["images"] = lof
    # print(pixel_table.shape)
    # print(pixel_table)
    pixel_table.to_csv(res_path+data_name+'_modifiedImage.csv')



# ----------------------------------------------------------------------------------------------

def Pca(pix_arr):
    '''
    Define the number of components to perform PCA, save the reconstituted PC images in matrix
    '''

    im_pca = PCA(n_components=15) # <-- perform PCA with 15 components
    # print(im_pca)
    pixel_new = im_pca.fit_transform(pix_arr) # <-- store reconstructed images for each PCs
    # print(pixel_new.shape)

    return im_pca, pixel_new



# ----------------------------------------------------------------------------------------------

def save_pcpict(pixels, flist, res_path, data_name):
    '''
    Save the reconstructed PC images table to result directory
    '''

    n_col = 15 # <-- determine number of columns
    columns_name = ['PC{}'.format(i+1) for i in range(n_col)] # <-- determine names of columns based on nb of columns
    principalDF = pd.DataFrame(data=pixels[:, :n_col], columns=columns_name) # <-- create the dataframe with images and column names
    principalDF["images"] = flist # <-- add the sample columns corresponding to the images
    # print(principalDF.shape)
    # print(principalDF)
    principalDF.to_csv(res_path+data_name+'_PCs.csv') # <-- save the dataframe to result folder with specific name



# ----------------------------------------------------------------------------------------------

def variances(pca, res_path, data_name):
    '''
    Save the variance explained for each PCs and the cumulative variance
    '''

    var_explained = pca.explained_variance_ratio_ # <-- gets the variance explained by PCs
    var_cum = pca.explained_variance_ratio_.cumsum() # <-- gets the cumulative variance explained 

    varDF = pd.DataFrame(var_explained) # <-- makes a dataframe of explained variance per PCs
    varDF.to_csv(res_path+data_name+'_var.csv') # <-- saves the dataframe of explained variance by PCs to result directory

    varcumDF = pd.DataFrame(var_cum) # <-- makes a dataframe of cumulative explained variance
    varcumDF.to_csv(res_path+data_name+'_varcum.csv') # <-- save the cumulative explained variance table to result diretory



# ----------------------------------------------------------------------------------------------

def plot_heatmap(b_m, rgb_m, pca, component, effect, res_path, data_name):
    '''
    Create the fish heatmap corresponding to the importance of each pixels in the PCs variations
    '''

    feat_imp = pca.components_ 

    pc_index = component - 1
    print(pc_index)


    if ((effect == "LAB") or (effect == "RGB") or (effect == "HSV")):
        im_PC = feat_imp[pc_index].reshape((b_m.sum(), 3), order='C') # 347361

    elif (effect == "AB"):
        im_PC = feat_imp[pc_index].reshape((b_m.sum(), 2), order='C')

    elif ((effect == "L") or (effect == "A") or (effect == "B")):
        im_PC = feat_imp[pc_index].reshape((b_m.sum(), 1), order='C')

    else:
        print("not implemented")


    im_PCscores = [v for v in np.linalg.norm(im_PC,axis = 1)]


    min_val, max_val = min(im_PCscores), max(im_PCscores)
    print(min_val, max_val)
    
    norm = matplotlib.colors.Normalize(0,0.006)
    
    boundaries = [0, 0.002, 0.004, 0.006]
    hex_colors = sns.color_palette("rocket", n_colors=10).as_hex()
    print(hex_colors)
    # hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 4)]
    # print(hex_colors)
    # colors=list(zip(boundaries, hex_colors))
    # print(colors)

    colors = [[norm(0), "grey"], #grey
          [norm(0.002), "yellow"], #yellow
          [norm(0.004), "orange"], #orange
          [norm(0.006), "darkred"]] #darkred

    # if abs(max_val)<abs(min_val):
    #     bounds = [min_val, -max_val, 0, max_val, -min_val]
    #     norm = cm.colors.Normalize(vmin=min_val, vmax=-min_val)
    # elif abs(max_val)>abs(min_val):
    #     bounds = [-max_val, -min_val, 0, min_val, max_val]
    #     norm = cm.colors.Normalize(vmin=-(max_val), vmax=max_val)
    

    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    mapper = cm.ScalarMappable(cmap=cmap,norm=norm)
    print(mapper)
    bounds = [0, 0.002, 0.004, 0.006] # <-- create min and max for the colorbar scale

    color_mask = np.array([(r, g, b) for r, g, b, a in mapper.to_rgba(im_PCscores)])
    
    rgb_im = np.array(rgb_m) # <-- transform the black and white input image as a numpy array
    rgb_im[b_m] = color_mask # <-- fill the positions of the black and white numpy image where the boolean mask is TRUE with the color values
    
    print(rgb_im.shape)

    img_cropped = rgb_im[200:700, 240:1350, :]

    fig, (ax, cax) = plt.subplots(nrows=2,figsize=(12,10), gridspec_kw={"height_ratios":[1, 0.05]})
    #fig, ax = plt.subplots(figsize=(12,8))
    ax.imshow(img_cropped)
    #cax = fig.add_axes([0.001, 0.01, 0.775, 0.05]) # <-- create ax for the colorbar scale
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks = bounds, orientation='horizontal') # <-- Create a colorbar axes
    #ax.despine(left=True, bottom=True) # <-- remove x and y lines
    #fig.patch.set_visible(False)
    ax.axis('off')
    cb.outline.set_visible(False) # <-- remove colorbar outline
    ax.set(xlabel=None) # <-- remove x label
    ax.set_xticklabels([]) # <-- remove x tick labels
    ax.set_yticklabels([]) # <-- remove y tick labels
    cb.ax.tick_params(labelsize=30) # <-- change the font of the axis numbers
    ax.set(xticks=[]) # <-- remove x ticks
    ax.set(yticks=[]) # <-- remove y ticks 
    fig.subplots_adjust(hspace = 0.01, wspace = 0.01)  # <-- Add space so the colorbar doesn't overlap the plot
    
    plt.margins(0,0)
    plt.savefig(res_path+data_name+"_PC"+str(component)+"_originalrescaled.png",bbox_inches='tight') # <-- save in appropriate figure folder with region id as file title



# ----------------------------------------------------------------------------------------------

def plot_abseigen(feature, res_path, data_name):
    '''
    Save in table the absolute eigenvalues for each PC
    '''

    n_row = 15 # <-- number of rows for the table
    row_name = ['PC{}'.format(i+1) for i in range(n_row)] # <-- set row names based on number
    eigvectDF = pd.DataFrame(data=feature, index=row_name) # <-- make dataframe with absolute eigenvalues
    eigvectDF.to_csv(res_path+data_name+'_evect_abs.csv') # <-- save the dataframe in result folder



# MAIN -----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def main():
    '''
    Executes all essential commands
    '''

    bool_mask, rgb_mask, mask_blur = blur_mask(mask)

    arr_modified, list_f = modify_image(path_images, bool_mask, mask_blur, color_space, path_figures)
    save_modifiedImage(arr_modified, bool_mask, list_f, color_space, path_results, dataset)

    pca_im, new_arr = Pca(arr_modified)
    save_pcpict(new_arr, list_f, path_results, dataset)

    variances(pca_im, path_results, dataset)
    
    plot_heatmap(bool_mask, rgb_mask, pca_im, 1, color_space, path_figures, dataset)
    plot_heatmap(bool_mask, rgb_mask, pca_im, 2, color_space, path_figures, dataset)
    plot_heatmap(bool_mask, rgb_mask, pca_im, 3, color_space, path_figures, dataset)
    plot_heatmap(bool_mask, rgb_mask, pca_im, 4, color_space, path_figures, dataset)
    plot_heatmap(bool_mask, rgb_mask, pca_im, 5, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 6, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 7, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 8, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 9, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 10, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 11, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 12, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 13, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 14, color_space, path_figures, dataset)
    # plot_heatmap(bool_mask, rgb_mask, pca_im, 15, color_space, path_figures, dataset)



# GET ARGUMENTS AND EXECUTE MAIN ---------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
    print ("The script is called %s" % (sys.argv[0])) # <-- name of the script running
    arguments = len(sys.argv) - 1 # <-- count the arguments
    print ("The script is called with %i arguments" % (arguments))
    
    path_images = sys.argv[1] # <-- give path to the aligned images
    print(path_images)

    mask = cv2.imread(sys.argv[2], 0) # <-- open the mask file
    print(mask)

    color_space = sys.argv[3] # <-- give the effect to use
    print(color_space)

    path_results = sys.argv[4] # <-- path to the result folder
    print(path_results)

    path_figures = sys.argv[5] # <-- path to the figure folder
    print(path_figures)

    mask_name = sys.argv[6] # <-- name of mask
    print(mask_name)
 
    dataset = sys.argv[7] # <-- name of dataset
    print(dataset)

    main()
