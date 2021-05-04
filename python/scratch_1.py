#!/usr/bin/python

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

# ----------------------------------------------------------------------------------------------

def blur_mask(m):
    '''
    Blur the mask
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

def modify_image(im_path, bo_mask, bl_mask, effect):
    '''
    Perform blur and color space modification on aligned images
    '''
    file_list = [f for f in os.listdir(im_path) if f.endswith('.png')]
    # print(len(file_list))
    
    if (effect == "LAB") or (effect == "RGB") or (effect == "HSV"):
        pixel_array = np.full((len(file_list), bo_mask.sum() * 3), np.NaN)
    elif (effect == "L") or (effect == "A") or (effect == "B"):
        pixel_array = np.full((len(file_list), 273217), np.NaN) #347361
    elif (effect == "AB"):
        pixel_array = np.full((len(file_list), 546434), np.NaN) #694722
    else:
        print("specify correct effect")

    # print(pixel_array.shape)

    for i, f in enumerate(file_list):
        # print(i, f)
        # print(path_images+f)
        image = cv2.imread(im_path+f)
        # print(image.shape)

        image_blurred = image.copy()
        # print(image_blurred.shape)

        image_blurred[bo_mask == 0] = 0
        # print(image_blurred.shape)

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

def save_modifiedImage(pix_arr, lof, effect):
    '''
    Save the flatten & modified images for analysis in a matrix
    '''
    if (effect == "LAB") or (effect == "RGB") or (effect == "HSV"):
        n_col = 819651 # 1042083
    elif (effect == "L") or (effect == "A") or (effect == "B"):
        n_col = 273217 # 347361
    elif (effect == "AB"):
        n_col = 546434 # 694722
    else:
        print("specify a correct effect")


    columns_name = ['{}'.format(i+1) for i in range(n_col)]
    pixel_table = pd.DataFrame(pix_arr, columns=columns_name)
    pixel_table["images"] = lof
    # print(pixel_table.shape)
    # print(pixel_table)
    pixel_table.to_csv('/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/left_54off_59on/'+effect+'_modifiedImage.csv')

# ----------------------------------------------------------------------------------------------

def Pca(pix_arr):
    '''
    Define the number of components and the new image matrix
    '''
    im_pca = PCA(n_components=15)
    # print(im_pca)
    pixel_new = im_pca.fit_transform(pix_arr)
    # print(pixel_new.shape)

    return im_pca, pixel_new

# ----------------------------------------------------------------------------------------------

def variances(pca, effect):
    '''
    Save the varaiance explained for each PCs and the cumulative variance
    '''
    var_explained = pca.explained_variance_ratio_
    var_cum = pca.explained_variance_ratio_.cumsum()
    varDF = pd.DataFrame(var_explained)
    varDF.to_csv('/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/left_54off_59on/'+effect+'_var.csv')
    varcumDF = pd.DataFrame(var_cum)
    varcumDF.to_csv('/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/left_54off_59on/'+effect+'_varcum.csv')

# ----------------------------------------------------------------------------------------------

def plot_heatmap(b_m, rgb_m, pca, component, effect):
    '''
    Create the fish heatmap corresponding to the importance of each pixels in the PCs variations
    '''
    feat_imp = pca.components_ 
    print(feat_imp)
    print(feat_imp[0])
    print(feat_imp.shape)
    count = 0 
    for i, number in enumerate(feat_imp[0]):
        if number < 0:
            count += 1
            # print(i, number)
    print(count)

    if (component == 1) and ((effect == "LAB") or (effect == "RGB") or (effect == "HSV")):
        im_PC = feat_imp[0].reshape((273217, 3), order='C') # 347361
        im_PCscores = [v for v in (im_PC[:,0] + im_PC[:,1] + im_PC[:,2])]

    elif (component == 2) and ((effect == "LAB") or (effect == "RGB") or (effect == "HSV")):
        im_PC = feat_imp[1].reshape((273217, 3), order='C')
        im_PCscores = [v for v in (im_PC[:,0] + im_PC[:,1] + im_PC[:,2])]

    elif (component == 3) and ((effect == "LAB") or (effect == "RGB") or (effect == "HSV")):
        im_PC = feat_imp[2].reshape((273217, 3), order='C')
        im_PCscores = [v for v in (im_PC[:,0] + im_PC[:,1] + im_PC[:,2])]

    elif (component == 1) and (effect == "AB"):
        im_PC = feat_imp[0].reshape((273217, 2), order='C')
        im_PCscores = [v for v in (im_PC[:,0] + im_PC[:,1])]

    elif (component == 2) and (effect == "AB"):
        im_PC = feat_imp[1].reshape((273217, 2), order='C')
        im_PCscores = [v for v in (im_PC[:,0] + im_PC[:,1])]

    elif (component == 3) and (effect == "AB"):
        im_PC = feat_imp[2].reshape((273217, 2), order='C')
        im_PCscores = [v for v in (im_PC[:,0] + im_PC[:,1])]
    
    elif (component == 1) and ((effect == "L") or (effect == "A") or (effect == "B")):
        im_PC = feat_imp[0].reshape((273217, 1), order='C')
        im_PCscores = [v for v in (im_PC[:,0])]

    elif (component == 2) and ((effect == "L") or (effect == "A") or (effect == "B")):
        im_PC = feat_imp[1].reshape((273217, 1), order='C')
        im_PCscores = [v for v in (im_PC[:,0])]
        
    elif (component == 3) and ((effect == "L") or (effect == "A") or (effect == "B")):
        im_PC = feat_imp[2].reshape((273217, 1), order='C')
        im_PCscores = [v for v in (im_PC[:,0])]
    
    else:
        print("not implemented")

    # print(im_PCscores) 
    count = 0
    for i, number in enumerate(im_PCscores):
        if number < 0:
            count += 1
            # print(i, number)  
    print(count)

    # im_PCscores = [abs(score) for score in im_PCscores] 
    print(len(im_PCscores))
    print(im_PCscores[0:10])

# =================
    # mapper = cm.ScalarMappable(cmap=cm.get_cmap('coolwarm', 10))
    # node_color = np.array([(r, g, b) for r, g, b, a in mapper.to_rgba(im_PCscores)])
    # print(node_color.shape)
    # rgb_m[b_m] = node_color
    # print(rgb_m[b_m])


    # heatmap = plt.pcolor(node_color)

    # plt.figure(figsize=(8,4))
    # plt.imshow(rgb_m) 
    # plt.colorbar(heatmap)
    # plt.show()
    # plt.savefig('/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/28.01.2021/'+effect+'_PC'+str(component)+'_heatmap.png')
# =================

    min_val, max_val = min(im_PCscores), max(im_PCscores)
    print(min_val, max_val)
    print(-max_val, -min_val)
    
    cmap = cm.coolwarm

    if abs(max_val)<abs(min_val):
        bounds = [min_val, -max_val, 0, max_val, -min_val]
        norm = cm.colors.Normalize(vmin=min_val, vmax=-min_val)
    elif abs(max_val)>abs(min_val):
        bounds = [-max_val, -min_val, 0, min_val, max_val]
        norm = cm.colors.Normalize(vmin=-(max_val), vmax=max_val)
    
    print(bounds)

    mapper = cm.ScalarMappable(cmap=cmap,norm=norm)
    # print(mapper)
    # color_list = cmap(im_PCscores)
    # print(len(color_list))
    # print(color_list[0:10])
    color_mask = np.array([(r, g, b) for r, g, b, a in mapper.to_rgba(im_PCscores)])
    print(len(color_mask))
    print(color_mask[0:10])
    rgb_m[b_m] = color_mask

    fig, ax = plt.subplots(figsize=(12,8))
    im = ax.imshow(rgb_m)
    cax = fig.add_axes([0.910, 0.100, 0.009, 0.775])
    cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks = bounds, orientation='vertical')
    cb.ax.tick_params(labelsize=5)
    # cb.ax.ticklabel_format(style='scientific',useOffset=True, useMathText=True)
    # ax.xaxis.set_major_formatter(cb.ScalarFormatter(useMathText=True))
    cb.ax.set_yticklabels(['{:e}'.format(x) for x in bounds])
    plt.show()


    return feat_imp

# ----------------------------------------------------------------------------------------------

def plot_abseigen(feature, effect):
    '''
    Save in table and plot the absolute eigenvalues for each PC
    '''
    n_row = 15
    row_name = ['PC{}'.format(i+1) for i in range(n_row)]
    eigvectDF = pd.DataFrame(data=feature, index=row_name)
    # print(eigvectDF.shape)
    # print(eigvectDF)
    eigvectDF.to_csv('/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/left_54off_59on/'+effect+'evect_abs.csv')

# ----------------------------------------------------------------------------------------------

def save_pcpict(pixels, flist, effect):

    n_col = 15
    columns_name = ['PC{}'.format(i+1) for i in range(n_col)]
    principalDF = pd.DataFrame(data=pixels[:, :n_col], columns=columns_name)
    principalDF["images"] = flist
    # print(principalDF.shape)
    # print(principalDF)
    principalDF.to_csv('/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/after_python/left_54off_59on/'+effect+'_PCs.csv')

# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def main():

    bool_mask, rgb_mask, mask_blur = blur_mask(mask)

    arr_modified, list_f = modify_image(path_images, bool_mask, mask_blur, "LAB")


    # save_modifiedImage(arr_modified, list_f, "LAB")
    
    pca_im, new_arr = Pca(arr_modified)

    variances(pca_im, "LAB")
    
    # pc1_im = plot_heatmap(bool_mask, rgb_mask, pca_im, 1, "AB")
    # pc2_im = plot_heatmap(bool_mask, rgb_mask, pca_im, 1, "AB")
    # pc3_im = plot_heatmap(bool_mask, rgb_mask, pca_im, 1, "AB")
    
    # save_pcpict(new_arr, list_f, "LAB")

# ----------------------------------------------------------------------------------------------






if __name__ == "__main__":
    # execute only if run as a script
    path_images = "/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/left_54off_59on/3-registred/Modalities/RGB/"
    mask = cv2.imread("/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/smallDatasetl1/mask1.tif", 0)
    main()