#!/usr/bin/python
# ----------------------------------------------------------------------------------------------
# usage : 
# pythonm3 python/image_pca_heatmap.py
#          /Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/left_54off_59on/3-registred/Modalities/RGB/all/
#          /Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/smallDatasetl1/full_mask.tif
#          AB
#          /Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/images/AB_fullm_left_54off_59on/
#          /Users/fco/Desktop/PhD/1_CHAPTER1/1_GENETICS/chapter1/figures/7_gxp/AB_fullm_left_54off_59on/
#          fullm
#          54off_59on
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

import numpy as np
import cv2 
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from image_pca_heatmap import blur_mask, modify_image, Pca
import os
import sys
import PIL
from PIL import Image
import pandas as pd
from sklearn.decomposition import PCA
import functools 
import skimage.io
import skimage.viewer



# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def images_vector(im_path, bo_mask, bl_mask):
    '''
    Create a vector of pixels
    with the mask pixels, removes the background pixels
    '''

    file_list = [f for f in os.listdir(im_path) if f.endswith('.png')]
    # print(file_list)

    pixel_array = np.full((113, 273217, 3), np.NaN) # 546434 694722
    # print(pixel_array)
    # print(pixel_array.shape)

    for i, f in enumerate(file_list):
        print(i, f)
        image = cv2.imread(im_path+f)
        # print(image)
        # print(image.shape)

        img2=cv2.cvtColor(image,cv2.COLOR_BGR2RGB)
        # print(img2)
        # print(img2.shape)

        image_mask = img2.copy()
        image_mask[bo_mask == 0] = (11, 102, 35)
        # image_mask = cv2.GaussianBlur(image_mask,(5,5),cv2.BORDER_DEFAULT)
        # image_mask[bo_mask == 0] = image_mask[bo_mask == 0] / bl_mask[bo_mask == 0]
        
        # image_mask = cv2.cvtColor(image_mask, cv2.COLOR_RGB2LAB)

        # print(image_mask.shape)
        # print(image_mask[bo_mask].shape)

        Z = image_mask[bo_mask].reshape((-1,3))
        # print(Z.shape)

        Y = np.float32(Z)
        # print(Y.shape)

        pixel_array[i, :, :] = Y


    df = pixel_array.reshape((-1,3))
    # print(df.shape)

    df1 = np.float32(df)
    # print(df1.shape)

    return file_list, image_mask, pixel_array, df1



# ----------------------------------------------------------------------------------------------

def kmean_clustering(pixel_array, pixel_vector, bo_mask, gmask_im, f_list, K, fig_path):
    '''
    '''

    # Perform k_means clustering
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 10, 1.0) # <-- set conditions for kmean analysis
    ret,label,center=cv2.kmeans(pixel_vector,K,None,criteria,10,cv2.KMEANS_RANDOM_CENTERS) # <-- kmean algorithm itself
    center = np.uint8(center) # <-- list of kmean centers
    res = center[label.flatten()] 
    res2 = res.reshape((pixel_vector.shape))
    res3 = res2.reshape((pixel_array.shape))
    print(res3.shape)
    print(res3[1].shape)
    
    # Loop through file list to recover mask picture and save it as png
    for i, f in enumerate(f_list):
        print(i, f)
        gmask_im[bo_mask] = res3[i]
        im2 = Image.fromarray(gmask_im)
        im2.save(fig_path+f, format="png")

    # Loop through 
    discrete_array = np.full((len(f_list), bo_mask.sum() * 3), np.NaN) 
    for i in range(len(res3)):
        x = res3[i].ravel()
        print(x.shape)
        discrete_array[i, :] = x

    return discrete_array, center



# ----------------------------------------------------------------------------------------------

def single_color_images(fig_path, f_list, c_list, gmask_im, bo_mask, K):
    '''
    Create binary image for specific color based on the centers from the color discretisation
    Save binary images as vector for each sample, matrix image x pixels
    Save pixel count per image corresponding to the specific color cluster of interest 
    '''

    for j, center in enumerate(c_list):
        print(j, center)

        center_list = []

        pixel_array = np.full((len(f_list), bo_mask.sum() * 3), np.NaN) # 546434 694722
        # print(pixel_array)
        # print(pixel_array.shape)

        for i, f in enumerate(f_list):
            print(i, f)
            image = cv2.imread(fig_path+f) # <-- open the discrete image

            img2=cv2.cvtColor(image,cv2.COLOR_BGR2RGB) # <-- since OpenCV reads images in BGR we have to reput the image in RGB 
            
            image_mask = img2.copy()

            Z = image_mask[bo_mask].reshape((-1,3)) # <-- reshape the image to have the mask region (part of the fish to consider)
                                                    #     of the picture as a vector of pixels (3D vector)
            print(Z)
            
            count = 0
            for x, pix in enumerate(Z):
                # print(x, pix)
                
                if (functools.reduce(lambda x, y : x and y, map(lambda p, q: p == q, pix, center), True)):
                    Z[x] = pix
                    count += 1
                else:
                    Z[x] = (220, 220, 220)
    
            print(Z)
            print(Z.shape)
            print(count)

            
            gmask_im[bo_mask] = Z
            im2 = Image.fromarray(gmask_im)
            im2.save(fig_path+"k"+str(K)+"_color"+str(j)+"_"+f, format="png")

            pixel_array[i, :] = gmask_im[bo_mask].ravel()

            center_list.append(count)
        

        print(pixel_array.shape)
        n_col = bo_mask.sum() * 3
        columns_name = ['{}'.format(i+1) for i in range(n_col)]
        pixel_table = pd.DataFrame(pixel_array, columns=columns_name)
        pixel_table["images"] = f_list
        # print(pixel_table.shape)
        # print(pixel_table)
        pixel_table.to_csv(fig_path+"k"+str(K)+"_color"+str(j)+"_BinaryImage.csv")


        print(center_list)
        columns_name = ['pixel_count'] # <-- determine names of columns
        principalDF = pd.DataFrame(data=center_list, columns=columns_name) # <-- create the dataframe with images and column names
        principalDF["images"] = f_list # <-- add the sample columns corresponding to the images
        # print(principalDF.shape)
        # print(principalDF)
        principalDF.to_csv(fig_path+"k"+str(K)+"_color"+str(j)+"_counts.csv") # <-- save the dataframe to result folder with specific name


        im_pca = PCA(n_components=15)
        pixel_new = im_pca.fit_transform(pixel_array)
        
        n_col = 15 # <-- determine number of columns
        columns_name = ['PC{}'.format(i+1) for i in range(n_col)] # <-- determine names of columns based on nb of columns
        pc_DF = pd.DataFrame(data=pixel_new[:, :n_col], columns=columns_name) # <-- create the dataframe with images and column names
        pc_DF["images"] = f_list # <-- add the sample columns corresponding to the images
        # print(principalDF.shape)
        # print(principalDF)
        pc_DF.to_csv(fig_path+"k"+str(K)+"_color"+str(j)+"_PCs.csv") # <-- save the dataframe to result folder with specific name

        var_explained = im_pca.explained_variance_ratio_ # <-- gets the variance explained by PCs
        var_cum = im_pca.explained_variance_ratio_.cumsum() # <-- gets the cumulative variance explained 

        varDF = pd.DataFrame(var_explained) # <-- makes a dataframe of explained variance per PCs
        varDF.to_csv(fig_path+"k"+str(K)+"_color"+str(j)+'_var.csv') # <-- saves the dataframe of explained variance by PCs to result directory

        varcumDF = pd.DataFrame(var_cum) # <-- makes a dataframe of cumulative explained variance
        varcumDF.to_csv(fig_path+"k"+str(K)+"_color"+str(j)+'_varcum.csv') # <-- save the cumulative explained variance table to result diretory



# ----------------------------------------------------------------------------------------------

def save_pcs(array, f_list, fig_path, K):
    '''
    Save the PC table sample x PC matrix 
    '''

    n_col = 15 # <-- determine number of columns
    columns_name = ['PC{}'.format(i+1) for i in range(n_col)] # <-- determine names of columns based on nb of columns
    principalDF = pd.DataFrame(data=array[:, :n_col], columns=columns_name) # <-- create the dataframe with images and column names
    principalDF["images"] = f_list # <-- add the sample columns corresponding to the images
    # print(principalDF.shape)
    # print(principalDF)
    principalDF.to_csv(fig_path+'k'+str(K)+'_discrete_PCs.csv') # <-- save the dataframe to result folder with specific name



# ----------------------------------------------------------------------------------------------

def create_variances(pca, fig_path, K):
    '''
    Calculate variances for each PCs and cumulative variances
    Save tables of variances and cumulative variances
    '''

    var_explained = pca.explained_variance_ratio_ # <-- gets the variance explained by PCs
    var_cum = pca.explained_variance_ratio_.cumsum() # <-- gets the cumulative variance explained 

    varDF = pd.DataFrame(var_explained) # <-- makes a dataframe of explained variance per PCs
    varDF.to_csv(fig_path+'k'+str(K)+'_discrete_var.csv') # <-- saves the dataframe of explained variance by PCs to result directory

    varcumDF = pd.DataFrame(var_cum) # <-- makes a dataframe of cumulative explained variance
    varcumDF.to_csv(fig_path+'k'+str(K)+'_discrete_varcum.csv') # <-- save the cumulative explained variance table to result diretory



# ----------------------------------------------------------------------------------------------

def test_k(pixel_vector, fig_path):


    K = range(1, 10)

    list_ret = []
    nb = []

    for i in K:
        print(i)
        nb.append(i)
        
        criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 10, 1.0) # <-- set conditions for kmean analysis
        ret,label,center=cv2.kmeans(pixel_vector,i,None,criteria,10,cv2.KMEANS_RANDOM_CENTERS) # <-- kmean algorithm itself
        list_ret.append(ret)

    plt.plot(nb, list_ret)
    plt.savefig(fig_path+"elbow.png")
    plt.close()
    


# ----------------------------------------------------------------------------------------------

def get_histogram(im_path, bo_mask, f_list, fig_path):
    '''
    Create color histogram
    '''

    for i, f in enumerate(f_list):
        print(i, f)

        im = skimage.io.imread(im_path+f)
        im2 = im[bo_mask]

        colors = ("red", "green", "blue")
        channel_ids = (0, 1, 2)

        plt.xlim([0, 256])
        for channel_id, c in zip(channel_ids, colors):
            histogram, bin_edges = np.histogram(im2[:, channel_id], bins=256, range=(0, 256))
            plt.plot(bin_edges[0:-1], histogram, color=c)
        plt.xlabel("Color value")
        plt.ylabel("Pixels")
        # plt.show(block=True)
        plt.savefig(fig_path+"hist_"+f)
        plt.close()



# ----------------------------------------------------------------------------------------------

def get_histogram_global(im_vect, fig_path):
    colors = ("red", "green", "blue")
    channel_ids = (0, 1, 2)

    plt.xlim([0, 256])
    for channel_id, c in zip(channel_ids, colors):
        histogram, bin_edges = np.histogram(im_vect[:, channel_id], bins=256, range=(0, 256))
        plt.plot(bin_edges[0:-1], histogram, color=c)
    plt.xlabel("Color value")
    plt.ylabel("Pixels")
    # plt.show(block=True)
    plt.savefig(fig_path+"hist_global.png")
    plt.close()


# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def main():
    '''
    Executes all essential commands
    '''

    path = '/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/left_54off_59on/3-registred/Modalities/RGB/all/'
    path_figure = '/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/left_54off_59on/3-registred/Modalities/RGB/rgb_discrete/'

    mask = cv2.imread('/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/smallDatasetl1/body_mask.tif', 0)
    bool_mask, rgb_mask, mask_blur = blur_mask(mask)

    flist, im_green_mask, pix_arr, pix_vect = images_vector(path, bool_mask, mask_blur)
    # test_k(pix_vect, path_figure)

    cluster_array, center_list = kmean_clustering(pix_arr, pix_vect, bool_mask, im_green_mask, flist, K, path_figure)

    pca_im, new_arr = Pca(cluster_array)
    
    save_pcs(new_arr, flist, path_figure, K)

    create_variances(pca_im, path_figure, K)

    single_color_images(path_figure, flist, center_list, im_green_mask, bool_mask, K)

    get_histogram(path, bool_mask, flist, path_figure)
    get_histogram_global(pix_vect, path_figure)



# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
    print ("The script is called %s" % (sys.argv[0])) # <-- name of the script running
    arguments = len(sys.argv) - 1 # <-- count the arguments
    print ("The script is called with %i arguments" % (arguments))
    
    K = int(sys.argv[1]) # <-- number of cluster for kmeans
    print(K)
    
    main()
    




# pythonm3 -m IPython