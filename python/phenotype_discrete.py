#!/usr/bin/python
# by: Floriane Coulmance: 11/05/2021
# usage : 
# pythonm3 python/phenotype_discrete.py <IMAGE_PATH> <MASK_PATH> <IMAGE_RESULT_PATH> <FIGURE_PATH> <CLUSTER_NUMBER> <COLOR_SPACE> <OUTPUT_PATH> <MASK_LABEL>
# ----------------------------------------------------------------------------------------------
# IMAGE_PATH : $BASE_DIR/ressources/images/left_54off_59on/3-registred/Modalities/RGB/all/
# MASK_PATH : $BASE_DIR/ressources/images/(body_mask.tif, full_mask.tif)
# IMAGE_RESULT_PATH : $BASE_DIR/ressources/images/left_54off_59on/3-registred/Modalities/RGB/discrete/<COLOR_SPACE>/<CLUSTER>
# OUTPUT_PATH : $BASE_DIR/images/discrete/(LAB, RGB)/
# FIGURE_PATH : $BASE_DIR/figures/7_gxp/discrete/<COLOR_SPACE>/<CLUSTER>/
# CLUSTER_NUMBER : number of clusters, 1, 2, 3, 4, 5, 6, 7
# COLOR_SPACE : RGB, LAB
# MASK_LABEL : bodym or fullm
# ----------------------------------------------------------------------------------------------



# ********** Import all necessary library *********
# -------------------------------------------------

import numpy as np
import cv2 
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, ListedColormap, LinearSegmentedColormap
from phenotype_continuous import blur_mask, Pca
import os
import sys
import errno
import PIL
from PIL import Image
import pandas as pd
from sklearn.decomposition import PCA
import functools 
import skimage.io
import skimage.viewer



# PRACTICAL FUNCTION ---------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def try_mkdir(dirname):

	try:
		os.mkdir(dirname)
	except OSError as exc:
		if exc.errno != errno.EEXIST:
			raise
		pass
	return 0



# FUNCTIONS ------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def images_vector(im_path, bo_mask, bl_mask, effect):
    '''
    Perform blur and color space modification on aligned images
    Creates a flatten modified images and store it into a vector
    Store all the flatten images into a table of flatten images
    '''

    file_list = [f for f in os.listdir(im_path) if f.endswith('.png')] # --> list all image files in the folder
    # print(file_list)

    # print(bo_mask.sum())

    # --> assuming we are working in RGB or LAB only,
    #     create an empty array of the shape of a 3 channels image
    if (effect == "LAB") or (effect == "RGB"):
        pixel_array = np.full((len(file_list), bo_mask.sum() * 3), np.NaN) # 546434 694722
    else:
        print("specify correct effect")
       
    # print(pixel_array.shape)

    # --> loop through image and modify according to parameter, 
    #      store modifications in array previously created
    for i, f in enumerate(file_list): 
        # print(i, f)
        image = cv2.imread(im_path+f)
        # print(image.shape)

        image_mask = image.copy()
        # print(image_mask.shape)

        image_mask[bo_mask == 0] = (11, 102, 35)
        # print(image_mask.shape)

        image_mask = cv2.GaussianBlur(image_mask,(5,5),cv2.BORDER_DEFAULT)
        # print(image_mask.shape)

        image_mask[bo_mask == 0] = image_mask[bo_mask == 0] / bl_mask[bo_mask == 0]
        # print(image_mask.shape)

        if (effect == "RGB"):
            im=cv2.cvtColor(image_mask,cv2.COLOR_BGR2RGB)
        elif (effect == "LAB"):
            im = cv2.cvtColor(image_mask, cv2.COLOR_RGB2LAB)
        else:
            print("specify a correct effect")

        # print(im.shape) 

        Z = im[bo_mask].reshape((-1,3))
        # print(Z.shape)

        Y = np.float32(Z)
        # print(Y.shape)

        pixel_array[i, :] = Y


    df = pixel_array.reshape((-1,3))
    # print(df.shape)

    df1 = np.float32(df)
    # print(df1.shape)

    return file_list, im, pixel_array, df1



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

def test_k(pixel_vector, fig_path):


    K = range(1, 10)

    list_ret = []
    nb = []

    for i in K:
        print(i)
        nb.append(i)
        
        criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 10, 1.0) # --> set conditions for kmean analysis
        ret,label,center=cv2.kmeans(pixel_vector,i,None,criteria,10,cv2.KMEANS_RANDOM_CENTERS) # --> kmean algorithm itself
        list_ret.append(ret)

    plt.plot(nb, list_ret)
    plt.savefig(fig_path+"elbow.png")
    plt.close()
    


# ----------------------------------------------------------------------------------------------

def kmean_clustering(pixel_array, pixel_vector, bo_mask, gmask_im, f_list, K, res_path):
    '''
    Perform Kmean clustering of colors on all images
    '''

    try_mkdir(res_path+"all/") # --> create directory to store image with all color clusters

    # --> Perform k_means clustering
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 10, 1.0) # --> set conditions for kmean analysis
    ret,label,center=cv2.kmeans(pixel_vector,K,None,criteria,10,cv2.KMEANS_RANDOM_CENTERS) # --> kmean algorithm itself
    center = np.uint8(center) # --> list of kmean centers
    res = center[label.flatten()] 
    res2 = res.reshape((pixel_vector.shape))
    res3 = res2.reshape((pixel_array.shape))
    print(res3.shape)
    print(res3[1].shape)
    
    # --> Loop through file list to recover mask picture and save it as png
    for i, f in enumerate(f_list):
        print(i, f)
        gmask_im[bo_mask] = res3[i]
        im2 = Image.fromarray(gmask_im)
        im2.save(res_path+"all/k"+str(K)+"_all_"+f, format="png")

    # --> Create clustered image vector of pixels
    discrete_array = np.full((len(f_list), bo_mask.sum() * 3), np.NaN) 
    for i in range(len(res3)):
        x = res3[i].ravel()
        print(x.shape)
        discrete_array[i, :] = x

    return discrete_array, center



# ----------------------------------------------------------------------------------------------

def save_pcs(array, f_list, res_path, m_name, K):
    '''
    Save the PC table sample x PC matrix 
    '''

    try_mkdir(res_path+"all/") # --> create directory to store image with all color clusters

    n_col = 15 # --> determine number of columns
    columns_name = ['PC{}'.format(i+1) for i in range(n_col)] # --> determine names of columns based on nb of columns
    principalDF = pd.DataFrame(data=array[:, :n_col], columns=columns_name) # --> create the dataframe with images and column names
    principalDF["images"] = f_list # --> add the sample columns corresponding to the images
    # print(principalDF.shape)
    # print(principalDF)
    principalDF.to_csv(res_path+"all/"+m_name+"_k"+str(K)+"_all_PCs.csv") # --> save the dataframe to result folder with specific name



# ----------------------------------------------------------------------------------------------

def create_variances(pca, res_path, m_name, K):
    '''
    Calculate variances for each PCs and cumulative variances
    Save tables of variances and cumulative variances
    '''

    try_mkdir(res_path+"all/") # --> create directory to store image with all color clusters

    var_explained = pca.explained_variance_ratio_ # --> gets the variance explained by PCs
    var_cum = pca.explained_variance_ratio_.cumsum() # --> gets the cumulative variance explained 

    varDF = pd.DataFrame(var_explained) # --> makes a dataframe of explained variance per PCs
    varDF.to_csv(res_path+"all/"+m_name+"_k"+str(K)+"_all_var.csv") # --> saves the dataframe of explained variance by PCs to result directory

    varcumDF = pd.DataFrame(var_cum) # --> makes a dataframe of cumulative explained variance
    varcumDF.to_csv(res_path+"all/"+m_name+"_k"+str(K)+"_allcolors_varcum.csv") # --> save the cumulative explained variance table to result diretory



# ----------------------------------------------------------------------------------------------

def single_color_images(im_path, out_path, f_list, c_list, gmask_im, bo_mask, m_name, K):
    '''
    Create binary image for specific color based on the centers from the color discretisation
    Save binary images as vector for each sample, matrix image x pixels
    Save pixel count per image corresponding to the specific color cluster of interest 
    '''

    d = {}
    for j, center in enumerate(c_list):
        print(j, center)

        try_mkdir(im_path+"color"+str(j)+"/") # --> create directory per color

        center_list = []

        pixel_array = np.full((len(f_list), bo_mask.sum() * 3), np.NaN) # 546434 694722
        # print(pixel_array)
        # print(pixel_array.shape)

        for i, f in enumerate(f_list):
            print(i, f)
            image = cv2.imread(im_path+"all/"+m_name+"_k"+str(K)+"_all_"+f) # --> open the discrete image

            img2=cv2.cvtColor(image,cv2.COLOR_BGR2RGB) # --> since OpenCV reads images in BGR we have to reput the image in RGB 
            
            image_mask = img2.copy()

            Z = image_mask[bo_mask].reshape((-1,3)) # --> reshape the image to have the mask region (part of the fish to consider)
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
            im2.save(im_path+"color"+str(j)+"/"+m_name+"_k"+str(K)+"_color"+str(j)+"_"+f, format="png")

            pixel_array[i, :] = gmask_im[bo_mask].ravel()

            center_list.append(count)
        

        print(pixel_array.shape)
        n_col = bo_mask.sum() * 3
        columns_name = ['{}'.format(i+1) for i in range(n_col)]
        pixel_table = pd.DataFrame(pixel_array, columns=columns_name)
        pixel_table["images"] = f_list
        # print(pixel_table.shape)
        # print(pixel_table)
        pixel_table.to_csv(out_path+"color"+str(j)+"/"+m_name+"_k"+str(K)+"_color"+str(j)+"_BinaryImage.csv")


        print(center_list)
        columns_name = ['pixel_count'] # --> determine names of columns
        principalDF = pd.DataFrame(data=center_list, columns=columns_name) # --> create the dataframe with images and column names
        principalDF["images"] = f_list # --> add the sample columns corresponding to the images
        # print(principalDF.shape)
        # print(principalDF)
        principalDF.to_csv(out_path+"color"+str(j)+"/"+m_name+"_k"+str(K)+"_color"+str(j)+"_counts.csv") # --> save the dataframe to result folder with specific name


        im_pca = PCA(n_components=15)
        pixel_new = im_pca.fit_transform(pixel_array)
        d["{0}".format(j)] = im_pca
        
        n_col = 15 # --> determine number of columns
        columns_name = ['PC{}'.format(i+1) for i in range(n_col)] # --> determine names of columns based on nb of columns
        pc_DF = pd.DataFrame(data=pixel_new[:, :n_col], columns=columns_name) # --> create the dataframe with images and column names
        pc_DF["images"] = f_list # --> add the sample columns corresponding to the images
        # print(principalDF.shape)
        # print(principalDF)
        pc_DF.to_csv(out_path+"color"+str(j)+"/"+m_name+"_k"+str(K)+"_color"+str(j)+"_PCs.csv") # --> save the dataframe to result folder with specific name


        var_explained = im_pca.explained_variance_ratio_ # --> gets the variance explained by PCs
        var_cum = im_pca.explained_variance_ratio_.cumsum() # --> gets the cumulative variance explained 

        varDF = pd.DataFrame(var_explained) # --> makes a dataframe of explained variance per PCs
        varDF.to_csv(out_path+"color"+str(j)+"/"+m_name+"_k"+str(K)+"_color"+str(j)+"_var.csv") # --> saves the dataframe of explained variance by PCs to result directory

        varcumDF = pd.DataFrame(var_cum) # --> makes a dataframe of cumulative explained variance
        varcumDF.to_csv(out_path+"color"+str(j)+"/"+m_name+"_k"+str(K)+"_color"+str(j)+"_varcum.csv") # --> save the cumulative explained variance table to result diretory

    return d



# ----------------------------------------------------------------------------------------------

def plot_heatmap(b_m, rgb_m, d_pca, c_list, component, fig_path, m_name):
    '''
    Create the fish heatmap corresponding to the importance of each pixels in the PCs variations
    '''

    for j, center in enumerate(c_list):
        print(j, center)
        
        try_mkdir(fig_path+"k"+str(K)+"_color"+str(j)+"/") # --> create directory per color

        print(d_pca["{0}".format(j)]) 

        j_pca = d_pca["{0}".format(j)] # --> retrieve the PCA corresponding to the color

        feat_imp = j_pca.components_ # --> get the eigenvector from PCA for each PCs
        print(feat_imp)

        pc_index = component - 1 # --> get the index in table corresponding the PC specified in parameters
        print(pc_index)


        im_PC = feat_imp[pc_index].reshape((b_m.sum(), 3), order='C') # --> reshape eigenvector of PC of interest


        im_PCscores = [v for v in np.linalg.norm(im_PC,axis = 1)] # --> get score/vector from norm of the 3 channels of the eigenvector


        min_val, max_val = min(im_PCscores), max(im_PCscores) # --> get max and min values
        print(min_val, max_val)
    
        # Determine color norm, color palette, color map and bounds for heatmap
        norm = matplotlib.colors.Normalize(0,0.015)
        colors = [[norm(0), "grey"],
              [norm(0.005), "yellow"],
            [norm(0.01), "orange"],
            [norm( 0.015), "darkred"]]

        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
        mapper = cm.ScalarMappable(cmap=cmap,norm=norm)
        bounds = [0, 0.005, 0.010, 0.015] # --> create min and max for the colorbar scale

        # Translate eigenvector norm scores to RGB colors 
        color_mask = np.array([(r, g, b) for r, g, b, a in mapper.to_rgba(im_PCscores)])

        rgb_m[b_m == False] = (255, 255, 255)
        rgb_im = np.array(rgb_m) # --> transform the black and white input image as a numpy array
        rgb_im[b_m] = color_mask # --> fill the positions of the black and white numpy image where the boolean mask is TRUE with the color values

        # plot the heatmap
        fig, ax = plt.subplots(figsize=(12,8))
        im = ax.imshow(rgb_im)
        cax = fig.add_axes([0.910, 0.100, 0.009, 0.775]) # --> create ax for the colorbar scale
        cb = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks = bounds, orientation='vertical') # --> Create a colorbar axes
        #ax.despine(left=True, bottom=True) # --> remove x and y lines
        #fig.patch.set_visible(False)
        ax.axis('off')
        ax.set(xlabel=None) # --> remove x label
        ax.set_xticklabels([]) # --> remove x tick labels
        ax.set_yticklabels([]) # --> remove y tick labels
        ax.set(xticks=[]) # --> remove x ticks
        ax.set(yticks=[]) # --> remove y ticks 
        fig.subplots_adjust(right=.9)  # --> Add space so the colorbar doesn't overlap the plot

        plt.savefig(fig_path+"color"+str(j)+"/"+m_name+"_k"+str(K)+"_color"+str(j)+"_PC"+str(component)+".png") # --> save in appropriate figure folder with region id as file title



# MAIN -----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

def main():
    '''
    Executes all essential commands
    '''

    # --> create rgb mask and boolean mask for image analysis
    bool_mask, rgb_mask, mask_blur = blur_mask(mask)

    # --> create a global vector of pixels (3D vector) encompassing all images within the mask area
    flist, im_green_mask, pix_arr, pix_vect = images_vector(path_images, bool_mask, mask_blur, color_space)
    
    # --> This part is exploratory and will help us decide the number of color clusters we want
    # get_histogram(path_images, bool_mask, flist, path_figure)
    # get_histogram_global(pix_vect, path_figure)
    # test_k(pix_vect, path_figures)

    # --> This part creates discrete color images for all the fish,
    # --> performs PCA on all the colors of the discrete images,
    # --> save PCA results and variances,
    # --> save discrete color images with all the colors
    cluster_array, center_list = kmean_clustering(pix_arr, pix_vect, bool_mask, im_green_mask, flist, K, imout_path)
    pca_im, new_arr = Pca(cluster_array) 
    save_pcs(new_arr, flist, output_path, mask_name, K)
    create_variances(pca_im, output_path, mask_name, K)

    # --> Create images, pca, variances and pixel counts for each color from kmean seperately
    pca_dict = single_color_images(imout_path, output_path, flist, center_list, im_green_mask, bool_mask, mask_name, K)

    # --> Create a pixel heatmap of pixel contribution to the first 15 Principal Components
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 1, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 2, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 3, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 4, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 5, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 6, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 7, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 8, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 9, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 10, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 11, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 12, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 13, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 14, path_figures, mask_name)
    plot_heatmap(bool_mask, rgb_mask, pca_dict, center_list, 15, path_figures, mask_name)

 

# GET ARGUMENTS AND EXECUTE MAIN ---------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

if __name__ == "__main__":
    
    print ("The script is called %s" % (sys.argv[0])) # --> name of the script running
    arguments = len(sys.argv) - 1 # --> count the arguments
    print ("The script is called with %i arguments" % (arguments))
    
    # path_images = "/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/left_54off_59on/3-registred/Modalities/RGB/"
    path_images = sys.argv[1] # --> give path to the aligned images
    print(path_images)

    # mask = cv2.imread('/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/smallDatasetl1/body_mask.tif', 0)
    mask = cv2.imread(sys.argv[2], 0) # --> open the mask file
    print(mask)

    imout_path = sys.argv[3] # --> path to the result folder
    print(imout_path)

    output_path = sys.argv[4] # --> get the path for output results
    print(output_path)

    # path_figures = '/Users/fco/Desktop/PhD/1_CHAPTER1/0_IMAGES/convert_png/left_54off_59on/3-registred/Modalities/RGB/rgb_discrete/'
    path_figures = sys.argv[5] # --> path to the figure folder
    print(path_figures)

    K = int(sys.argv[6]) # --> number of cluster for kmeans
    print(K)

    color_space = sys.argv[7] # --> give the color space to use
    print(color_space)

    mask_name = sys.argv[8] # <-- name of mask
    print(mask_name)

    main()
    


# pythonm3 -m IPython