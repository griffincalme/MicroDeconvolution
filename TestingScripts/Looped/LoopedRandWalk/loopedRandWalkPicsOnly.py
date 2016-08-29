#Griffin Calme
#This program will output DAB, Hematoxylin, and GBI red images for every picture in a directory

#Hematoxylin (Basic blue)= binds to nuclei
#Cytokines are GBI Red Chromagen
#CD3 (T cells) are DAB


import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.exposure import rescale_intensity
from skimage.segmentation import random_walker
from skimage.color import separate_stains
from skimage.color import rgb2grey
from pyamg import *


import os

#Enter the master directory
PictureDirectory = r'C:\Users\Griffin\Documents\+School\UROP\LumLab\Images\AllUsablePics'


def get_filepaths(directory):
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.


# Run the above function and store its results in a variable.
full_file_paths = get_filepaths(PictureDirectory)


# Normalized optical density matrix
# Hematoxylin(0), Red(1), DAB(2)
rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                         [0.0326, 0.873, 0.487],
                         [0.270, 0.562, 0.781]])
hrd_from_rgb = linalg.inv(rgb_from_hrd)


# Rescale signals so that intensity ranges from 0 to 1
# ihc_hrd[:, :, (0,1, or 2 -- is the color channel)]
def stainspace_to_2d_array(ihc_xyz, channel):
    rescale = rescale_intensity(ihc_xyz[:, :, channel], out_range=(0,1))
    stain_array = np.dstack((np.zeros_like(rescale), rescale, rescale))
    grey_array = rgb2grey(stain_array)

    return grey_array


#Get markers for random walk
def get_markers(grey_array, bottom_thresh, top_thresh):

    markers = np.zeros_like(grey_array)
    markers[grey_array < bottom_thresh] = 1
    markers[grey_array > top_thresh] = 2

    return markers


# Main loop of the image analysis
for image_path in full_file_paths:
    if image_path.endswith(".jpg"):
        try:
            pic_file_name = os.path.basename(image_path)
            file_path = image_path

            #Import picture
            ihc_rgb = imread(file_path)


            #Stain space conversion
            ihc_hrd = separate_stains(ihc_rgb, hrd_from_rgb)


            DAB_Grey_Array = stainspace_to_2d_array(ihc_hrd, 2)
            Hema_Gray_Array = stainspace_to_2d_array(ihc_hrd, 0)
            GBIred_Gray_Array = stainspace_to_2d_array(ihc_hrd, 1)

            # perform Random Walker, fills in positive regions
            DAB_segmentation = random_walker(DAB_Grey_Array, get_markers(DAB_Grey_Array, .3, .5), beta=130, mode='cg_mg')
            Hema_segmentation = random_walker(Hema_Gray_Array, get_markers(Hema_Gray_Array, .2, .4), beta=130, mode='cg_mg')
            GBIred_segmentation = random_walker(GBIred_Gray_Array, get_markers(GBIred_Gray_Array, .4, .5), beta=130, mode='cg_mg')


            #Plot images
            fig, axes = plt.subplots(2, 2, figsize=(12, 11))
            ax0, ax1,  ax2, ax3 = axes.ravel()

            ax0.imshow(ihc_rgb, cmap=plt.cm.gray)
            ax0.set_title("Original")

            ax1.imshow(DAB_segmentation, cmap=plt.cm.gray)
            ax1.set_title("DAB")

            ax2.imshow(GBIred_segmentation, cmap=plt.cm.gray)
            ax2.set_title("GBI red")

            ax3.imshow(Hema_segmentation, cmap=plt.cm.gray)
            ax3.set_title("Hematoxylin")

            for ax in axes.ravel():
                ax.axis('off')

            fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

            SaveFileName = str(pic_file_name[:-4])
            plt.savefig(SaveFileName)

            print(pic_file_name)
            plt.close('all')


            #plt.show()

        except Exception as e:
            pic_file_name = os.path.basename(image_path)
            text_file = open("ErrorImages.txt", "w")
            text_file.write("\n error " + str(e) + " in image " + pic_file_name)
            text_file.close()

            pass