#Griffin Calme
#This program will output DAB, Hematoxylin, and GBI red

#Hematoxylin (Basic blue)= binds to nuclei
#Cytokines are GBI Red Chromagen
#CD3 (T cells) are DAB

import numpy as np
from numpy import linalg

import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.color import rgb2grey, separate_stains, combine_stains
import os

#Color deconvolution vector
#Hematoxylin(0), Red(1), DAB(2)
rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                         [0.0326, 0.873, 0.487],
                         [0.270, 0.562, 0.781]])
hrd_from_rgb = linalg.inv(rgb_from_hrd)


source_picture_directory = r'/home/griffin/Desktop/OnlyGoodImages'
#source_file_path = ('/home/griffin/Desktop/MicroDeconvolution/TestingScripts/SamplePics/SK111 818 B_2_CD3_IL 10_Set 4_20X_Take 2.jpg')
target_file_path = ('/home/griffin/Desktop/MicroDeconvolution/TestingScripts/SamplePics/SK108 3554 28_1_CD3_IL10_Set 1_20X_Take 1.jpg')



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
full_file_paths = get_filepaths(source_picture_directory)


# Rescale signals so that intensity ranges from 0 to 1
# ihc_hrd[:, :, (0,1, or 2 -- is the color channel)]
def stainspace_to_2d_array(ihc_xyz, channel):
    #rescale = rescale_intensity(ihc_xyz[:, :, channel], out_range=(0,1))
    #stain_array = np.dstack((np.zeros_like(rescale), rescale, rescale))

    #try to not reverse engineer rescale right now
    stain_array = ihc_xyz[:, :, channel]
    #plt.imshow(stain_array)
    gray_array = rgb2grey(stain_array)
    #plt.imshow(gray_array)
    return gray_array


def normalize_mapping(source_array, source_mean, target_mean):
    for x in source_array:
        x[...] = x/source_mean
        x[...] = x * target_mean
    return source_array


#Import and decovolute target
target_ihc_rgb = imread(target_file_path)
target_ihc_hrd = separate_stains(target_ihc_rgb, hrd_from_rgb)

target_Hema_Gray_Array = stainspace_to_2d_array(target_ihc_hrd, 0)
target_red_Gray_Array = stainspace_to_2d_array(target_ihc_hrd, 1)
target_DAB_Gray_Array = stainspace_to_2d_array(target_ihc_hrd, 2)


for image_path in full_file_paths:
    if image_path.endswith(".jpg"):
        try:
            pic_file_name = os.path.basename(image_path)
            file_path = image_path

            #Import picture
            source_ihc_rgb = imread(file_path)

            #Color deconvolution
            #Hematoxylin(0) & DAB(2)
            rgb_from_hrd = np.array([[0.65, 0.70, 0.29],
                                     [0.1, 0.95, 0.95],
                                     [0.27, 0.57, 0.78]])
            hrd_from_rgb = linalg.inv(rgb_from_hrd)

            #Stain space conversion
            source_ihc_hrd = separate_stains(source_ihc_rgb, hrd_from_rgb)

            source_Hema_Gray_Array = stainspace_to_2d_array(source_ihc_hrd, 0)
            source_red_Gray_Array = stainspace_to_2d_array(source_ihc_hrd, 1)
            source_DAB_Gray_Array = stainspace_to_2d_array(source_ihc_hrd, 2)

            normalized_source_Hema = normalize_mapping(source_Hema_Gray_Array, np.mean(source_Hema_Gray_Array), np.mean(target_Hema_Gray_Array))
            normalized_source_red = normalize_mapping(source_red_Gray_Array, np.mean(source_red_Gray_Array), np.mean(target_red_Gray_Array))
            normalized_source_DAB = normalize_mapping(source_DAB_Gray_Array, np.mean(source_DAB_Gray_Array), np.mean(source_DAB_Gray_Array))


            normalized_ihc_hrd = np.dstack((normalized_source_Hema, normalized_source_red, normalized_source_DAB))

            normalized_ihc_rgb = combine_stains(normalized_ihc_hrd, rgb_from_hrd)


            #Plot images
            fig, axes = plt.subplots(2, 2, figsize=(12, 11))
            #ax0 = axes.ravel()
            ax0, ax1,  ax2, ax3 = axes.ravel()

            ax0.imshow(source_ihc_rgb, cmap=plt.cm.gray, interpolation='nearest')
            ax0.set_title("Original source")

            ax1.imshow(target_ihc_rgb, cmap=plt.cm.gray, interpolation='nearest')
            ax1.set_title("original target")

            ax2.imshow(normalized_ihc_rgb, cmap=plt.cm.gray)
            ax2.set_title("normalized source")

            ax3.imshow(normalized_source_DAB, cmap=plt.cm.gray)
            ax3.set_title("normalized source DAB")

            for ax in axes.ravel():
                ax.axis('off')

            fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)


            SaveFileName = str(pic_file_name[:-4])
            plt.savefig(SaveFileName)

            print(pic_file_name)
            plt.close('all')

        except Exception as e:
            pic_file_name = os.path.basename(image_path)
            text_file = open("ErrorImages.txt", "w")
            text_file.write("\n error " + str(e) + " in image " + pic_file_name)
            text_file.close()

            pass