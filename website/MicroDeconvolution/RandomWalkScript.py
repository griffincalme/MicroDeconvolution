import numpy as np
from numpy import linalg

import matplotlib.pyplot as plt

from skimage.exposure import rescale_intensity
from skimage.segmentation import random_walker
from skimage.color import separate_stains
from skimage.color import rgb2grey
from skimage.io import imread

import time

from pyamg import *


#Color deconvolution
#Hematoxylin(0), Red(1), DAB(2)
rgb_from_hrd = np.array([[0.65, 0.70, 0.29],
                         [0.1, 0.95, 0.95],
                         [0.27, 0.57, 0.78]])
hrd_from_rgb = linalg.inv(rgb_from_hrd)



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



def random_walk_segmentation(input_image, output_folder):

    input_image = imread(input_image)
    ihc_hrd = separate_stains(input_image, hrd_from_rgb)

    DAB_Grey_Array = stainspace_to_2d_array(ihc_hrd, 2)
    Hema_Gray_Array = stainspace_to_2d_array(ihc_hrd, 0)
    GBIred_Gray_Array = stainspace_to_2d_array(ihc_hrd, 1)

    #Perform Random Walker, fills in positive regions
    DAB_segmentation = random_walker(DAB_Grey_Array, get_markers(DAB_Grey_Array, .3, .5), beta=130, mode='cg_mg')
    Hema_segmentation = random_walker(Hema_Gray_Array, get_markers(Hema_Gray_Array, .2, .4), beta=130, mode='cg_mg')
    GBIred_segmentation = random_walker(GBIred_Gray_Array, get_markers(GBIred_Gray_Array, .4, .5), beta=130,
                                        mode='cg_mg')

    '''Compute and Output'''
    #Compute and output percentages of pixels stained by each chromagen
    pic_dimensions = np.shape(DAB_segmentation)  # both arrays same shape
    total_pixels = pic_dimensions[0] * pic_dimensions[1]

    #Change negative pixel values from 1 -> 0, positives 2 -> 1
    subtrahend_array = np.ones_like(DAB_segmentation)
    DAB_segmentation = np.subtract(DAB_segmentation, subtrahend_array)
    Hema_segmentation = np.subtract(Hema_segmentation, subtrahend_array)
    GBIred_segmentation = np.subtract(GBIred_segmentation, subtrahend_array)

    #Count positive pixels
    DAB_pixels = np.count_nonzero(DAB_segmentation)
    Hema_pixels = np.count_nonzero(Hema_segmentation)
    red_pixels = np.count_nonzero(GBIred_segmentation)


    #Percent of image covered by positive staining
    DAB_coverage_percent = (round((DAB_pixels / total_pixels * 100), 1))
    Hema_coverage_percent = (round((Hema_pixels / total_pixels * 100), 1))

    #An overlay of the DAB and Hematoxylin segmented images, for total cellular area
    total_cell_array = np.add(DAB_segmentation, Hema_segmentation)
    #Number of pixels covered by cellular area
    total_cell_pixels = np.count_nonzero(total_cell_array)

    #Percent of image covered by cellular area (DAB OR Hematoxylin)
    total_cell_percent = (round((total_cell_pixels / total_pixels * 100), 1))

    #The percentage of DAB/CD3+ cells out of the total number of cells
    percent_pos_cells = (round((DAB_pixels / total_cell_pixels * 100), 1))


    #The percentage of the image covered by cytokines
    Red_coverage_percent = (round((red_pixels / total_pixels * 100), 1))

    red_plus_total_array = np.add(total_cell_array, GBIred_segmentation)
    red_plus_total_pixels = np.count_nonzero(red_plus_total_array)

    #The percentage of the area covered by cytokines, with non-cellular regions subtracted
    adjusted_red_coverage_percent = (round((red_pixels / red_plus_total_pixels * 100), 1))




    # Plot images
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    ax0, ax1, ax2, ax3 = axes.ravel()

    ax0.imshow(input_image, cmap=plt.cm.gray, interpolation='nearest')
    ax0.set_title("Original")

    ax1.imshow(DAB_segmentation, cmap=plt.cm.gray, interpolation='nearest')
    ax1.set_title("DAB")

    ax2.imshow(GBIred_segmentation, cmap=plt.cm.gray)
    ax2.set_title("GBI red")

    ax3.imshow(Hema_segmentation, cmap=plt.cm.gray)
    ax3.set_title("Hematoxylin")

    for ax in axes.ravel():
        ax.axis('off')

    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)


    output_filename = 'output' + time.strftime("%Y-%m-%d %H:%M:%S")
    plt.savefig(output_folder + output_filename)

    #do a save csv here, maybe delete return statement after this comment

    return output_filename#, DAB_coverage_percent, Hema_coverage_percent, total_cell_percent, percent_pos_cells, Red_coverage_percent, adjusted_red_coverage_percent


#--- Test ---
#file_path = '/home/griffin/Desktop/MicroDeconvolution/TestingScripts/SamplePics/TestImage.jpg'
#save_directory = '/home/griffin/Desktop/MicroDeconvolution/website/media/images/output/'
#random_walk_segmentation(file_path, save_directory)
