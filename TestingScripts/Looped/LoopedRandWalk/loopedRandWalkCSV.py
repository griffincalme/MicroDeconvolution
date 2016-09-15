#Griffin Calme
#This program will output DAB, Hematoxylin, and GBI red analysis into csv

#Hematoxylin (Basic blue)= binds to nuclei
#Cytokines are GBI Red Chromagen
#CD3 (T cells) are DAB

import numpy as np
from numpy import linalg

from skimage.io import imread
from skimage.exposure import rescale_intensity
from skimage.segmentation import random_walker
from skimage.color import separate_stains
from skimage.color import rgb2grey
from pyamg import *

import os
import csv


#Enter the master directory
PictureDirectory = r'/media/griffin/Seagate Backup Plus Drive/UROP Backup/MasterImages/UsableBrCaPics'

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

print(full_file_paths)

# CSV Table
with open('ProgOutputRW.csv', 'a', newline='') as fp:
    a = csv.writer(fp, delimiter=',')
    firstrow = [['file name', '% DAB pixels', '% Hematoxylin pixels', '% DAB & Hema', '% CD3+', 'proportion score', '% cytokine', '% cytokine (non-cellular regions subtracted)']]
    a.writerows(firstrow)

# Main loop
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


            """PRINTING OUTPUT"""

            print('-------------------------')
            print(' ')
            print('Picture directory: ' + str(file_path))
            print(' ')

            '''Compute and Output'''
            #Compute and output percentages of pixels stained by each chromagen
            pic_dimensions = np.shape(DAB_segmentation) #both arrays same shape
            total_pixels = pic_dimensions[0] * pic_dimensions[1]


            # change negative pixel values from 1 -> 0, positives 2 -> 1
            subtrahend_array = np.ones_like(DAB_segmentation)
            DAB_segmentation = np.subtract(DAB_segmentation, subtrahend_array)
            Hema_segmentation = np.subtract(Hema_segmentation, subtrahend_array)
            GBIred_segmentation = np.subtract(GBIred_segmentation, subtrahend_array)


            #Count positive pixels diaminobenzidine
            DAB_pixels = np.count_nonzero(DAB_segmentation)
            #Count positive pixels hematoxylin
            Hema_pixels = np.count_nonzero(Hema_segmentation)
            #Count positive pixels red
            red_pixels = np.count_nonzero(GBIred_segmentation)


            #Cells
            DAB_coverage_percent = (round((DAB_pixels / total_pixels * 100), 1))
            print("The percentage of the image covered by DAB is: " + str(DAB_coverage_percent) + "%")

            Hema_coverage_percent = (round((Hema_pixels / total_pixels * 100), 1))
            print("The percentage of the image covered by Hematoxylin is: " + str(Hema_coverage_percent) + "%")

            total_cell_array = np.add(DAB_segmentation, Hema_segmentation)
            total_cell_pixels = np.count_nonzero(total_cell_array)
            total_cell_percent = (round((total_cell_pixels / total_pixels * 100), 1))
            print("The percentage of the image covered by DAB & Hematoxylin is: " + str(total_cell_percent) + "%")

            percent_pos_cells = (round((DAB_pixels / total_cell_pixels * 100), 1))
            print("The percentage of CD3+ cells out of the total number of cells is: " + str(percent_pos_cells) + "%")

            PercentPos = percent_pos_cells

            if PercentPos >= 81:
                print('Proportion Score: 5')
                proportion_score = 5
            elif PercentPos >= 61:
                print('Proportion Score: 4')
                proportion_score = 4
            elif PercentPos >= 41:
                print('Proportion Score: 3')
                proportion_score = 3
            elif PercentPos >= 21:
                print('Proportion Score: 2')
                proportion_score = 2
            elif PercentPos >= 5:
                print('Proportion Score: 1')
                proportion_score = 1
            elif PercentPos >= 0:
                print('Proportion Score: 0')
                proportion_score = 0
            else:
                print('error, proportion score below zero?')
                proportion_score = -1


            #Cytokine
            print(" ")

            Red_coverage_percent = (round((red_pixels / total_pixels * 100), 1))
            print("The percentage of the image covered by cytokines is: " + str(Red_coverage_percent) + "%")

            red_plus_total_array = np.add(total_cell_array, GBIred_segmentation)
            red_plus_total_pixels = np.count_nonzero(red_plus_total_array)
            adjusted_red_coverage_percent = (round((red_pixels / red_plus_total_pixels * 100), 1))
            print("The percentage of the area covered by cytokines, with non-cellular regions subtracted is: " + str(adjusted_red_coverage_percent) + "%")


            #Append CSV with program output
            with open('ProgOutputRW.csv', 'a', newline='') as fp:
                a = csv.writer(fp, delimiter=',')
                newrow = [[str(pic_file_name), str(DAB_coverage_percent), str(Hema_coverage_percent), str(total_cell_percent), str(percent_pos_cells), str(proportion_score), str(Red_coverage_percent), str(adjusted_red_coverage_percent)]]
                a.writerows(newrow)


        except:
            with open('ProgOutputRW.csv', 'a', newline='') as fp:
                a = csv.writer(fp, delimiter=',')
                newrow = [["error in picture file"]]
                a.writerows(newrow)

            pass
