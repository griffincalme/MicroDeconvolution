#Hematoxylin (Basic blue)= binds to nuclei
#Cytokines are GBI Red Chromagen
#CD3 (TIL) are DAB

import numpy as np
import matplotlib.pyplot as plt

from skimage import data
from skimage.color import separate_stains
from skimage.exposure import rescale_intensity
from skimage.feature import blob_dog
from skimage.color import rgb2grey
from numpy import linalg
import csv

import os

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
full_file_paths = get_filepaths(r"C:\Users\Griffin\Documents\+School\UROP\LumLab\Images\AllUsablePics")

print(full_file_paths)

# CSV Table
with open('ProgOutputThresh.csv', 'a', newline='') as fp:
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
            ihc_rgb = data.imread(file_path)

            # Normalized optical density matrix
            # Hematoxylin(0), Red(1), DAB(2)
            rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                                     [0.0326, 0.873, 0.487],
                                     [0.270, 0.562, 0.781]])
            hrd_from_rgb = linalg.inv(rgb_from_hrd)

            #Stain space conversion
            ihc_hrd = separate_stains(ihc_rgb, hrd_from_rgb)


            '''DAB'''
            #Rescale signals
            #[:, :, 012 color]
            DAB_rescale = rescale_intensity(ihc_hrd[:, :, 2], out_range=(0, 1))
            DAB_array = np.dstack((np.zeros_like(DAB_rescale), DAB_rescale, DAB_rescale))

            #set pixel values below threshold to zero
            DAB_Threshold = rgb2grey(DAB_array)
            DAB_Threshold[DAB_Threshold < 0.5] = 0

            #Count positive pixels and divide by total pixels
            DAB_pixels = np.count_nonzero(DAB_Threshold)


            '''Hematoxylin'''
            #Rescale signals
            #[:, :, 012 color]
            Hema_rescale = rescale_intensity(ihc_hrd[:, :, 0], out_range=(0, 1))
            Hema_array = np.dstack((np.zeros_like(Hema_rescale), Hema_rescale, Hema_rescale))

            #set pixel values below threshold to zero
            Hema_Threshold = rgb2grey(Hema_array)
            Hema_Threshold[Hema_Threshold < 0.4] = 0

            #Count positive pixels
            Hema_pixels = np.count_nonzero(Hema_Threshold)


            '''GBI red'''
            #Rescale signals? - might help, might not
            #[:, :, 012 color]
            GBIred_rescale = rescale_intensity(ihc_hrd[:, :, 1], out_range=(0, 1))
            GBIred_array = np.dstack((np.zeros_like(GBIred_rescale), GBIred_rescale, GBIred_rescale))

            #set pixel values below threshold to zero
            Red_Threshold = rgb2grey(GBIred_array)
            Red_Threshold[Red_Threshold < 0.475] = 0

            #Count positive pixels
            red_pixels = np.count_nonzero(Red_Threshold)


            print('-------------------------')
            print(' ')
            print('Picture directory: ' + str(file_path))
            print(' ')

            '''Compute and Output'''
            #Compute and output percentages of pixels stained by each chromagen
            pic_dimensions = np.shape(DAB_Threshold) #both arrays same shape
            total_pixels = pic_dimensions[0] * pic_dimensions[1]

            #Cells
            DAB_coverage_percent = (round((DAB_pixels / total_pixels * 100), 1))
            print("The percentage of the image covered by DAB is: " + str(DAB_coverage_percent) + "%")

            Hema_coverage_percent = (round((Hema_pixels / total_pixels * 100), 1))
            print("The percentage of the image covered by Hematoxylin is: " + str(Hema_coverage_percent) + "%")

            total_cell_array = np.add(DAB_Threshold, Hema_Threshold)
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

            red_plus_total_array = np.add(total_cell_array, Red_Threshold)
            red_plus_total_pixels = np.count_nonzero(red_plus_total_array)
            adjusted_red_coverage_percent = (round((red_pixels / red_plus_total_pixels * 100), 1))
            print("The percentage of the area covered by cytokines, with non-cellular regions subtracted is: " + str(adjusted_red_coverage_percent) + "%")


            #Append CSV with program output
            with open('ProgOutputThresh.csv', 'a', newline='') as fp:
                a = csv.writer(fp, delimiter=',')
                newrow = [[str(pic_file_name), str(DAB_coverage_percent), str(Hema_coverage_percent), str(total_cell_percent), str(percent_pos_cells), str(proportion_score), str(Red_coverage_percent), str(adjusted_red_coverage_percent)]]
                a.writerows(newrow)


        except:
            with open('ProgOutputThresh.csv', 'a', newline='') as fp:
                a = csv.writer(fp, delimiter=',')
                newrow = [["error in picture file"]]
                a.writerows(newrow)

            pass
