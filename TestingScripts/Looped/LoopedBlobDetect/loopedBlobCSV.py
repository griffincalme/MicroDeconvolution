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
from math import sqrt

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
with open('ProgOutputBlobDetect.csv', 'a', newline='') as fp:
    a = csv.writer(fp, delimiter=',')
    firstrow = [['file name', '# DAB blobs', '# Hema blobs', '% DAB/Hema', 'Score']]
    a.writerows(firstrow)

# Main loop
for image_path in full_file_paths:
    if image_path.endswith(".jpg"):
        try:
            pic_file_name = os.path.basename(image_path)
            file_path = image_path

            # Import picture
            ihc_rgb = data.imread(file_path)

            # Normalized optical density matrix
            # Hematoxylin(0), Red(1), DAB(2)
            rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                                     [0.0326, 0.873, 0.487],
                                     [0.270, 0.562, 0.781]])
            # conv_matrix
            hrd_from_rgb = linalg.inv(rgb_from_hrd)
            print(hrd_from_rgb)

            hrd_from_rgb = linalg.inv(rgb_from_hrd)

            # Stain space conversion
            ihc_hrd = separate_stains(ihc_rgb, hrd_from_rgb)

            '''DAB'''
            # Rescale signals
            # [:, :, 012 color]
            dab_rescale = rescale_intensity(ihc_hrd[:, :, 2], out_range=(0, 1))
            dab_array = np.dstack((np.zeros_like(dab_rescale), dab_rescale, dab_rescale))

            # Blob detection
            image2d = rgb2grey(dab_array)

            blobs_DoG_DAB = blob_dog(image2d, min_sigma=1, max_sigma=25, threshold=.3, overlap=0.9)
            blobs_DoG_DAB[:, 2] = blobs_DoG_DAB[:, 2] * sqrt(2)


            '''Hematoxylin'''
            # Rescale signals
            # [:, :, 012 color]
            hema_rescale = rescale_intensity(ihc_hrd[:, :, 0], out_range=(0, 1))
            hema_array = np.dstack((np.zeros_like(hema_rescale), hema_rescale, hema_rescale))

            # Blob detection
            image2d = rgb2grey(hema_array)

            blobs_DoG_Hema = blob_dog(image2d, min_sigma=10, max_sigma=15, threshold=.1, overlap=0.9)
            blobs_DoG_Hema[:, 2] = blobs_DoG_Hema[:, 2] * sqrt(2)


            num_blobsD = len(blobs_DoG_DAB)
            num_blobsH = len(blobs_DoG_Hema)
            print('----------')
            print('The picture was from directory: ' + str(file_path))
            print('Number of DAB (CD3+) blobs detected: ' + str(num_blobsD))
            print('Number of Hematoxylin (Total Nuclei) blobs detected: ' + str(num_blobsH))

            PercentPos = (100 * num_blobsD / num_blobsH) // 1
            print('Percentage (DAB/Hematoxylin): ' + str(PercentPos) + '%')

            if PercentPos >= 81:
                proportionScore = 5
            elif PercentPos >= 61:
                proportionScore = 4
            elif PercentPos >= 41:
                proportionScore = 3
            elif PercentPos >= 21:
                proportionScore = 2
            elif PercentPos >= 5:
                proportionScore = 1
            elif PercentPos < 5:
                proportionScore = 0
            else: proportionScore = 'percentage error'


            print(proportionScore)

            #Append CSV with program output
            with open('ProgOutputBlobDetect.csv', 'a', newline='') as fp:
                a = csv.writer(fp, delimiter=',')
                newrow = [
                    [str(pic_file_name), str(num_blobsD), str(num_blobsH), str(PercentPos), str(proportionScore)]]
                a.writerows(newrow)


        except:
            with open('ProgOutputBlobDetect.csv', 'a', newline='') as fp:
                a = csv.writer(fp, delimiter=',')
                newrow = [["error in picture file"]]
                a.writerows(newrow)

            pass

