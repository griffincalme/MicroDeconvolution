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

# Main loop of the image analysis
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
            dab_rescale = rescale_intensity(ihc_hrd[:, :, 2], out_range=(0, 1))
            dab_array = np.dstack((np.zeros_like(dab_rescale), dab_rescale, dab_rescale))

            #Blob detection
            image2d = rgb2grey(dab_array)

            blobs_DoG_DAB = blob_dog(image2d, min_sigma=1, max_sigma=25, threshold=.3, overlap=0.9)
            blobs_DoG_DAB[:, 2] = blobs_DoG_DAB[:, 2] * sqrt(2)

            blobs = [blobs_DoG_DAB]
            colors = ['red']
            titles = ['Difference of Gaussian: DAB']
            sequence = zip(blobs, colors, titles)

            for blobs, color, title in sequence:
                fig, ax = plt.subplots(1, 1)
                ax.set_title(title)
                ax.imshow(image2d, interpolation='nearest', cmap=plt.cm.gray)
                for blob in blobs:
                    y, x, r = blob
                    c = plt.Circle((x, y), r, color=color, linewidth=1, fill=False)
                    ax.add_patch(c)

            SaveFileName = str(pic_file_name[:-4]) + '_DAB'
            plt.savefig(SaveFileName)


            plt.close('all')


            '''Hematoxylin'''
            # Rescale signals
            # [:, :, 012 color]
            hema_rescale = rescale_intensity(ihc_hrd[:, :, 0], out_range=(0, 1))
            hema_array = np.dstack((np.zeros_like(hema_rescale), hema_rescale, hema_rescale))

            # Blob detection
            image2d = rgb2grey(hema_array)

            blobs_DoG_Hema = blob_dog(image2d, min_sigma=10, max_sigma=15, threshold=.1, overlap=0.9)
            blobs_DoG_Hema[:, 2] = blobs_DoG_Hema[:, 2] * sqrt(2)

            blobs = [blobs_DoG_Hema]
            colors = ['red']
            titles = ['Difference of Gaussian: Hematoxylin']
            sequence = zip(blobs, colors, titles)

            for blobs, color, title in sequence:
                fig, ax = plt.subplots(1, 1)
                ax.set_title(title)
                ax.imshow(image2d, interpolation='nearest', cmap=plt.cm.gray)
                for blob in blobs:
                    y, x, r = blob
                    c = plt.Circle((x, y), r, color=color, linewidth=1, fill=False)
                    ax.add_patch(c)

            SaveFileName = str(pic_file_name[:-4]) + '_Hema'
            plt.savefig(SaveFileName)


            plt.close('all')


            num_blobsD = len(blobs_DoG_DAB)
            num_blobsH = len(blobs_DoG_Hema)
            print('----------')
            print(pic_file_name)
            print('The picture was from directory: ' + str(file_path))
            print('Number of DAB (CD3+) blobs detected: ' + str(num_blobsD))
            print('Number of Hematoxylin (Total Nuclei) blobs detected: ' + str(num_blobsH))

            PercentPos = (100 * num_blobsD / num_blobsH) // 1
            print('Percentage (DAB/Hematoxylin): ' + str(PercentPos) + '%')

            if PercentPos < 5:
                print('Proportion Score: 0')
            elif PercentPos >= 5:
                print('Proportion Score: 1')
            elif PercentPos >= 21:
                print('Proportion Score: 2')
            elif PercentPos >= 41:
                print('Proportion Score: 3')
            elif PercentPos >= 61:
                print('Proportion Score: 4')
            elif PercentPos >= 81:
                print('Proportion Score: 5')


        except Exception as e:
            pic_file_name = os.path.basename(image_path)
            text_file = open("ErrorImages.txt", "w")
            text_file.write("\n error " + str(e) + " in image " + pic_file_name)
            text_file.close()

            pass