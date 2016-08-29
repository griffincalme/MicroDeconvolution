#Hematoxylin (Basic blue)= binds to nuclei
#Cytokines are GBI Red Chromagen
#CD3 (TIL) are DAB

import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.color import separate_stains
from skimage.exposure import rescale_intensity
from skimage.color import rgb2grey

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
            ihc_rgb = imread(file_path)

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


            '''Hematoxylin'''
            #Rescale signals
            #[:, :, 012 color]
            Hema_rescale = rescale_intensity(ihc_hrd[:, :, 0], out_range=(0, 1))
            Hema_array = np.dstack((np.zeros_like(Hema_rescale), Hema_rescale, Hema_rescale))

            #set pixel values below threshold to zero
            Hema_Threshold = rgb2grey(Hema_array)
            Hema_Threshold[Hema_Threshold < 0.4] = 0


            '''GBI red'''
            #Rescale signals? - might help, might not
            #[:, :, 012 color]
            GBIred_rescale = rescale_intensity(ihc_hrd[:, :, 1], out_range=(0, 1))
            GBIred_array = np.dstack((np.zeros_like(GBIred_rescale), GBIred_rescale, GBIred_rescale))

            #set pixel values below threshold to zero
            Red_Threshold = rgb2grey(GBIred_array)
            Red_Threshold[Red_Threshold < 0.475] = 0


            #Compute and output percentages of pixels stained by each chromagen
            pic_dimensions = np.shape(DAB_Threshold) #both arrays same shape
            total_pixels = pic_dimensions[0] * pic_dimensions[1]


            total_cell_array = np.add(DAB_Threshold, Hema_Threshold)

            #Plot images
            fig, axes = plt.subplots(2, 2, figsize=(12, 11))
            ax0, ax1,  ax2, ax3 = axes.ravel()

            ax0.imshow(ihc_rgb, cmap=plt.cm.gray)
            ax0.set_title("Original")

            ax1.imshow(DAB_Threshold, cmap=plt.cm.gray)
            ax1.set_title("DAB")

            ax2.imshow(Red_Threshold, cmap=plt.cm.gray)
            ax2.set_title("GBI red")

            ax3.imshow(Hema_Threshold, cmap=plt.cm.gray)
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