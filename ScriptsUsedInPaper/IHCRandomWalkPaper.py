#Griffin Calme
#This program will output DAB, Hematoxylin, and perm red random walk segmentation

#Hematoxylin (Basic blue)= binds to nuclei
#Cytokines are perm Red Chromagen
#CD3 (T cells) are DAB

import numpy as np
from numpy import linalg

import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.exposure import rescale_intensity
from skimage.segmentation import random_walker
from skimage.color import separate_stains
from skimage.color import rgb2grey

#from pyamg import *


#Color deconvolution
#Normalized optical density matrix
#see Ruifrok AC, Johnston DA. Quantification of histological staining by color deconvolution.
# R      G      B
# X      X      X   Hematoxylin(0)
# X      X      X   Red(1)
# X      X      X   DAB(2)
#Hematoxylin(0), Red(1), DAB(2)
rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                         [0.0326, 0.873, 0.487],
                         [0.270, 0.562, 0.781]])
#conv_matrix
hrd_from_rgb = linalg.inv(rgb_from_hrd)

ihc_rgb = imread(r'TestImage.jpg')


# Rescale signals so that intensity ranges from 0 to 1
# ihc_hrd[:, :, (0,1, or 2 -- is the color channel)]
def stainspace_to_2d_array(ihc_xyz, channel):
    rescale = rescale_intensity(ihc_xyz[:, :, channel], out_range=(0,1))
    stain_array = np.dstack((np.zeros_like(rescale), rescale, rescale))
    grey_array = rgb2grey(stain_array)

    return grey_array


#Stain space conversion
ihc_hrd = separate_stains(ihc_rgb, hrd_from_rgb)


DAB_Grey_Array = stainspace_to_2d_array(ihc_hrd, 2)
Hema_Gray_Array = stainspace_to_2d_array(ihc_hrd, 0)
permRed_Gray_Array = stainspace_to_2d_array(ihc_hrd, 1)


#Get markers for random walk
def get_markers(grey_array, bottom_thresh, top_thresh):

    markers = np.zeros_like(grey_array)
    markers[grey_array < bottom_thresh] = 1
    markers[grey_array > top_thresh] = 2

    return markers


#perform Random Walker, fills in positive regions
DAB_segmentation = random_walker(DAB_Grey_Array, get_markers(DAB_Grey_Array, .3, .5), beta=130, mode='cg')
Hema_segmentation = random_walker(Hema_Gray_Array, get_markers(Hema_Gray_Array, .2, .4), beta=130, mode='cg')
permRed_segmentation = random_walker(permRed_Gray_Array, get_markers(permRed_Gray_Array, .4, .5), beta=130, mode='cg')


"""PRINTING OUTPUT"""

print('-------------------------')


'''Compute and Output'''
# Compute and output percentages of pixels stained by each chromagen
pic_dimensions = np.shape(DAB_segmentation)  # both arrays same shape
total_pixels = pic_dimensions[0] * pic_dimensions[1]


#change negative pixel values from 1 -> 0, positives 2 -> 1
subtrahend_array = np.ones_like(DAB_segmentation)
DAB_segmentation = np.subtract(DAB_segmentation, subtrahend_array)
Hema_segmentation = np.subtract(Hema_segmentation, subtrahend_array)
permRed_segmentation = np.subtract(permRed_segmentation, subtrahend_array)

#count positive pixels
DAB_pixels = np.count_nonzero(DAB_segmentation)
Hema_pixels = np.count_nonzero(Hema_segmentation)
red_pixels = np.count_nonzero(permRed_segmentation)


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

"""
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
"""

# Cytokines
print(" ")

Red_coverage_percent = (round((red_pixels / total_pixels * 100), 1))
print("The percentage of the image covered by cytokines is: " + str(Red_coverage_percent) + "%")

red_plus_total_array = np.add(total_cell_array, permRed_segmentation)
red_plus_total_pixels = np.count_nonzero(red_plus_total_array)

adjusted_red_coverage_percent = (round((red_pixels / red_plus_total_pixels * 100), 1))

print("The percentage of the area covered by cytokines, with non-cellular regions subtracted is: " + str(
    adjusted_red_coverage_percent) + "%")

print('image shape:')
print(np.shape(ihc_rgb))

"""PLOTTING IMAGES"""

#Plot images
fig, axes = plt.subplots(3, 1, figsize=(12, 11))
#ax0 = axes.ravel()
ax1,  ax2, ax3 = axes.ravel()
'''ax0,'''
#ax0.imshow(ihc_rgb, cmap=plt.cm.gray, interpolation='nearest')
#ax0.set_title("Original")

ax1.imshow(DAB_segmentation, cmap=plt.cm.gray, interpolation='nearest')
ax1.set_title("DAB")

ax2.imshow(permRed_segmentation, cmap=plt.cm.gray)
ax2.set_title("Permanent Red")

ax3.imshow(Hema_segmentation, cmap=plt.cm.gray)
ax3.set_title("Hematoxylin")

for ax in axes.ravel():
    ax.axis('off')

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)


plt.show()
