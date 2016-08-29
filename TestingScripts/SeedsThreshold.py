import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.color import separate_stains, rgb2grey
from skimage.exposure import rescale_intensity
from skimage.segmentation import random_walker
from skimage import morphology
from skimage.filters import sobel

#ihc_rgb = imread(r'C:\Users\griffin\Desktop\MicroDeconvolution\TestingScripts\SamplePics\SK108 3554 28_1_CD3_IL10_Set 1_20X_Take 1.jpg')
ihc_rgb = imread(r'TestImage.jpg')

# Color deconvolution
# Hematoxylin(0) & DAB(2)
rgb_from_hrd = np.array([[0.65, 0.70, 0.29],
                         [0.1, 0.95, 0.95],
                         [0.27, 0.57, 0.78]])
hrd_from_rgb = linalg.inv(rgb_from_hrd)

# Stain space conversion
ihc_hrd = separate_stains(ihc_rgb, hrd_from_rgb)


# Rescale signals so that intensity ranges from 0 to 1
# ihc_hrd[:, :, (0,1, or 2 -- is the color channel)]
def stainspace2array(ihc_xyz, channel):
    rescale = rescale_intensity(ihc_xyz[:, :, channel], out_range=(0,1))
    stain_array = np.dstack((np.zeros_like(rescale), rescale, rescale))
    grey_array = rgb2grey(stain_array)

    return grey_array


DAB_Grey_Array = stainspace2array(ihc_hrd, 2)
Hema_Grey_Array = stainspace2array(ihc_hrd, 0)
red_Grey_Array = stainspace2array(ihc_hrd, 1)

#Get markers for random walk
def get_markers(grey_array, bottom_thresh, top_thresh):

    markers = np.zeros_like(grey_array)
    markers[grey_array < bottom_thresh] = 1
    markers[grey_array > top_thresh] = 2

    return markers


DAB_markers = get_markers(DAB_Grey_Array, .3, .5)
Hema_markers = get_markers(Hema_Grey_Array, .2, .4)
red_markers = get_markers(red_Grey_Array, .4, .5)


#Plot images
fig, axes = plt.subplots(2, 2, figsize=(12, 11))
#ax0 = axes.ravel()
ax0, ax1,  ax2, ax3 = axes.ravel()

ax0.imshow(ihc_rgb, cmap=plt.cm.gray, interpolation='nearest')
ax0.set_title("Original")

ax1.imshow(DAB_markers, cmap=plt.cm.spectral, interpolation='nearest')
ax1.set_title("DAB Markers")

ax2.imshow(Hema_markers, cmap=plt.cm.spectral)
ax2.set_title("Hematoxylin Markers")

ax3.imshow(red_markers, cmap=plt.cm.spectral)
ax3.set_title("Permanent Red Markers")

for ax in axes.ravel():
    ax.axis('off')

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

plt.show()