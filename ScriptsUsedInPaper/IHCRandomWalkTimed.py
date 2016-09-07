'''
To Do
Show speed comparison between my hand annotation and
the computer output with 6600k @ 4.5 GHz
and 6200U as well
with and without library speedup cg_mg

for human counting time for just a count and also time
for calculations of percent coverage. using pixels? or ratio of inches?
hand annotate in paint

Change title for publication
"Color Deconvolution and Random Walker Segmentation
for Automating Immunohistochemical Characterization in Breast Cancer"

Do bland-altman plot comparing mean and differences for computer vs human

If the regression of counts to patient outcome doesnt work explain why not

'''



import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.color import separate_stains, rgb2grey
from skimage.exposure import rescale_intensity
from skimage.segmentation import random_walker
from skimage import morphology
from skimage.filters import sobel

#from pyamg import *

import timeit

def timedRun():

    #ihc_rgb = imread(r'C:\Users\griffin\Desktop\MicroDeconvolution\TestingScripts\SamplePics\SK108 3554 28_1_CD3_IL10_Set 1_20X_Take 1.jpg')
    ihc_rgb = imread(r'TestImage.jpg')

    # Color deconvolution
    # Hematoxylin(0) & DAB(2)
    rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                             [0.0326, 0.873, 0.487],
                             [0.270, 0.562, 0.781]])
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


    #Get markers for random walk
    markers = np.zeros_like(DAB_Grey_Array)
    markers[DAB_Grey_Array < .3] = 1
    markers[DAB_Grey_Array > .5] = 2

    #Sobel transformation for Watershed
    elevation_map = sobel(DAB_Grey_Array)


    #perform Random Walker, fills in positive regions
    #can either use raw grayscale array or sobel image, not sure which is better, raw is about .2 secs faster (~ 5 secs to run each)

    #segmentation = random_walker(elevation_map, markers)
    RW_segmented = random_walker(DAB_Grey_Array, markers)


    #perform watershed, fills in positive regions
    watershed_segmented = morphology.watershed(DAB_Grey_Array, markers)


    #Plot images
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    #ax0 = axes.ravel()
    ax0, ax1,  ax2, ax3 = axes.ravel()

    ax0.imshow(ihc_rgb, cmap=plt.cm.gray, interpolation='nearest')
    ax0.set_title("Original")

    ax1.imshow(markers, cmap=plt.cm.spectral, interpolation='nearest')
    ax1.set_title("Segmentation Markers")

    ax2.imshow(RW_segmented, cmap=plt.cm.gray)
    ax2.set_title("Random Walk Segmentation")

    ax3.imshow(watershed_segmented, cmap=plt.cm.gray)
    ax3.set_title("Watershed Segmentation")

    for ax in axes.ravel():
        ax.axis('on')

    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

    #plt.show()

timeit.Timer(timedRun()).timeit(number=)
