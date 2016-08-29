import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.color import separate_stains
from skimage.exposure import rescale_intensity
from skimage.color import rgb2grey
from skimage.segmentation import random_walker


#ihc_rgb = imread(r'C:\Users\griffin\Desktop\MicroDeconvolution\TestingScripts\SamplePics\SK108 3554 28_1_CD3_IL10_Set 1_20X_Take 1.jpg')
ihc_rgb = imread(r'/home/griffin/Desktop/MicroDeconvolution/TestingScripts/SamplePics/TestImage.jpg')

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


#perform Random Walker, fills in positive regions
#can either use raw grayscale array or sobel image, not sure which is better, raw is about .2 secs faster (~ 5 secs to run each)

#segmentation = random_walker(elevation_map, markers)
segmentation = random_walker(DAB_Grey_Array, markers)


#Plot images
fig, axes = plt.subplots(2, 2, figsize=(12, 11))
#ax0 = axes.ravel()
ax0, ax1,  ax2, ax3 = axes.ravel()

ax0.imshow(ihc_rgb, cmap=plt.cm.gray, interpolation='nearest')
ax0.set_title("Original")

ax1.imshow(markers, cmap=plt.cm.spectral, interpolation='nearest')
ax1.set_title("Markers")

ax2.imshow(segmentation, cmap=plt.cm.gray)
ax2.set_title("Random Walk Segmentation")

ax3.imshow(ihc_rgb, cmap=plt.cm.gray)
ax3.set_title("original")

for ax in axes.ravel():
    ax.axis('off')

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

plt.show()


