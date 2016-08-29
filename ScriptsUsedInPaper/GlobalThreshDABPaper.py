#Hematoxylin (Basic blue)= binds to nuclei
#Cytokines are Permanent Red Chromogen
#CD3 (TIL) are DAB
#This program will output the DAB global threshold with intensity conservation

import numpy as np
from numpy import linalg

import matplotlib.pyplot as plt

from skimage import data
from skimage.color import separate_stains, rgb2grey
from skimage.exposure import rescale_intensity


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

#Import picture
ihc_rgb = data.imread(r'TestImage.jpg')

#Stain space conversion
ihc_hrd = separate_stains(ihc_rgb, hrd_from_rgb)


# Rescale signals so that intensity ranges from 0 to 1
# ihc_hrd[:, :, (0,1, or 2 -- is the color channel)]
DAB_rescale = rescale_intensity(ihc_hrd[:, :, 2], out_range=(0, 1))
DAB_array = np.dstack((np.zeros_like(DAB_rescale), DAB_rescale, DAB_rescale))

#set pixel values below threshold to zero
DAB_Threshold = rgb2grey(DAB_array)
DAB_Threshold[DAB_Threshold < 0.5] = 0

#Count positive pixels and divide by total pixels
DAB_pixels = np.count_nonzero(DAB_Threshold)


#Plot images
fig, axes = plt.subplots(1, 2, figsize=(12, 11))
ax0, ax1 = axes.ravel()

ax0.imshow(ihc_rgb, cmap=plt.cm.gray)
ax0.set_title("Original")

ax1.imshow(DAB_Threshold, cmap=plt.cm.gray)
ax1.set_title("DAB Global Threshold")

for ax in axes.ravel():
    ax.axis('on')

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)


plt.show()






