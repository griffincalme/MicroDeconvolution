#Griffin Calme
#This program will output DAB, Hematoxylin, and perm red

#Hematoxylin (Basic blue)= binds to nuclei
#Cytokines are perm Red Chromagen
#CD3 (T cells) are DAB

import numpy as np
from numpy import linalg

import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.exposure import rescale_intensity
from skimage.color import rgb2grey #, separate_stains

from tkinter.filedialog import askopenfilename
from skimage.util import dtype
import os


#Normalized optical density matrix
#see Ruifrok AC, Johnston DA. Quantification of histological staining by color deconvolution.
# R      G      B
# X      X      X   Hematoxylin(0)
# X      X      X   Red(1)
# X      X      X   DAB(2)
#Color deconvolution
#Hematoxylin(0), Red(1), DAB(2)
rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                         [0.0326, 0.873, 0.487],
                         [0.270, 0.562, 0.781]])
#conv_matrix
hrd_from_rgb = linalg.inv(rgb_from_hrd)
print(hrd_from_rgb)


#Import picture
ihc_rgb = imread('TestImage.jpg')



#from skimage.color
def separate_stains(rgb, color_deconv_vector):
    rgb = dtype.img_as_float(rgb, force_copy=True)
    print(type(rgb))
    print(rgb[0])
    print(np.shape(rgb))
    rgb += 2
    print('-----')
    print(rgb[0])
    print('-----')
    print(-np.log(rgb[0]))
    stains = np.dot(np.reshape(-np.log(rgb), (-1, 3)), color_deconv_vector)
    print('-----')
    print()
    print(np.shape(stains))
    print(stains[0])
    return np.reshape(stains, rgb.shape)


# Rescale signals so that intensity ranges from 0 to 1
# ihc_hrd[:, :, (0,1, or 2 -- is the color channel)]
def stainspace_to_2d_array(ihc_xyz, channel):
    rescale = rescale_intensity(ihc_xyz[:, :, channel], out_range=(0,1))
    stain_array = np.dstack((np.zeros_like(rescale), rescale, rescale))
    grey_array = rgb2grey(stain_array)
    return grey_array


#Stain space conversion
ihc_hrd = separate_stains(ihc_rgb, hrd_from_rgb)

print(np.shape(ihc_hrd))
Hema_Gray_Array = stainspace_to_2d_array(ihc_hrd, 0)
permred_Gray_Array = stainspace_to_2d_array(ihc_hrd, 1)
DAB_Grey_Array = stainspace_to_2d_array(ihc_hrd, 2)



#Plot images
fig, axes = plt.subplots(2, 2, figsize=(12, 11))
#ax0 = axes.ravel()
ax0, ax1,  ax2, ax3 = axes.ravel()

ax0.imshow(ihc_rgb, cmap=plt.cm.gray, interpolation='nearest')
ax0.set_title("Original")

ax1.imshow(DAB_Grey_Array, cmap=plt.cm.gray, interpolation='nearest')
ax1.set_title("DAB")

ax2.imshow(permred_Gray_Array, cmap=plt.cm.gray)
ax2.set_title("perm red")

ax3.imshow(Hema_Gray_Array, cmap=plt.cm.gray)
ax3.set_title("Hematoxylin")

for ax in axes.ravel():
    ax.axis('off')

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)


#plt.show()
