#This program will output permanent red blob counts

#Hematoxylin (Basic blue)= binds to nuclei
#CD3 (T cells) are DAB

import numpy as np
from numpy import linalg
from math import sqrt

import matplotlib.pyplot as plt

from skimage import data
from skimage.color import separate_stains, rgb2grey
from skimage.exposure import rescale_intensity
from skimage.feature import blob_dog


#Color deconvolution
#DAB and perm red(1)
rgb_from_drx = np.array([[0.270, 0.562, 0.781],
                         [0.0326, 0.873, 0.487],
                         [0.0, 0.0, 0.0]])
rgb_from_drx[2, :] = np.cross(rgb_from_drx[0, :], rgb_from_drx[1, :])
drx_from_rgb = linalg.inv(rgb_from_drx)

#Import picture
ihc_rgb = data.imread(r'TestImage.jpg')


#Stain space conversion
#ihc_hax = separate_stains(ihc_rgb, hax_from_rgb)
ihc_drx = separate_stains(ihc_rgb, drx_from_rgb)


#Rescale signals? - might help, might not
#[:, :, 012 color]
permred_rescale = rescale_intensity(ihc_drx[:, :, 1], out_range=(0, 1))
permred_array = np.dstack((np.zeros_like(permred_rescale), permred_rescale, permred_rescale))


#Blob detection
image2d = rgb2grey(permred_array)


blobs_dog = blob_dog(image2d, min_sigma=1, max_sigma=40, threshold=.5, overlap=0.1)
blobs_dog[:, 2] = blobs_dog[:, 2] * sqrt(2)

blobs = [blobs_dog]
colors = ['red']
titles = ['Determinant of Gaussian: IFN-$\gamma$']
sequence = zip(blobs, colors, titles)

for blobs, color, title in sequence:
    fig, ax = plt.subplots(1, 1)
    ax.set_title(title)
    ax.imshow(image2d, interpolation='nearest', cmap=plt.cm.gray)
    for blob in blobs:
        y, x, r = blob
        c = plt.Circle((x, y), r, color=color, linewidth=1, fill=False)
        ax.add_patch(c)

num_blobs = len(blobs_dog)
print('----------')
print('Number of blobs detected: ' + str(num_blobs))


plt.show()

