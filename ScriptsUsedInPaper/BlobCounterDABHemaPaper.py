#This program will output DAB and Hematoxylin blob counts

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
#Hematoxylin(0) perm red(1) & DAB(2)
rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                         [0.0326, 0.873, 0.487],
                         [0.270, 0.562, 0.781]])
hrd_from_rgb = linalg.inv(rgb_from_hrd)

#Import picture
ihc_rgb = data.imread(r'TestImage.jpg')

#Stain space conversion
ihc_hrd = separate_stains(ihc_rgb, hrd_from_rgb)


'''DAB'''
#Rescale signals
#[:, :, 012 color]
dab_rescale = rescale_intensity(ihc_hrd[:, :, 2], out_range=(0, 1))
dab_array = np.dstack((np.zeros_like(dab_rescale), dab_rescale, dab_rescale))


#Blob detection
image2dD = rgb2grey(dab_array)


blobs_DoG_DAB = blob_dog(image2dD, min_sigma=1, max_sigma=25, threshold=.3, overlap=0.9)
blobs_DoG_DAB[:, 2] = blobs_DoG_DAB[:, 2] * sqrt(2)

blobs = [blobs_DoG_DAB]
colors = ['orange']
titles = ['Difference of Gaussian: DAB']
sequence = zip(blobs, colors, titles)

for blobs, color, title in sequence:
    fig, ax = plt.subplots(1, 1)
    ax.set_title(title)
    ax.imshow(image2dD, interpolation='nearest', cmap=plt.cm.gray)
    for blob in blobs:
        y, x, r = blob
        c = plt.Circle((x, y), r, color=color, linewidth=1, fill=False)
        ax.add_patch(c)




'''Hematoxylin'''
#Rescale signals
#[:, :, 012 color]
hema_rescale = rescale_intensity(ihc_hrd[:, :, 0], out_range=(0, 1))
hema_array = np.dstack((np.zeros_like(hema_rescale), hema_rescale, hema_rescale))


#Blob detection
image2dH = rgb2grey(hema_array)


blobs_DoG_Hema = blob_dog(image2dH, min_sigma=10, max_sigma=15, threshold=.1, overlap=0.9)
blobs_DoG_Hema[:, 2] = blobs_DoG_Hema[:, 2] * sqrt(2)

blobs = [blobs_DoG_Hema]
colors = ['blue']
titles = ['Difference of Gaussian: Hematoxylin']
sequence = zip(blobs, colors, titles)

for blobs, color, title in sequence:
    fig, ax = plt.subplots(1, 1)
    ax.set_title(title)
    ax.imshow(image2dH, interpolation='nearest', cmap=plt.cm.gray)
    for blob in blobs:

        y, x, r = blob
        c = plt.Circle((x, y), r, color=color, linewidth=1, fill=False)
        ax.add_patch(c)



num_blobsD = len(blobs_DoG_DAB)
num_blobsH = len(blobs_DoG_Hema)
print('----------')
print('Number of DAB (CD3+) blobs detected: ' + str(num_blobsD))
print('Number of Hematoxylin (Total Nuclei) blobs detected: ' + str(num_blobsH))

PercentPos = (100 * num_blobsD/num_blobsH)//1
print('Percentage (DAB/Hematoxylin): ' + str(PercentPos)+'%')

if PercentPos >= 81:
    print('Proportion Score: 5')
elif PercentPos >= 61:
    print('Proportion Score: 4')
elif PercentPos >= 41:
    print('Proportion Score: 3')
elif PercentPos >= 21:
    print('Proportion Score: 2')
elif PercentPos >= 5:
    print('Proportion Score: 1')
elif PercentPos < 5:
    print('Proportion Score: 0')
else: print('percentage error')



plt.show()
