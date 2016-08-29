import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt

from skimage import morphology
from skimage.io import imread
from skimage.color import separate_stains
from skimage.exposure import rescale_intensity
from skimage.color import rgb2grey
from skimage.filters import sobel


ihc_rgb = imread(r'TestImage.jpg')

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

'''DAB'''
# Rescale signals
# [:, :, (0,1,2--color)]
DAB_rescale = rescale_intensity(ihc_hrd[:, :, 2], out_range=(0, 1))
DAB_array = np.dstack((np.zeros_like(DAB_rescale), DAB_rescale, DAB_rescale))

DAB_Grey_Array = rgb2grey(DAB_array)




#Sobel transformation
elevation_map = sobel(DAB_Grey_Array)

fig, ax = plt.subplots(figsize=(4, 3))
ax.imshow(elevation_map, cmap=plt.cm.gray, interpolation='nearest')
ax.axis('off')
ax.set_title('elevation_map')


#Get markers for watershed
markers = np.zeros_like(DAB_Grey_Array)
markers[DAB_Grey_Array < .3] = 1
markers[DAB_Grey_Array > .5] = 2

fig, ax = plt.subplots(figsize=(4, 3))
ax.imshow(markers, cmap=plt.cm.spectral, interpolation='nearest')
ax.axis('off')
ax.set_title('markers')


#perform watershed, fills in positive regions
segmentation = morphology.watershed(elevation_map, markers)

fig, ax = plt.subplots(figsize=(4, 3))
ax.imshow(segmentation, cmap=plt.cm.gray, interpolation='nearest')
ax.axis('off')
ax.set_title('segmentation')


plt.show()

