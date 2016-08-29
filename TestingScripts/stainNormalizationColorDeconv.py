#Griffin Calme
#This program will output DAB, Hematoxylin, and GBI red

#Hematoxylin (Basic blue)= binds to nuclei
#Cytokines are GBI Red Chromagen
#CD3 (T cells) are DAB

import numpy as np
from numpy import linalg

import matplotlib.pyplot as plt

from skimage.io import imread
from skimage.exposure import rescale_intensity
from skimage.color import rgb2grey, separate_stains, combine_stains
from skimage.util import dtype


#Color deconvolution vector
#Hematoxylin(0), Red(1), DAB(2)
rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                         [0.0326, 0.873, 0.487],
                         [0.270, 0.562, 0.781]])
hrd_from_rgb = linalg.inv(rgb_from_hrd)

source_file_path = ('/home/griffin/Desktop/MicroDeconvolution/TestingScripts/SamplePics/SK111 818 B_2_CD3_IL 10_Set 4_20X_Take 2.jpg')
target_file_path = ('/home/griffin/Desktop/MicroDeconvolution/TestingScripts/SamplePics/SK108 3554 28_1_CD3_IL10_Set 1_20X_Take 1.jpg')

#Import picture
source_ihc_rgb = imread(source_file_path)
target_ihc_rgb = imread(target_file_path)


# Rescale signals so that intensity ranges from 0 to 1
# ihc_hrd[:, :, (0,1, or 2 -- is the color channel)]
def stainspace_to_2d_array(ihc_xyz, channel):
    #rescale = rescale_intensity(ihc_xyz[:, :, channel], out_range=(0,1))
    #stain_array = np.dstack((np.zeros_like(rescale), rescale, rescale))

    #try to not reverse engineer rescale right now
    stain_array = ihc_xyz[:, :, channel]
    #plt.imshow(stain_array)
    gray_array = rgb2grey(stain_array)
    #plt.imshow(gray_array)
    return gray_array


#Stain space conversion
source_ihc_hrd = separate_stains(source_ihc_rgb, hrd_from_rgb)
target_ihc_hrd = separate_stains(target_ihc_rgb, hrd_from_rgb)


source_Hema_Gray_Array = stainspace_to_2d_array(source_ihc_hrd, 0)
source_red_Gray_Array = stainspace_to_2d_array(source_ihc_hrd, 1)
source_DAB_Gray_Array = stainspace_to_2d_array(source_ihc_hrd, 2)

target_Hema_Gray_Array = stainspace_to_2d_array(target_ihc_hrd, 0)
target_red_Gray_Array = stainspace_to_2d_array(target_ihc_hrd, 1)
target_DAB_Gray_Array = stainspace_to_2d_array(target_ihc_hrd, 2)


'''
def get_stats(input_1D_array):
    mean = np.mean(input_1D_array)
    fifth_percentile = np.percentile(input_1D_array, 5)
    ninetyfifth_percentile = np.percentile(input_1D_array, 95)
    return [mean, fifth_percentile, ninetyfifth_percentile]


source_statistical_map_vector = np.array([get_stats(source_Hema_Gray_Array),
                                          get_stats(source_red_Gray_Array),
                                          get_stats(source_DAB_Grey_Array)])

target_statistical_map_vector = np.array([get_stats(target_Hema_Gray_Array),
                                          get_stats(target_red_Gray_Array),
                                          get_stats(target_DAB_Grey_Array)])

normalized_source_map_vector = target_statistical_map_vector

print(source_statistical_map_vector)
print(source_statistical_map_vector[0,0])
print(target_statistical_map_vector)
'''


def normalize_mapping(source_array, source_mean, target_mean):

    for x in source_array:
        x[...] = x/source_mean
        x[...] = x * target_mean

        #x[...] = x/source_stat_vector[0,0]
        #x[...] = x * target_stat_vector[0,0]

    return source_array

normalized_source_Hema = normalize_mapping(source_Hema_Gray_Array, np.mean(source_Hema_Gray_Array), np.mean(target_Hema_Gray_Array))
normalized_source_red = normalize_mapping(source_red_Gray_Array, np.mean(source_red_Gray_Array), np.mean(target_red_Gray_Array))
normalized_source_DAB = normalize_mapping(source_DAB_Gray_Array, np.mean(source_DAB_Gray_Array), np.mean(source_DAB_Gray_Array))


normalized_ihc_hrd = np.dstack((normalized_source_Hema, normalized_source_red, normalized_source_DAB))

normalized_ihc_rgb = combine_stains(normalized_ihc_hrd, rgb_from_hrd)


#Plot images
fig, axes = plt.subplots(2, 2, figsize=(12, 11))
#ax0 = axes.ravel()
ax0, ax1,  ax2, ax3 = axes.ravel()

ax0.imshow(source_ihc_rgb, cmap=plt.cm.gray, interpolation='nearest')
ax0.set_title("Original source")

ax1.imshow(target_ihc_rgb, cmap=plt.cm.gray, interpolation='nearest')
ax1.set_title("original target")

ax2.imshow(normalized_ihc_rgb, cmap=plt.cm.gray)
ax2.set_title("normalized ")

ax3.imshow(normalized_source_DAB, cmap=plt.cm.gray)
ax3.set_title("normalized source DAB")

for ax in axes.ravel():
    ax.axis('off')

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)

plt.show()


#try also normalizing 5th and 95th percentile, as in "A Nonlinear Mapping Approach to Stain Normalization"