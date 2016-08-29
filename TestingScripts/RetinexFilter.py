#multiscale retinex color restoration?
#adapted from https://gist.github.com/shunsukeaihara/4603234

import numpy as np
from skimage.io import imread
import matplotlib.pyplot as plt


#imgArray = imread('/home/griffin/Desktop/MicroDeconvolution/TestingScripts/SamplePics/SK107 1960 B11_2_CD3_IFN G_Set 2_20X_Take 2.jpg')
#imgArray = imread(r'/home/griffin/Desktop/Screenshot at 2016-07-25 19:46:26.png')
imgArray = imread(r'/home/griffin/Desktop/house.png')
print(type(imgArray))
#plt.imshow(imgArray)

def retinex(nimg):
    #rgb array -> splits into three 2d color arrays (one for each of R, G, and B)
    nimg = nimg.transpose(2, 0, 1).astype(np.uint32)

    max_intens_green = nimg[1].max()
    nimg[0] = np.minimum(nimg[0]*(max_intens_green/float(nimg[0].max())),255)
    nimg[2] = np.minimum(nimg[2]*(max_intens_green/float(nimg[2].max())),255)
    return nimg.transpose(1, 2, 0).astype(np.uint8)

def retinex_adjust(nimg):
    """
    from 'Combining Gray World and Retinex Theory for Automatic White Balance in Digital Photography'
    """
    nimg = nimg.transpose(2, 0, 1).astype(np.uint32)
    sum_r = np.sum(nimg[0])
    sum_r2 = np.sum(nimg[0]**2)
    max_r = nimg[0].max()
    max_r2 = max_r**2
    sum_g = np.sum(nimg[1])
    max_g = nimg[1].max()
    coefficient = np.linalg.solve(np.array([[sum_r2,sum_r],[max_r2,max_r]]),
                                  np.array([sum_g,max_g]))
    nimg[0] = np.minimum((nimg[0]**2)*coefficient[0] + nimg[0]*coefficient[1],255)
    sum_b = np.sum(nimg[1])
    sum_b2 = np.sum(nimg[1]**2)
    max_b = nimg[1].max()
    max_b2 = max_r**2
    coefficient = np.linalg.solve(np.array([[sum_b2,sum_b],[max_b2,max_b]]),
                                             np.array([sum_g,max_g]))
    nimg[1] = np.minimum((nimg[1]**2)*coefficient[0] + nimg[1]*coefficient[1],255)
    return nimg.transpose(1, 2, 0).astype(np.uint8)



retinexImage = retinex(imgArray)
retinexAdjustedImage = retinex_adjust(retinexImage)
r_adj_w_out = retinex_adjust(imgArray)


"""PLOTTING IMAGES"""

#Plot images
fig, axes = plt.subplots(2, 2, figsize=(12, 11))
#ax0 = axes.ravel()
ax0, ax1,  ax2, ax3 = axes.ravel()

ax0.imshow(imgArray, cmap=plt.cm.gray, interpolation='nearest')
ax0.set_title("Original")

ax1.imshow(retinexImage, cmap=plt.cm.gray, interpolation='nearest')
ax1.set_title("retinex image")

ax2.imshow(retinexAdjustedImage, cmap=plt.cm.gray)
ax2.set_title("retinex image adjusted")

ax3.imshow(r_adj_w_out, cmap=plt.cm.gray)
ax3.set_title("retinex adjust without retinex function")

for ax in axes.ravel():
    ax.axis('off')

fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)


plt.show()