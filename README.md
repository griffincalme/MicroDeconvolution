# MicroDeconvolution

Tumor microenvironment color deconvolution and segmentation.

Quantifies staining area for cells, cytokines, etc. on microscope slide images.
Can differentiate up to 3 different colors.

#ScriptsUsedInPaper
These are for replicating results in the paper.

###BlobCounterDABHemaPaper.py 
Outputs Determinant of Gaussian blob detection for DAB and Hematoxylin

###BlobCounterRedPaper.py
Same thing as above, but for permanent red-stained cytokines in the extracellular matrix.
This is an example of something NOT to do, as the cytokine staining (IL-10 and IFN-gamma in this paper)
is often diffuse and amorphic.

###ColorDeconvolutionPaper.py
Color deconvolution for DAB, hematoxylin and permanent red, the basis of all other scripts.

###GlobalThreshDABPaper.py
An example of a global threshold for the DAB stain intensity matrix
first unmixes stain channels, then applies threshold at T > 0.5

###IHCRandomWalkPaper.py
For producing part of the flowchart

###IHCRandomWalkResults.py
The paper's proposed algorithm for separating, segmenting, and quantifying staining.
Performs color deconvolution and then a random walk using seeds that are
determined by thresholds on the stain intensity matrices.

###RW_Watershed_Comparison_Paper.py
Compares the random walker to watershed segmentation. Illustrates the author's choice
for the random walker segmentation algorithm.

###TestImage.jpg
A cropped de-identified test image for use in these scripts.


#TestingScripts
##Where I have been developing new ideas.
###calculate_P_hat_from_256RGB.py
This script requires the RGB values (0-255) for a stain
This script outputs decimal RGB values for the stain,
optical density (OD) of the stain,
and most importantly, the normilized OD.

The normalized OD is one row of the color deconvolution matrix as seen here:

```python
rgb_from_hrd = np.array([[0.644, 0.710, 0.285],
                         [0.0326, 0.873, 0.487],
                         [0.270, 0.562, 0.781]])
```

This example outputs the values for permanent red.


###Looped
These scripts iterate over a folder of IHC images and either output processed images or a csv of the data.


#Website
A website that will allow users to test out the algorithms used in the 
paper in a user-friendly manner. Currently it doesn't do much and is in 
need of further development.
