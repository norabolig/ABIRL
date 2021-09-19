# ABIRL
#
# A. C. Boley 15 Sept 2021
#
#

A Basic Image Reduction Library
===============================

ABIRL is a simple implementation of astropy and related libraries to reduce
astronomical observations. It is also a work in progress, and this repo is
to help me stay organized. 

In `python`, load the library as normal:
> >import ABIRL as ab

You can always type
> help(ab) 
to see your options

The following provides a list of the main tasks.

headerSearch
------------

`headerSearch(filelist, SEARCHKEY="OBJECT", SEARCHVAL="bias",NAMEONLY=False,LISTALL=False)`

* Scans list of file names to display exposure times and object types from HDU header.
* Will also show all values in header itels if listall=True.
* NAMEONLY just prints the filename if a match exists.

addCoords
---------

`addCoords(fout,filelist,pxu=1e-3)`

* Use to get approximate coords if you trust the plate scale. But be careful and check it.
* This is not an astrometric solution.
* pxu converts the pixelsize to the same units used for the platescale. Default assumes pixel is in um while platescale is mm/".
* fout is the prefix of the output files.

medianFilter
------------

`medianFilter(fout,filelist,BLOCKX=1,BLOCKY=1)`
    
* Applies a median filter over the entire image in that it runs a window over the original data and replaces each pixel with the median for a window centred on that pixel. 
* BLOCKY and BLOCKX control the size of the rectangular window (+/- BLOCK* from the centre pixel).
* fout is the prefix of the output files.

chunkQuartileFlatten
--------------------

`chunkQuartileFlatten(fout,filelist,BLOCKX=32,BLOCKY=32,SKIP=1)`

* Subtract the median value from all pixels within a chunk of the image.
* The median value is based on a window of width 2 BLOCKX +1 and 2 BLOCKY +1.
* The setting SKIP sets the number of cells in each image chunk, which is used to speed things up.
* SKIP=1 does not use image chunks. Larger chunks are faster, but can create pixelation.
* But this is better than the deceptive pixelation that could be created by too small of BLOCKX and BLOCKY.
* Ultimately, this is for removing gradients in the images.
* fout is the prefix of the output files.

medianFlatten
-------------

`medianFlatten(fout,filelist,BLOCKX=32,BLOCKY=32)`

* Subtracts the median within a window from each pixel.
* The window is centred on the pixel of interest, with size 2 BLOCK* + 1.
* Careful of bright features, which can create image artefacts.
* The purpose is to flatten an image.  If this is too slow, consider using chunkQuartileFlatten
* fout is the prefix of the output files.

chunkMedianSubtract
-------------------

`chunkMedianSubtract(fout,filelist,BLOCKY=16,BLOCKX=8)`

* Super fast way to flatten an image by directly subtracting medians
* based on a supergrid set by BLOCKY and BLOCKX subdivisions of the
* image in X and Y, respectively. This can leave sudden jumps between
* regions of different brightnesses. If this is the case, consider chunkQuartileFlatten
* fout is the prefix of the output files.

addBase
-------

`addBase(fout,filelist,base=100)`

* Sometimes it could be advantageous to add a flat background adu. Not standard.
* fout is the prefix of the output files.

subtractImg
-----------

`subtractImg(fout, fitsToSubtract, filelist,catch=1)`

* Subtracts one image from all images in a list. 
* fout is the prefix of the output files.

divideImg
---------

`divideImg(fout, fitsDenominator, filelist)`

* Divides all images in list by another image. 
* fout is the prefix of the output files.

meanCombine
-----------

`meanCombine(fout, filelist)`

* Mean combine files in file list takes.
* fout is the prefix of the output files.

medianCombine
-------------

`meanCombine(fout, filelist)`

* Median combine files in file list.
* fout is the prefix of the output files.

skyMedianCombine
----------------

`skyMedianCombine(fout, filelist)`

* Median combine sky flat files in file list.
* Ensures normalization before combing, as exposure times vary
* fout is the prefix of the output files.

getMedian
---------

`getMedian(filename)`

* Median of one file 

simpleAperture
--------------

`simpleAperture(fout,filelist,aperturelist)`

* You need to first create an aperture list that contains the pixel X Y coordinates, the inner radius of the aperture, and the outer radius of the
 aperture. 
* The filelist should have the same order as the aperture list.
* fout is the prefix of the output files.

getMags
-------

`getMags(fout,photlist,rad=0.0014)`

* Get magnitudes for objects in aperture file and compare with Gaia catalogue.
* Requires apiEx and assumes images are astrometrically calibrated.
* photlist is the list of aperture files.
* fout is the prefix of the output files.

regionLengths
-------------

`regionLengths(regionfile)`

* get region width in arcseconds. 
* Assumes ds9 regionfile output of a polygon. 
* Units need to be in degrees for the region file

rollImage
---------

`rollImage(fout,fin,RISE,RUN)`

* Roll the image following RISE and RUN in pixels. 
* Used for aligning satellite streaks along columns or rows. Experimental.








