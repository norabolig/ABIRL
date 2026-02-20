# 15 Sept 2021
#
# A Basic Image Reduction Library using astropy and related tools. 
# Has also been used as a base for processing for NEOSSAT data.
#
#
#
from astropy.io import fits as pyfits
import sys
import numpy as np
import os

def set_bitpix(hdul,BITPIX):
    if BITPIX==16: hdul[0].scale('int16')
    elif BITPIX==32: hdul[0].scale('int32')
    elif BITPIX==-32: hdul[0].scale('float32')
    elif BITPIX==-64: hdul[0].scale('float64')
    else:
        print("Invalid selection for BITPIX {}".format(BITPIX))

def headerSearch(filelist, SEARCHKEY="OBJECT", SEARCHVAL="bias",NAMEONLY=False,LISTALL=False):
    """
    Scans list of file names to display exposure times and object types from HDU header.
    Will also show all values in header itels if listall=True
    """
    fh = open(filelist, "r")

    for filename in fh:
        hdul = pyfits.open(filename.rstrip())  # rstrip is to remove carriage return
        head = hdul[0].header
#        print(filename,head[SEARCHKEY])

        def printheader(s):
            if NAMEONLY:
               print(filename.strip()) # just used to get file names that match
            else:
               print(s)

        if LISTALL:
            for item in head: print("{} = {}".format(item,head[item]))
        else:
            try:
                if head[SEARCHKEY].find(SEARCHVAL)>-1:
                    printheader(" file = {}, OBJECT = {}, EXPOSURE = {}".format(filename.rstrip(), head["OBJECT"],
                                                                          head["EXPTIME"]))
                else: 
                    print("The SEARCHVAL {} was not found in file {}.".format(SEARCHVAL,filename.rstrip()))
                    print("Here is what I found: {} = {}".format(SEARCHKEY,head[SEARCHKEY]))
            except:
                print("The SEARCHKEY {} was probably not found. Perhaps list header to see what is there.".format(SEARCHKEY))
                print("If the SEARCHKEY is there, then check the defaul print statement,")
                print("as some of your header keywords might be different from what is expected ('OBJECT' and 'EXPTIME')")

        hdul.close()

    fh.close()

def addCoords(fout,filelist,pxu=1e-3,BITPIX=-32):
    """IF YOU KNOW THE APPROXIMATE PLATE SCALE AND COORDS, YOU CAN USE THIS. BUT BE CAREFUL. IT MUST BE CHECKED MANUALLY, ALTHOUGH IT SHOULD BE OK FOR DAO IMAGES""" 
    fh = open(filelist, "r")
    print(fh)
    iter = 0
    for line in fh:
        fpath = line.rstrip()
        hdul = pyfits.open(fpath)

        head = hdul[0].header

        pixsize = head["PIXSIZE"]*pxu # put pixel units to be compatible with pltscal
        pltscal = head["PLTSCALE"]

        pixang = pltscal*pixsize/3600 # pixang in deg

# note how some of the axes are flipped around.
        data = hdul[0].data
        hdul[0].data=np.rot90(data)
        NAX1 = head['NAXIS1']
        NAX2 = head['NAXIS2']
        print(NAX1,NAX2)
        head['NAXIS1']=(NAX2,'E-W')
        head['NAXIS2']=(NAX1,'N-S')
        print(head['NAXIS1'],head['NAXIS2'])

        head.set('CDELT1',-pixang*head['XBIN'],'deg per x pixel')
        head.set('CDELT2',pixang*head['YBIN'],'deg per y pixel')
        head.set('CTYPE1','RA---TAN')
        head.set('CTYPE2','DEC--TAN')
        head.set('CUNIT1','deg')
        head.set('CUNIT2','deg')
        head.set('CRPIX1',head['NAXIS1']//2+1,'field centre')
        head.set('CRPIX2',head['NAXIS2']//2+1,'field centre')

        ra = head['RA'].split(":")
        dec = head['DEC'].split(":")

        ra[0]=float(ra[0])
        ra[1]=float(ra[1])
        ra[2]=float(ra[2])

        dec[0]=float(dec[0])
        dec[1]=float(dec[1])
        dec[2]=float(dec[2])
        

        radeg = 360/24*(ra[0] + ra[1]/60 + ra[2]/3600)

        sgn=1.
        if dec[0] < 0: sgn = -1.
        decdeg = dec[0] + sgn*(dec[1]/60 + dec[2]/3600)

        head.set('CRVAL1',radeg,'RA at equinox (deg)')
        head.set('CRVAL2',decdeg,'DEC at equinox (deg)')

        head.set('comment', 'Set coords from file {}'.format(fpath))


        set_bitpix(hdul,BITPIX)
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

def medianFilter(fout,filelist,BLOCKX=1,BLOCKY=1,BITPIX=-32):
    """ 
    Applies a median filter over the entire image in that it runs a window over
    the original data and replaces each pixel with the median for a window centred
    on that pixel. BLOCKY and BLOCKX control the size of the rectangular window
    (+/- BLOCK* from the centre pixel).
    """

    from scipy.ndimage import median_filter
    fh = open(filelist, "r")
    iter = 0
    for line in fh:
        fpath = line.rstrip()
        print(fpath)
        hdul = pyfits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data
        print(imgdata.dtype.name)

        Ndata = imgdata.size

        NY, NX = imgdata.shape

        corrected = median_filter(imgdata,size=(BLOCKY*2+1,BLOCKX*2+1))
         

        set_bitpix(hdul,BITPIX)
        head.set('comment', 'Median filtered (despiked) from {}'.format(fpath))
        hdul[0].data = corrected
        

        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1


def chunkQuartileFlatten(fout,filelist,BLOCKX=32,BLOCKY=32,SKIP=1,BITPIX=-32):
    """ 
    Subtract the median value from all pixels within a chunk of the image.
    The median value is based on a window of width 2 BLOCKX +1 and 2 BLOCKY +1. 
    The setting SKIP sets the number of cells in each image chunk, which is used to speed things up. 
    SKIP=1 does not use image chunks. Larger chunks are faster, but can create pixelation. 
    But this is better than the deceptive pixelation that could be created by too small of BLOCKX and BLOCKY. 
     Ultimately, this is for removing gradients in the images.
    """

    from scipy.ndimage import median_filter
    fh = open(filelist, "r")
    iter = 0
    for line in fh:
        fpath = line.rstrip()
        print(fpath)
        hdul = pyfits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data
        print(imgdata.dtype.name)

        Ndata = imgdata.size

        NY, NX = imgdata.shape

        corrected = np.zeros([NY,NX])

        for j in range(SKIP,NY,SKIP*2):
          print("Working on column {}".format(j))
          jstart = j-BLOCKY
          jend = j+BLOCKY
          if jstart < 0: 
            jend -=jstart
            jstart = 0
          if jend > NY-1:
            jstart-=jend-NY+1
            jend = NY-1
          for i in range(SKIP,NX,SKIP*2):
             istart = i-BLOCKX
             iend = i+BLOCKX
             if istart < 0:
               iend -=istart
               istart = 0
             if iend > NX-1:
               istart-=iend-NX+1
               iend = NX-1

             #iq0 = int ( (BLOCKY+1)*(BLOCKX+1) * 1 / 4)
             #iq1 = int ( (BLOCKY+1)*(BLOCKX+1) * 3 / 4)
             #vsort=np.sort(np.ravel(imgdata[jstart:jend+1,istart:iend+1]))
             #med = np.median(vsort[iq0:iq1])
             med = np.median(imgdata[jstart:jend+1,istart:iend+1])
             corrected[j-SKIP:j+SKIP+1,i-SKIP:i+SKIP+1]=imgdata[j-SKIP:j+SKIP+1,i-SKIP:i+SKIP+1]-med


        head.set('comment', 'Quartile flattened from {}'.format(fpath))
        set_bitpix(hdul,BITPIX)
        hdul[0].data = corrected
        

        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

 
def medianFlatten(fout,filelist,BLOCKX=32,BLOCKY=32,BITPIX=-32):
    """ 
    Subtracts the median within a window from each pixel. 
    The window is centred on the pixel of interest, with size 2 BLOCK* + 1. 
    Careful of bright features, which can create image artefacts. 
    The purpose is to flatten an image.  If this is too slow, consider using chunkQuartileFlatten
    """
    from scipy.ndimage import median_filter
    fh = open(filelist, "r")
    iter = 0
    for line in fh:
        fpath = line.rstrip()
        print(fpath)
        hdul = pyfits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data
        print(imgdata.dtype.name)

        Ndata = imgdata.size

        NY, NX = imgdata.shape

        corrected = np.zeros([NY,NX])

        for j in range(NY):
          jstart = j-BLOCKY
          jend = j+BLOCKY
          if jstart < 0: 
            jend -=jstart
            jstart = 0
          if jend > NY-1:
            jstart-=jend-NY+1
            jend = NY-1
          for i in range(NX):
             istart = i-BLOCKX
             iend = i+BLOCKX
             if istart < 0:
               iend -=istart
               istart = 0
             if iend > NX-1:
               istart-=iend-NX+1
               iend = NX-1

             med = np.median(imgdata[jstart:jend+1,istart:iend+1])
             corrected[j][i]=imgdata[j][i]-med


        head.set('comment', 'Median flattened from {}'.format(fpath))
        set_bitpix(hdul,BITPIX)
        hdul[0].data = corrected
        

        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

   


def chunkMedianSubtract(fout,filelist,BLOCKY=16,BLOCKX=8,BITPIX=-32):
    """ 
    Super fast way to flatten an image by directly subtracting medians 
    based on a supergrid set by BLOCKY and BLOCKX subdivisions of the 
    image in X and Y, respectively. This can leave sudden jumps between 
    regions of different brightnesses. If this is the case, consider chunkQuartileFlatten
    """

    fh = open(filelist, "r")
    iter = 0
    for line in fh:
        fpath = line.rstrip()
        hdul = pyfits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data
        print(imgdata.dtype.name)

        Ndata = imgdata.size

        NY, NX = imgdata.shape

        CHUNKY = NY//BLOCKY
        CHUNKX = NX//BLOCKX

        print(NX,NY)
        print(CHUNKX,CHUNKY) 

        corrected = np.zeros((NY,NX))

        medArray = np.zeros((NY,NX))
        centres=[]

        for j in range(BLOCKY):
            jin = CHUNKY*j
            jout = min(CHUNKY*(j+1),NY)
            centres.append([])
            for i in range(BLOCKX):
               iin = CHUNKX*i
               iout = min(CHUNKX*(i+1),NX)

               medArray[j][i] = np.percentile(imgdata[jin:jout+1,iin:iout+1],50)
               centres[j].append( [ (jout+jin)/2, (iout+iin)/2 ] )

               corrected[jin:jout+1,iin:iout+1] = imgdata[jin:jout+1,iin:iout+1]-medArray[j][i]

        print(corrected.dtype.name)
        head.set('comment', 'Subtracted median chunks from {}'.format(fpath))
        hdul[0].data = corrected


        set_bitpix(hdul,BITPIX)
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1



def addBase(fout,filelist,base=100,BITPIX=-32):
    """Sometimes it could be advantageous to add a flat background adu. Not standard"""

    fh = open(filelist, "r")
    iter = 0
    for line in fh:
        fpath = line.rstrip()
        hdul = pyfits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data

        Ndata = imgdata.size

        NY, NX = imgdata.shape

        corrected = np.zeros((NY,NX))
        for j in range(NY):
            for i in range(NX):
                  corrected[j][i] = imgdata[j][i] + base

        head.set('comment', 'Added {} to file {}'.format(base, fpath))
        hdul[0].data = corrected

        set_bitpix(hdul,BITPIX)
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1


def subtractImg(fout, fitsToSubtract, filelist,catch=1,BITPIX=-32):
    """Subtracts one image from all images in a list. takes(fout,fitsToSubtract,filelist)"""

    hdulref = pyfits.open(fitsToSubtract)
    head = hdulref[0].header
    refdata = hdulref[0].data*1.
    print(refdata)
    hdulref.close()

    fh = open(filelist, "r")
    iter = 0
    for line in fh:
        fpath = line.rstrip()
        hdul = pyfits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data

        Ndata = imgdata.size

        Nref = refdata.size

        if not Ndata == Nref:
            print("Mismatch between data {} and reference {}".format(Ndata, Nref))
            sys.exit()

        NY, NX = imgdata.shape

        corrected = np.zeros((NY,NX)) # defining because I was getting some weird errors on occasion. 
        corrected = imgdata - refdata

        head.set('comment', 'Subtracted {} from file {}'.format(fitsToSubtract, fpath))
        hdul[0].data = corrected*1.

        set_bitpix(hdul,BITPIX)
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

def divideImg(fout, fitsDenominator, filelist,BITPIX=-32):
    """Divides all images in list by another image. takes(fout,fitsDenominator,filelist)"""

    hdulref = pyfits.open(fitsDenominator)
    head = hdulref[0].header
    refdata = hdulref[0].data
    hdulref.close()

    fh = open(filelist, "r")
    iter = 0
    for line in fh:
        fpath = line.rstrip()
        hdul = pyfits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data

        Ndata = imgdata.size

        Nref = refdata.size

        if not Ndata == Nref:
            print("Mismatch between data {} and reference {}".format(Ndata, Nref))
            sys.exit()

        NY, NX = imgdata.shape

        corrected = np.zeros((NY,NX))
        flags = refdata>0
        corrected[flags]=imgdata[flags]/refdata[flags]

        head.set('comment', 'Divided by {} into file {}'.format(fitsDenominator, fpath))
        hdul[0].data = corrected

        set_bitpix(hdul,BITPIX)
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

def meanCombine(fout, filelist, BITPIX=-32):
    """Mean combine files in file list takes (fout,filelist)"""

    fh = open(filelist, "r")
    NFILES = 0
    for line in fh: NFILES += 1
    fh.seek(0)
    print("Number of files is {}".format(NFILES))

    hdularr = []

    hdul = pyfits.open(fh.readline().rstrip())
    hdularr.append(hdul)

    imgdata0 = hdul[0].data
    NY, NX = imgdata0.shape
    print(NY, NX)

    imgavg = np.zeros([NY,NX])

    i = 1
    names = []
    for line in fh:
        names.append(line.rstrip())
        hdul = pyfits.open(line.rstrip())
        hdularr.append(hdul)
        data = hdul[0].data
        print(data.shape)
        imgavg[0:][0:] += data[0:][0:]
    fh.close()

    imgavg=imgavg/NFILES

    hdulout = hdularr[NFILES // 2]
    head = hdulout[0].header

    hdulout[0].data = imgavg
    head.set('comment',
             'Mean combined from {} files. Header from central image ({}).'.format(NFILES, names[NFILES // 2]))

    set_bitpix(hdulout,BITPIX)
    hdulout.writeto(fout)
    for i in range(NFILES): hdularr[i].close()


def medianCombine(fout, filelist,BITPIX=-32):
    """Median combine files in file list takes (fout,filelist)"""

    fh = open(filelist, "r")
    NFILES = 0
    for line in fh: NFILES += 1
    fh.seek(0)
    print("Number of files is {}".format(NFILES))

    hdularr = []

    hdul = pyfits.open(fh.readline().rstrip())
    hdularr.append(hdul)

    imgdata0 = hdul[0].data
    NY, NX = imgdata0.shape
    print(NY, NX)

    imgcube = np.zeros([NFILES, NY, NX])
    print(imgcube[0].shape)

    i = 1
    names = []
    for line in fh:
        names.append(line.rstrip())
        hdul = pyfits.open(line.rstrip())
        hdularr.append(hdul)
        data = hdul[0].data
        print(data.shape)
        imgcube[i][0:][0:] = data[0:][0:]
        i += 1
    fh.close()

    # check lengths
    Ndata0 = imgcube.size
    for i in range(1, NFILES):
        Ndatai = imgcube.size
        print(Ndatai)
        if not Ndata0 == Ndatai:
            print("Data length mismatch between file 0 and file {}".format(i))
            sys.exit()

    hdulout = hdularr[NFILES // 2]
    head = hdulout[0].header
    meddata = np.zeros((NY,NX))

    meddata = np.median(imgcube,axis=0)
    print("size of median array is {}".format(meddata.shape))

    hdulout[0].data = meddata
    head.set('comment',
             'Median combined from {} files. Header from central image ({}).'.format(NFILES, names[NFILES // 2]))

    set_bitpix(hdulout,BITPIX)
    hdulout.writeto(fout)
    for i in range(NFILES): hdularr[i].close()


def skyMedianCombine(fout, filelist,BITPIX=-32):
    """Median combine sky flat files in file list takes (fout,filelist).
       Normalizes files before combine.
    """

    fh = open(filelist, "r")
    NFILES = 0
    for line in fh: NFILES += 1
    fh.seek(0)
    print("Number of files is {}".format(NFILES))

    hdularr = []

    hdul = pyfits.open(fh.readline().rstrip())
    hdularr.append(hdul)

    imgdata0 = hdul[0].data
    NY, NX = imgdata0.shape
    print(NY, NX)

    imgcube = np.zeros([NFILES, NY, NX])
    print(imgcube[0].shape)


    med=np.median(imgdata0)

    imgcube[0,0:,0:]=imgdata0[0:,0:]/med

    i = 1
    names = []
    for line in fh:
        names.append(line.rstrip())
        hdul = pyfits.open(line.rstrip())
        hdularr.append(hdul)
        data = hdul[0].data
        med = np.median(data)
        data = data/med
        print(data.shape)
        imgcube[i][0:][0:] = data[0:][0:]*1
        i += 1
    fh.close()

    # check lengths
    Ndata0 = imgcube.size
    for i in range(1, NFILES):
        Ndatai = imgcube.size
        print(Ndatai)
        if not Ndata0 == Ndatai:
            print("Data length mismatch between file 0 and file {}".format(i))
            sys.exit()

    hdulout = hdularr[NFILES // 2]
    head = hdulout[0].header

    meddata=(np.median(imgcube,axis=0))

    set_bitpix(hdulout,BITPIX)
    hdulout[0].data = meddata
    head.set('comment',
             'Median combined from {} files. Header from central image ({}).'.format(NFILES, names[NFILES // 2]))

    hdulout.writeto(fout)
    hdulout.close()
    for i in range(NFILES): hdularr[i].close()

def getMedian(filename):
    """Median of one file takes (filename)"""

    hdul = pyfits.open(filename)

    imgdata0 = hdul[0].data
    NY, NX = imgdata0.shape
    print(NY, NX)
    print("Median {} ".format(np.median(imgdata0)))
    hdul.close()

def simpleAperture(outfile,filelist,aperturelist):
    """You need to first create an aperture list that contains the pixel X Y 
       coordinates, the inner radius of the aperture, and the outer radius of the 
       aperture. The filelist should have the same order as the aperture list.  
       takes (outfile,filelist,aperturelist)
    """

    fh = open(filelist, "r")
    ah = open(aperturelist,"r")

    from astropy.wcs import WCS 


    iter = 0
    for line in fh:
        fpath = line.rstrip()
        hdul = pyfits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data

        NY, NX = imgdata.shape


        fhout = open(outfile+ "-" + repr(iter).zfill(6), "w")
        fhout.write("#{}\n".format(fpath))
        fcent = open(ah.readline().rstrip(),"r")
        for locals in fcent:
            X0,Y0,APERTURE,BAPERTURE=locals.split()
            X0=int(float(X0))
            Y0=int(float(Y0))            
            APERTURE=float(APERTURE)
            BAPERTURE=float(BAPERTURE)

            BLOOK = int(BAPERTURE) + 2
            bg = []
            cells = 0
            sum = 0.
            for j in range(Y0 - BLOOK, min(Y0 + BLOOK,NY-1)):
                for i in range(X0 - BLOOK, min(X0 + BLOOK,NX-1)):
                    r = np.sqrt((j - Y0) ** 2 + (i - X0) ** 2)
                    if r < BAPERTURE:

                        if r > APERTURE:
                            bg.append(imgdata[j][i])
                        else:
                            sum += imgdata[j][i]
                            cells += 1

            medbg = np.median(np.array(bg))
            adu = sum - cells * medbg

            sky=(WCS(head).pixel_to_world_values(X0,Y0))

            fhout.write(" {} {} {} {} {} {}\n ".format(X0,Y0,sky[0],sky[1],adu,-2.5*np.log10(adu)))
        fhout.close()

        fcent.close()
        hdul.close()
        iter += 1
    fh.close()
    ah.close()

def simpleApertureDeg(outfile,filelist,apertureFile,PSKIP=4,PANNULUS=4,ifstart=0,islist=True,Constant=0.,centering=False,astrometry=False,plot=False):
    """You need to first create an aperture list that contains the pixel X Y 
       coordinates, the inner radius of the aperture, and the outer radius of the 
       aperture. The filelist should have the same order as the aperture list.  
       takes (outfile,filelist,apertureFile)
    """

    ah = open(apertureFile,"r")

    from astropy.wcs import WCS 
    from astropy.coordinates import SkyCoord


    apertures=[]
    for line in ah:
       if line[0:6]=="circle":
          info,junk=line.split('")')
          val=info.lstrip("circle(").split(",")
          print(len(val))
          x0=float(val[0])
          y0=float(val[1])
          r0=float(val[2])
          r1=float(val[2])
          apertures.append({'x0':x0,'y0':y0,'r0':r0,'r1':r1})
    ah.close()
    for dicts in apertures: print(dicts)


    if islist: 
        fh = open(filelist, "r")
    else:
        fh = [filelist+"\n"]
    iter = 0
    for line in fh:
        fpath = line.rstrip()
        hdul = pyfits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data
        w=WCS(head)

        NY, NX = imgdata.shape

        try:
            arcsec_pix1 = np.sqrt(head['CD1_1']**2+head['CD1_2']**2)*3600
            arcsec_pix2 = np.sqrt(head['CD2_1']**2+head['CD2_2']**2)*3600
            arcsec_pix = (arcsec_pix1+arcsec_pix2)/2
        except:
            print("CD's not found in header. Using defined dx ={} and dy={} in arcsec".format(dx,dy))
            arcsec_pix = (pixdx+pixdy)/2
        print("Plate scale is set to {} ''/pixel".format(arcsec_pix))



        fhout = open(outfile+ "-" + repr(iter+ifstart).zfill(6), "w")
        fhout.write("#{}\n".format(fpath))
        for dicts in apertures:
            x0=dicts["x0"]
            y0=dicts["y0"]
            sky=SkyCoord(ra=x0,dec=y0,frame='icrs',unit='deg')
            print(sky)
            X0,Y0=w.world_to_pixel(sky)
            X0=int(np.round(X0))
            #X0=int(X0)-1
            Y0=int(np.round(Y0))
            #Y0=int(Y0)-1
            print(X0,Y0)


            APERTURE=dicts["r0"]/arcsec_pix
            BAPERTURE=dicts["r0"]/arcsec_pix+PANNULUS+PSKIP
            BLOOK = int(BAPERTURE) + 2

            if centering:
                sum=0.
                xsum=0.
                ysum=0.
                bg = []
                for j in range(Y0 - BLOOK, min(Y0 + BLOOK,NY-1)):
                    for i in range(X0 - BLOOK, min(X0 + BLOOK,NX-1)):
                        r = np.sqrt((j - Y0) ** 2 + (i - X0) ** 2)
                        if r < BAPERTURE:

                            if r > APERTURE+PSKIP:
                                bg.append(imgdata[j][i])

                medbg = np.median(np.array(bg))
#                xproj=[]
#                yproj=[]
#                for j in range(Y0 - BLOOK, min(Y0 + BLOOK,NY-1)):
#                    xsum=0
#                    for i in range(X0 - BLOOK, min(X0 + BLOOK,NX-1)):
#                        r = np.sqrt((j - Y0) ** 2 + (i - X0) ** 2)
#                        if r < APERTURE:
#                            xsum+=imgdata[j][i]-medbg
#                    yproj.append(xsum)
#                for i in range(X0 - BLOOK, min(X0 + BLOOK,NX-1)):
#                    ysum=0
#                    for j in range(Y0 - BLOOK, min(Y0 + BLOOK,NY-1)):
#                        r = np.sqrt((j - Y0) ** 2 + (i - X0) ** 2)
#                        if r < APERTURE:
#                            ysum+=imgdata[j][i]-medbg
#                    xproj.append(ysum)
#
#                xproj=np.array(xproj)
#                yproj=np.array(yproj)
#                print(xproj)
#                print(yproj)
 #               print(xproj.argmax())
 #               X0=X0-BLOOK+xproj.argmax()
 #               Y0=Y0-BLOOK+yproj.argmax()
                       
                for j in range(Y0 - BLOOK, min(Y0 + BLOOK,NY-1)):
                    for i in range(X0 - BLOOK, min(X0 + BLOOK,NX-1)):
                        r = np.sqrt((j - Y0) ** 2 + (i - X0) ** 2)
                
                        if r < APERTURE:
                            sum += (imgdata[j][i]-medbg)
                            xsum += (i)*(imgdata[j][i]-medbg)
                            ysum += (j)*(imgdata[j][i]-medbg)

                X0 = (xsum/sum)
                Y0 = (ysum/sum)

                print(X0,Y0)
                sky1=w.pixel_to_world_values(X0,Y0)
                X0=int(np.round(X0))
                Y0=int(np.round(Y0))
                print(X0,Y0)
                print(sky1[0],sky1[1])
                hrs=sky1[0]/15
                hr =int(hrs)
                minutes=(hrs-hr)*60
                minute =int(minutes)
                secs = (minutes-minute)*60
                sec = np.round(secs,2)
                deg = int(sky1[1])
                aminutes = abs(sky1[1]-deg)*60
                aminute = int(aminutes)
                asecs = (aminutes-aminute)*60
                asec = np.round(asecs,1)
                astrometryString="RA {}:{}:{}  DEC {}:{}:{}".format(hr,minute,sec,deg,aminute,asec)
                print(astrometryString)


            bg = []
            cells = 0
            sum = 0.
            for j in range(Y0 - BLOOK, min(Y0 + BLOOK,NY-1)):
                for i in range(X0 - BLOOK, min(X0 + BLOOK,NX-1)):
                    r = np.sqrt((j - Y0) ** 2 + (i - X0) ** 2)
                    if r < BAPERTURE:

                        if r > APERTURE+PSKIP:
                            bg.append(imgdata[j][i])
                        else:
                            sum += imgdata[j][i]
                            cells += 1

            medbg = np.median(np.array(bg))
            adu = sum - cells * medbg
            print("S/N: {}".format(adu/np.sqrt(sum)))
            if plot:
                  apcount=[]
                  rcount=[]
                  for j in range(Y0 - BLOOK, min(Y0 + BLOOK,NY-1)):
                      for i in range(X0 - BLOOK, min(X0 + BLOOK,NX-1)):
                          r = np.sqrt((j - Y0) ** 2 + (i - X0) ** 2)
                
                          if r < BAPERTURE:
                              apcount.append(imgdata[j][i]-medbg)
                              rcount.append(r)
                  import matplotlib.pylab as plt
                  plt.figure()
                  plt.xlabel("Radius pixels")
                  plt.ylabel("Pixel ADU")
                  plt.scatter(rcount,apcount)
                  plt.show()
   



            sky=(WCS(head).pixel_to_world_values(X0,Y0))

            fhout.write(" {} {} {} {} {} {}\n ".format(X0,Y0,sky[0],sky[1],adu,-2.5*np.log10(adu)+Constant))
            if astrometry: fhout.write(astrometryString+"\n")
        fhout.close()

        hdul.close()
        iter += 1
    if islist:fh.close()
    fhout.close()


def getMags(fout,photlist,rad=0.0014,fltr='g',tol=1e-4):
    """ Get magnitudes for objects in aperture file and compare with Gaia catalogue. 
        Requires apiEx and assumes images are astrometrically calibrated"""
    import apiEx
    phname=open(photlist,"r")

    for iname,line in enumerate(phname):
      fh=open(line.rstrip(),"r")

      print(fh)
      magConst=0
      iter=0
      out="# x y ra dec adu magInst G GBP GRP g magConst \n"
      const=[]

      wsum=0.
      adu_sum=0.
      for item in fh:
         if item[0]=="#":continue
         if len(item)<=2:continue

         try: x,y,ra,dec,adu,magInst=item.split()
         except: continue
         magInst=float(magInst)
         gaia_json=apiEx.gaiaDR3_cone_search(ra,dec,rad,pgsize=1000,page=1)
         gaia_table=apiEx.mast_json2table(gaia_json)
         if len(gaia_table)>1:
           id=np.argsort(gaia_table[:]["phot_g_mean_mag"])
           row=gaia_table[id[0]]
         else: row=gaia_table[0]
         G  =row["phot_g_mean_mag"]
         GBP=row["phot_bp_mean_mag"]
         GRP=row["phot_rp_mean_mag"]

         col = GBP-GRP

         if fltr=='g':
             #Gminusg=0.13518 - 0.46245*(col) -0.25171*col**2 + 0.021349*col**3
             Gminusg=0.2199  -0.6365 *(col) -0.1548 *col**2 + 0.0064*col**3
             magout = G-Gminusg
             const.append(magout-magInst)
             magConst+=magout-magInst
  
         elif fltr=='I':
             GminusI= 0.01753 + 0.76*(col) -0.0991*col**2 
             magout = G-GminusI
             const.append(magout-magInst)
             magConst+=magout-magInst


         elif fltr=='R':
             # fix at some point for correct power multiplication
             GminusR= -0.02275 + col*(0.3961 + col*(-0.1243 + col*( -0.01396 + col*0.003775)))
             magout = G-GminusR
             const.append(magout-magInst)
             magConst+=magout-magInst

         elif fltr=='V':
             GminusV= -0.02704 + 0.01424*(col) -0.2156*col**2 + 0.01426*col**3
             magout = G-GminusV
             const.append(magout-magInst)
             magConst+=magout-magInst

         elif fltr=='B':
             GminusV= -0.02704 + 0.01424*(col) -0.2156*col**2 + 0.01426*col**3
             diff=999
             BminusV0=GminusV*1
             while diff > tol:
                 BminusV  = -(GminusV + 0.04749 +0.2901*BminusV0**2 - 0.02008*BminusV0**3)/0.0124
                 diff = np.abs(BminusV-BminusV0)
                 BminusV0=BminusV*1

             magout = G-GminusV + BminusV

             const.append(magout-magInst)
             magConst+=magout-magInst

         else: 
             raise ValueError('Invalid fiter')

         adu=float(adu)
         adu_sum+=adu
         wsum+=(magout-magInst)*adu

         out+="{} {} {} {} {} {} {} {} {} {} {}\n".format(x,y,ra,dec,adu,magInst,G,GBP,GRP,magout,magout-magInst)
         iter+=1

      fh.close()
      const=np.array(const)
      fhout=open(fout+ "-" + repr(iname).zfill(6), "w")
      out+="#Weighted m-magInst {} average {} median {} and std {}\n".format(wsum/adu_sum,magConst/iter,np.median(const),np.std(const))
      print(out)   
      fhout.write(out)
      fhout.close()
    phname.close()
            

def regionLengths(regionfile):
    """ get region width in arcseconds. Assumes four corners and degree in for reading """
    fh = open(regionfile,"r")
    for line in fh:
       if line[0:7]=="polygon":
          val=line.rstrip(")\n").lstrip("polygon(").split(",")
          print(len(val))
          for i in range(2,len(val)+1,2):
            x0=float(val[i-2])
            y0=float(val[i-1])
            if i < len(val):
              x1=float(val[i])
              y1=float(val[i+1])
            else:
              x1=float(val[i-len(val)])
              y1=float(val[i+1-len(val)])
            d = np.sqrt( (x0-x1)**2*(np.cos(y0*np.pi/180)+np.cos(y1*np.pi/180))**2/4+(y0-y1)**2)
            print("Segment length {} arcsec".format(d*3600))
          
    fh.close()

def rollImage(fout,fin,RISE,RUN,BITPIX=-32):
     hdul = pyfits.open(fin)
     head = hdul[0].header
     imgdata = hdul[0].data
     NY, NX = imgdata.shape

     shifted = np.zeros((NY,NX))

     for iy in range(NY):
       xvals = imgdata[iy,:]
       SHIFT = int(iy*RUN/RISE)
       xvals_roll = np.roll(xvals,SHIFT)
       for ix in range(NX):shifted[iy,ix]=xvals_roll[ix]

     line = np.zeros(NX)
     for iy in range(NY):
       for ix in range(NX):
         line[ix]+=shifted[iy,ix]

     import matplotlib.pylab as plt
     plt.figure()
     plt.plot(line)
     plt.show()


     set_bitpix(hdul,BITPIX)
     hdul[0].data =shifted
     hdul.writeto(fout)


def gradientMap(foutresid,fin,BITPIX=-32):
     """ Shows gradients """
     hdul = pyfits.open(fin)
     head = hdul[0].header
     imgdata = hdul[0].data
     NY, NX = imgdata.shape

     resid = np.zeros((NY,NX))

     for iy in range(NY):
      ylow = max(0,iy-1)
      yhi  = min(NY-1,iy+1)
      for ix in range(NX):

        xlow = max(0,ix-1)
        xhi =  min(NX-1,ix+1)

        gradavg  = abs( imgdata[iy,ix]-imgdata[ylow,ix])
        gradavg += abs( imgdata[iy,ix]-imgdata[yhi,ix])
        gradavg += abs( imgdata[iy,ix]-imgdata[iy,xlow])
        gradavg += abs( imgdata[iy,ix]-imgdata[iy,xhi])
        
        resid[iy,ix]=gradavg*0.25

     set_bitpix(hdul,BITPIX)
     hdul[0].data = resid
     hdul.writeto(foutresid)


def gradientMask(fout,foutresid,fin,sig=3,frac=0.1,limit=2000,BITPIX=-32):
     """ Picks out gradients. Intended for eventual masking for artefact removal"""
     hdul = pyfits.open(fin)
     head = hdul[0].header
     imgdata = hdul[0].data
     NY, NX = imgdata.shape

     resid = np.zeros((NY,NX))+imgdata
     mask = np.zeros((NY,NX))

     for iy in range(NY):
      for ix in range(NX):

        if imgdata[iy,ix] < limit: continue
        maxval = imgdata[iy,ix]
        gradmax = maxval/sig**2*np.exp(-0.5/sig**2) # assumes gradient over 1 pixel. Just magnitude
      
        gradavg=0
        iavg=0
        for j in range(iy-1,iy+2):
          if j > NY-1 or j < 0: continue
          for i in range(ix-1,ix+2):
             if i > NX-1 or i < 0: continue
             gradavg = abs(maxval-imgdata[j,i])
             iavg+=1

        if iavg==0:continue
        gradavg/=iavg
        meds=[]
        if gradavg > 2*gradmax: # found artefact
          for j in range(iy-5,iy+6):
            if j > NY-1 or j < 0: continue
            for i in range(ix-5,ix+6):
               if i > NX-1 or i < 0: continue
               r = np.sqrt( (iy-j)**2+(ix-i)**2 )
               mask[j][i]=imgdata[j,i]

               if r > 3 and r < 5: meds.append(resid[j,i])

          meds=np.array(meds)
          med=np.median(meds)
          for j in range(iy-5,iy+6):
           if j > NY-1 or j < 0: continue
           for i in range(ix-5,ix+6):
              if i > NX-1 or i < 0: continue
              r = np.sqrt( (iy-j)**2+(ix-i)**2 )
              if r <= 3:
                 resid[j,i]=med 

 
     set_bitpix(hdul,BITPIX)
     hdul[0].data = mask
     hdul.writeto(fout)
  
     set_bitpix(hdul,BITPIX) 
     hdul[0].data = resid
     hdul.writeto(foutresid)

def gaussMask(fout,foutresid,fin,sig=3,frac=0.1,limit=50,bkgrnd=1000,BITPIX=-32):
     "In progress"
     from scipy.optimize import curve_fit
     import matplotlib.pylab as plt
     hdul = pyfits.open(fin)
     head = hdul[0].header
     imgdata = hdul[0].data
     NY, NX = imgdata.shape

     resid = np.zeros((NY,NX))+imgdata-bkgrnd
     stars = np.zeros((NY,NX))

     maxval=10*limit

     def gfunc(x,sig):
           return np.exp(-0.5*(x/sig)**2)

     while maxval > limit:

        iy,ix = np.unravel_index(np.argmax(resid),(NY,NX))
        maxval = np.max(resid)
        print(maxval,iy,ix)
        ylow = max(0,iy-5*sig)
        yhi  = min(NY-1,iy+5*sig)
        xlow = max(0,ix-5*sig)
        xhi  = min(NX-1,ix+5*sig)
 
        ydatay = resid[ylow:yhi+1,ix]/maxval
        ydatax = resid[iy,xlow:xhi+1]/maxval
        xdatay = np.linspace(ylow-iy,yhi-iy,(yhi-ylow)+1)
        xdatax = np.linspace(xlow-ix,xhi-ix,(xhi-xlow)+1)

        print(xdatax)
        print(ydatay)

        OK=1
        try:
          popty,pcovy = curve_fit(gfunc,xdatay,ydatay)
          poptx,pcovx = curve_fit(gfunc,xdatax,ydatax)

          print(popty,poptx)
          avgsig = 0.5*np.sqrt((popty[0]**2+poptx[0]**2))
          if popty[0]<sig-1 or poptx[0]<sig-1: OK=0
        except: 
          OK=0
          avgsig=sig

        #plt.figure()
        #plt.plot(xdatax,ydatax)
        #plt.plot(xdatay,ydatay)
        #plt.plot(xdatax,gfunc(xdatax,sig=avgsig),lw=2)
        #plt.show()

        if not OK:
           ylow = max(0,iy-3*sig)
           yhi  = min(NY-1,iy+3*sig)
           xlow = max(0,ix-3*sig)
           xhi  = min(NX-1,ix+3*sig)
           med = np.median(resid[ylow:yhi+1,xlow:xhi+1])

           resid[ylow:yhi+1,xlow:xhi+1] = med

        else: 

           ylow = max(0,iy-5*sig)
           yhi  = min(NY-1,iy+5*sig)
           xlow = max(0,ix-5*sig)
           xhi  = min(NX-1,ix+5*sig)

           for y in range(ylow,yhi+1):
             for x in range(xlow,xhi+1):

                val = maxval*np.exp( -0.5 * ( (y-iy)**2 + (x-ix)**2 )/avgsig**2 )
                stars[y,x] += val
                resid[y,x] -= val

 
     set_bitpix(hdul,BITPIX)
     hdul[0].data = stars
     hdul.writeto(fout)
  
     set_bitpix(hdul,BITPIX)
     hdul[0].data = resid+bkgrnd
     hdul.writeto(foutresid)

def gaussFilter(fout,fin,sig=3,bkgrnd=1000,BITPIX=-32):
     """ In progress. """
     hdul = pyfits.open(fin)
     head = hdul[0].header
     imgdata = hdul[0].data
     NY, NX = imgdata.shape

     mask = np.zeros((NY,NX))
     stars = np.zeros((NY,NX))

     for y in range(NY):
       for x in range(NX):
          mask[y,x]=np.exp(-0.5* ( (y-NY/2)**2 + (x-NX/2)**2 )/sig**2 )

     maskfft = np.fft.fft2(mask)
     imgfft = np.fft.fft2(imgdata-bkgrnd)

     V = imgfft*maskfft

     cc = np.fft.ifftshift(np.fft.ifft2(V)).real

     set_bitpix(hdul,BITPIX)
     hdul[0].data = cc
     hdul.writeto(fout)
  

def gaussRemove(fout,foutresid,fin,sig=3,frac=0.1,limit=1020,BITPIX=-32):
     """ Experimental """
     hdul = pyfits.open(fin)
     head = hdul[0].header
     imgdata = hdul[0].data
     NY, NX = imgdata.shape

     resid = np.zeros((NY,NX))+imgdata
     stars = np.zeros((NY,NX))

     maxval=10*limit

     while maxval > limit:

        iy,ix = np.unravel_index(np.argmax(resid),(NY,NX))
        maxval = np.max(resid)
        print(maxval,iy,ix)
        for j in range(iy-5*sig,iy+5*sig+1):
          if j > NY-1 or j < 0: continue
          for i in range(ix-5*sig,ix+5*sig+1):
             if i > NX-1 or i < 0: continue
             loc=maxval*frac*np.exp( -0.5* ( (j-iy)**2 + (i-ix)**2 )/sig**2 )
             stars[j][i]+=loc
             resid[j][i]-=loc
                
 
     set_bitpix(hdul,BITPIX)
     hdul[0].data = stars
     hdul.writeto(fout)
  
     set_bitpix(hdul,BITPIX)
     hdul[0].data = resid
     hdul.writeto(foutresid)
  
def autocorrelate(fout,infile,BITPIX=-32):
        """ Attempt to clip data in frequency space. Work in progress"""

        iter=0
        import matplotlib.pylab as plt
        import scipy.signal
        hdul = pyfits.open(infile)

        head = hdul[0].header
        imgdata = hdul[0].data
         
        F = np.fft.fftshift(np.fft.fft2(imgdata))

        S = F*F
        #S = F*F.conjugate()

        V=np.fft.ifft2(np.fft.ifftshift(S)).real
        #V = np.correlate(imgdata,imgdata)
        #V = scipy.signal.correlate2d(imgdata,imgdata)

        hdul[0].data=V
        
        if BITPIX==16: hdul[0].scale('int16')
        elif BITPIX==32: hdul[0].scale('int32')
        elif BITPIX==-32: hdul[0].scale('float32')
        elif BITPIX==-64: hdul[0].scale('float64')
        else:
            print("Invalid selection for BITPIX {}".format(BITPIX))
        
        set_bitpix(hdul,BITPIX)
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()


def fft(fout,infile,BITPIX=-32):
        """ Attempt to clip data in frequency space. Work in progress"""

        iter=0
        import matplotlib.pylab as plt
        hdul = pyfits.open(infile)

        head = hdul[0].header
        imgdata = hdul[0].data
         
        F = np.fft.fftshift(np.fft.fft2(imgdata))

        Freal = F.real
        Fimg = F.imag

        NY,NX = F.shape

        F[NY//2,NX//2]=0 # clip centre value with the idea to remove constants in image space.

        V=np.fft.ifft2(np.fft.ifftshift(F)).real

        hdul[0].data=V
        
        set_bitpix(hdul,BITPIX) 
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()

 
def Quantiles(infile,q=0.5):
        """ get count quantiles"""

        iter=0
        import matplotlib.pylab as plt
        hdul = pyfits.open(infile)

        head = hdul[0].header
        imgdata = np.ravel(hdul[0].data)
        hdul.close()
        return np.quantile(imgdata,q)

def Sharpness(infile):
        """ get count quantiles for determining image sharpness"""
        from scipy.ndimage import median_filter
        iter=0
        import matplotlib.pylab as plt
        hdul = pyfits.open(infile)

        head = hdul[0].header
        imgdata = hdul[0].data
        hdul.close()
        filtered = median_filter(imgdata,size=(3,3))
        dx,dy=np.gradient(filtered)
        dnorm = np.average(np.sqrt(dx*dx+dy*dy))
        return dnorm

def traceStreak(imagefile,regionfile,outtable="trace.tbl",fout=None,width=8,overwidth=16,sampleWidth=40,BITPIX=-32,rate=1,pixdx=1,pixdy=1,exposure=1,calMag=0,PLOT=False):
    """ region file in pixels """

    hdul = pyfits.open(imagefile)
    head = hdul[0].header
    imgdata = hdul[0].data
    NY, NX = imgdata.shape

    try:
        arcsec_pix1 = np.sqrt(head['CD1_1']**2+head['CD1_2']**2)*3600
        arcsec_pix2 = np.sqrt(head['CD2_1']**2+head['CD2_2']**2)*3600
        arcsec_pix = (arcsec_pix1+arcsec_pix2)/2
    except:
        print("CD's not found in header. Using defined dx ={} and dy={} in arcsec".format(dx,dy))
        arcsec_pix = (pixdx+pixdy)/2
    print("Plate scale is set to {} ''/pixel".format(arcsec_pix))

    try:
        exptime=head['EXPTIME']
    except:
        print("EXPTIME not in header. Using defined exposure = {} s".format(exposure))
        exptime=exposure
    print("Exposure time used for calcs is {} s".format(exptime))


    fh = open(regionfile,"r")
    for line in fh:
       if line[0:4]=="line":
          info,junk=line.split(")")
          val=info.lstrip("line(").split(",")
          print(len(val))
          x0=float(val[0])-1
          y0=float(val[1])-1
          x1=float(val[2])-1
          y1=float(val[3])-1

          Dx=x1-x0
          Dy=y1-y0 # get slope terms

          alpha = np.pi/2+np.arctan2( Dy, Dx )
          dx=np.cos(alpha)
          dy=np.sin(alpha)

          #xu0 = x0 + dx; xu1 = x1 + dx
          #yu0 = y0 + dy; yu1 = y1 + dy

          #xl0 = x0 - dx; xl1 = x1 - dx
          #yl0 = y0 - dy; yl1 = y1 - dy

    mask = np.zeros( (NY,NX) )
    mask_bg = np.zeros( (NY,NX) )

   
    lengthMask = int(np.sqrt(Dx**2 + Dy**2))
    print(x0,y0,x1,y1,lengthMask)
    alpha90 = alpha-np.pi/2

    idx_source =[]
    idx_bg =[]

    twoWidth=int(2*width)
    for l in range(-twoWidth,twoWidth+1):
    #for l in range(0,1):
        for w in range(lengthMask):
             xs0 = int(x0 + (l-0.5)*dx*0.5 + w*np.cos(alpha90))
             ys0 = int(y0 + (l-0.5)*dy*0.5 + w*np.sin(alpha90)) 

             try:
                 mask[ys0,xs0] = 1
             except:continue

    for j in range(NY):
        for i in range(NX):
            if mask[j,i]>0: idx_source.append([i,j])

    for l in range(int(-2*overwidth),-twoWidth):
    #for l in range(0,1):
        for w in range(lengthMask):
             xs0 = int(x0 + (l-0.5)*dx*0.5 + w*np.cos(alpha90))
             ys0 = int(y0 + (l-0.5)*dy*0.5 + w*np.sin(alpha90)) 

             try:
                 mask_bg[ys0,xs0] = 1
             except:continue

    for l in range(twoWidth+1,int(2*overwidth+1)):
    #for l in range(0,1):
        for w in range(lengthMask):
             xs0 = int(x0 + (l-0.5)*dx*0.5 + w*np.cos(alpha90))
             ys0 = int(y0 + (l-0.5)*dy*0.5 + w*np.sin(alpha90)) 

             try:
                 mask_bg[ys0,xs0] = 1
             except:continue

    for j in range(NY):
        for i in range(NX):
            if mask_bg[j,i]>0: idx_bg.append([i,j])


    #idx_source = np.array(idx_source)
    #idx_bg     = np.array(idx_bg)

    idx_source_rot = []
    idx_bg_rot     = []

    print("length of idx_source = {}".format(len(idx_source)))

    def rotate(x,y,a):
        return x*np.cos(a)-y*np.sin(a), x*np.sin(a) + y*np.cos(a)

    for i in range(len(idx_source)):
        x,y=idx_source[i]
        idx_source_rot.append(rotate(x-x0,y-y0,-alpha90))
    for i in range(len(idx_bg)):
        x,y=idx_bg[i]
        idx_bg_rot.append(rotate(x-x0,y-y0,-alpha90))

    idx_source_rot = np.array(idx_source_rot)
    idx_bg_rot = np.array(idx_bg_rot)
    xs=[];ys=[]
    xb=[];yb=[]
    for xy in idx_source_rot:
        xs.append(xy[0])
        ys.append(xy[1])
    for xy in idx_bg_rot:
        xb.append(xy[0])
        yb.append(xy[1])

    xs=np.array(xs)
    ys=np.array(ys)
    xb=np.array(xb)
    yb=np.array(yb)

    i_xs_sorted = np.argsort(xs)
    i_xb_sorted = np.argsort(xb)

    ix_last=0
    isample=0
    line_flux=[]
    line_r=[]
    line_c=[]
    flxs_median=[]
    flxs_mean=[]
    while True:
        xblock_start=isample*sampleWidth
        xblock_end=(isample+1)*sampleWidth
        if xblock_end > xs[i_xs_sorted[-1]]:break

        s=ix_last*1
        xval = xs[i_xs_sorted[s]]
        flux = 0
        count=0
        flxs=[]
        while xval < xblock_end:
            u,v = idx_source[i_xs_sorted[s]]
            flxs.append(imgdata[v,u])
            flux+=imgdata[v,u]
            count+=1
            s+=1
            xval = xs[i_xs_sorted[s]]
        ix_last=s*1
        flxs=np.sort(np.array(flxs))
        flxs_median.append(np.median(flxs[-5:]))
        flxs_mean.append(np.mean(flxs[-5:]))
        line_flux.append(flux)
        line_r.append(xblock_start+sampleWidth*0.5)
        line_c.append(count)
        isample+=1

    ix_last=0
    isample=0
    bg_median=[]
    while True:
        xblock_start=isample*sampleWidth
        xblock_end=(isample+1)*sampleWidth
        if xblock_end > xs[i_xs_sorted[-1]]:break

        s=ix_last*1
        xval = xb[i_xb_sorted[s]]
        bgs=[]
        count=0
        while xval < xblock_end:
            u,v = idx_bg[i_xb_sorted[s]]
            bgs.append(imgdata[v,u])
            count+=1
            s+=1
            xval = xb[i_xb_sorted[s]]
        bgs=np.array(bgs)
        bg_median.append(np.median(bgs))
        isample+=1

    line_flux=np.array(line_flux)
    line_r=np.array(line_r)
    line_c=np.array(line_c)
    bg_median=np.array(bg_median)
    flxs_median=np.array(flxs_median)

    line_nflux=(line_flux-bg_median*line_c)/(sampleWidth)
    line_reduced_flux = line_nflux*exptime*rate/arcsec_pix
    line_mag=-2.5*np.log10(line_reduced_flux) + calMag
    line_surfbavg = -2.5*np.log10(line_nflux*sampleWidth/(line_c*arcsec_pix**2)) + calMag
    line_surfb_md = -2.5*np.log10((flxs_median-bg_median)/arcsec_pix**2) + calMag
    line_surfb_mn = -2.5*np.log10((flxs_mean-bg_median)/arcsec_pix**2) + calMag

    if PLOT==True:
        import matplotlib.pylab as plt
        plt.figure()
        plt.scatter(xs,ys,s=0.2)
        plt.scatter(xb,yb,color='red',s=0.2)

        plt.figure()
        plt.title("Length normalized line flux")
        plt.xlabel('Pixels Along Streak')
        plt.ylabel('Pixel count divided by box width')
        plt.plot(line_r,line_nflux)

        plt.figure()
        plt.xlabel('Pixels Along Streak')
        plt.ylabel('Pixels in box')
        plt.title("Pixel counts per sample")
        plt.plot(line_r,line_c)

        plt.figure()
        plt.title("Surface Brightness per Sq. Arcsec")
        plt.xlabel('Pixels Along Streak')
        plt.ylabel('Mag per sq. arcsec')
        plt.plot(line_r,line_surfb_md)
        plt.plot(line_r,line_surfb_mn)
        #plt.plot(line_r,line_surfbavg)

        plt.figure()
        plt.title("Background Counts")
        plt.plot(line_r,bg_median)

        plt.figure()
        plt.title("Running magnitude")
        plt.xlabel('Pixels Along Streak')
        plt.ylabel('Mag')
        plt.plot(line_r,line_mag)

        plt.show()

    fo=open(outtable,"w")
    fo.write("#pixel, surfb_md, surfb_mn, mag_on_line\n")
    for i in range(len(line_r)):
        fo.write("{} {} {} {}\n".format(line_r[i],line_surfb_md[i],line_surfb_mn[i],line_mag[i]))
    fo.close()

    if fout is not None:
       mask_tot = mask + mask_bg
       mask_tot[mask_tot>1]=1
       V = imgdata*mask_tot
    

       hdul[0].data=V
        
       set_bitpix(hdul,BITPIX) 
       hdul.writeto(fout)
       hdul.close()

    fh.close()




        

         
