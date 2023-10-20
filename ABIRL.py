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

def addCoords(fout,filelist,pxu=1e-3):
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

        head['BITPIX']=-32
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

def medianFilter(fout,filelist,BLOCKX=1,BLOCKY=1):
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
         

        head['BITPIX']=-32
        print(head['BITPIX'])
        head.set('comment', 'Median filtered (despiked) from {}'.format(fpath))
        hdul[0].data = corrected
        

        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1


def chunkQuartileFlatten(fout,filelist,BLOCKX=32,BLOCKY=32,SKIP=1):
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


        head['BITPIX']=-32
        print(head['BITPIX'])
        head.set('comment', 'Quartile flattened from {}'.format(fpath))
        hdul[0].data = corrected
        

        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

 
def medianFlatten(fout,filelist,BLOCKX=32,BLOCKY=32):
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


        head['BITPIX']=-32
        print(head['BITPIX'])
        head.set('comment', 'Median flattened from {}'.format(fpath))
        hdul[0].data = corrected
        

        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

   


def chunkMedianSubtract(fout,filelist,BLOCKY=16,BLOCKX=8):
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

        head['BITPIX']=-32
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1



def addBase(fout,filelist,base=100):
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

        head['BITPIX']=-32
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1


def subtractImg(fout, fitsToSubtract, filelist,catch=1):
    """Subtracts one image from all images in a list. takes(fout,fitsToSubtract,filelist)"""

    hdulref = pyfits.open(fitsToSubtract)
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

        corrected = np.zeros((NY,NX)) # defining because I was getting some weird errors on occasion. 
        corrected = imgdata - refdata

        head.set('comment', 'Subtracted {} from file {}'.format(fitsToSubtract, fpath))
        hdul[0].data = corrected

        head['BITPIX']=-32
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

def divideImg(fout, fitsDenominator, filelist):
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

        head['BITPIX']=-32
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()
        iter += 1

def meanCombine(fout, filelist):
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
        imgavg[0:-1][0:-1] += data[0:-1][0:-1]
    fh.close()

    imgavg=imgavg/NFILES

    hdulout = hdularr[NFILES // 2]
    head = hdulout[0].header

    hdulout[0].data = imgavg
    head.set('comment',
             'Mean combined from {} files. Header from central image ({}).'.format(NFILES, names[NFILES // 2]))

    hdulout.writeto(fout)
    for i in range(NFILES): hdularr[i].close()



def medianCombine(fout, filelist):
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
        imgcube[i][0:-1][0:-1] = data[0:-1][0:-1]
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

    head['BITPIX']=-32
    hdulout.writeto(fout)
    for i in range(NFILES): hdularr[i].close()


def skyMedianCombine(fout, filelist):
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


    mean=np.mean(imgdata0)

    imgcube[0,0:-1,0:-1]=imgdata0[0:-1,0:-1]/mean

    i = 1
    names = []
    for line in fh:
        names.append(line.rstrip())
        hdul = pyfits.open(line.rstrip())
        hdularr.append(hdul)
        data = hdul[0].data
        mean = np.mean(data)
        data = data/mean
        print(data.shape)
        imgcube[i][0:-1][0:-1] = data[0:-1][0:-1]
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

    hdulout[0].data = meddata
    head.set('comment',
             'Median combined from {} files. Header from central image ({}).'.format(NFILES, names[NFILES // 2]))

    head['BITPIX']=0
    hdulout.writeto(fout)
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
    fhout.close()

def getMags(fout,photlist,rad=0.0014):
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

         #Gminusg=0.13518 - 0.46245*(col) -0.25171*col**2 + 0.021349*col**3
         Gminusg=0.2199  -0.6365 *(col) -0.1548 *col**2 + 0.0064*col**3
     
         g = G-Gminusg
         const.append(g-magInst)
         magConst+=g-magInst

         out+="{} {} {} {} {} {} {} {} {} {} {}\n".format(x,y,ra,dec,adu,magInst,G,GBP,GRP,g,g-magInst)
         
         iter+=1

      fh.close()
      const=np.array(const)
      fhout=open(fout+ "-" + repr(iname).zfill(6), "w")
      out+="#Average g-magInst {} median {} and std {}\n".format(magConst/iter,np.median(const),np.std(const))
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

def rollImage(fout,fin,RISE,RUN):
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

     head['BITPIX']=-32
     hdul[0].data =shifted
     hdul.writeto(fout)


def gradientMap(foutresid,fin):
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

     head['BITPIX']=-32
     hdul[0].data = resid
     hdul.writeto(foutresid)


def gradientMask(fout,foutresid,fin,sig=3,frac=0.1,limit=2000):
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

 
     head['BITPIX']=-32
     hdul[0].data = mask
     hdul.writeto(fout)
  
     hdul[0].data = resid
     hdul.writeto(foutresid)

def gaussMask(fout,foutresid,fin,sig=3,frac=0.1,limit=50,bkgrnd=1000):
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

 
     head['BITPIX']=-32
     hdul[0].data = stars
     hdul.writeto(fout)
  
     hdul[0].data = resid+bkgrnd
     hdul.writeto(foutresid)

def gaussFilter(fout,fin,sig=3,bkgrnd=1000):
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

     head['BITPIX']=-32
     hdul[0].data = cc
     hdul.writeto(fout)
  

def gaussRemove(fout,foutresid,fin,sig=3,frac=0.1,limit=1020):
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
                
 
     head['BITPIX']=-32
     hdul[0].data = stars
     hdul.writeto(fout)
  
     hdul[0].data = resid
     hdul.writeto(foutresid)
 
def fft(fout,infile):
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
        
        head['BITPIX']=-32
        
        hdul.writeto(fout + "-" + repr(iter).zfill(6) + ".fits")
        hdul.close()


