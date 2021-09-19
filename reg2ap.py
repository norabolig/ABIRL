# Aaron Boley, 15 Sep 2021. 
# Takes DS9 reg files in pixel coordinates for spheres and puts them into a format for simpleAperture in reductLib, assuming an outer aperture of +4 cells (can be manually altered, of course). 
import numpy as np

filelist=open("reg2ap.list","r")
foutPre="aperture"

iter=0
for line in filelist:
    fh = open(line.rstrip(),"r")
    fout = open(foutPre+"-"+repr(iter).zfill(6)+".ap","w")
    for row in fh:
        if not row[0:6]=="circle":continue
        cln1 = row.replace("circle(","")
        cln2 = cln1.rstrip().rstrip(")")
        x,y,ap1 = cln2.split(",")
        fout.write("{} {} {} {}\n".format(x,y,ap1,float(ap1)+4))
    fh.close()
    iter+=1

filelist.close()
