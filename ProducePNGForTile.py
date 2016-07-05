#!/usr/bin/env python
"""
 The simple program to create PNG images from stamps
 stored in multy extension fits file created by MakeCoaddCutouts.py program.
 It is also an example of how to work with such file.
 Input parameters are: -i <input fits file with tile images>
                       -p <PNG parameters file>
 An exzample of the PNG parameters file is as follow:
sCutMin=0    # number of background sigmas subtracted from image
sCutMax=500  # not used
satLevel=0.8 # saturation level of resulting image
enhLevel=0.2 # enhansement of the blue color factor
Qval=3.0     # Q factor of the scaling function
alpha=0.03   # alpha value of the scaling function
Rcor=0.8     # correction factors on separate bands
Gcor=1.4
Icor=1.
MakeTag=1    # if 1 PNG image will contain a tag, if 0 - no tags

 By N. Kuropatkin  06/21 2016
"""
import os
import sys
import fitsio
import exceptions
import string
import getopt
from MakeHumViImage import MakeHumViImage

class ProducePNGforTile(exceptions.Exception):

    def __init__(self, parfile):
        self.parfile = parfile
        self.MakeTag = 0
        self.inpar = {}
        self.headers = {}
    " Reads input parameter file like PNG.conf "    
    def getInput(self):
        f = open(self.parfile,"r")
        while 1:
            line = f.readline()
            if not line:
                break
            line = line[0:len(line)-1]
            tokens = string.split(line, "=", 1)
            key = tokens[0]
            value = tokens[1]
            self.inpar[key] = value
        self.sCutMin = self.inpar["sCutMin"]
        self.sCutMax = self.inpar["sCutMax"]
        self.satLevel = self.inpar["satLevel"]
        self.enhLevel = self.inpar["enhLevel"]
        self.Qval = self.inpar["Qval"]
        self.alpha = self.inpar["alpha"]
        self.Rcor = self.inpar["Rcor"]
        self.Gcor = self.inpar["Gcor"]
        self.Icor = self.inpar["Icor"]
        self.MakeTag = string.atoi(self.inpar["MakeTag"])
 
        f.close()   
        return(0)
    " Read fits file, extract header parameters and images "    
    def produce(self,infile):
        fitsfile = os.path.normpath(infile)
        fits = fitsio.FITS(fitsfile,'rw')

        prihdr = fits[0].read_header()
        self.ncol = prihdr["NAXIS1"]
        self.nrow = prihdr["NAXIS2"]
        self.objects = []
        self.images = {}
        self.objext = {}
        
        " loop on all hdus and create a list of objects and its extensions"
        for hdu in fits:
            header = hdu.read_header()
            objectN = header["Object"]
            if objectN not in self.objects:
                self.objects.append(objectN)
                self.objext[objectN] = {'g':0,'r':0,'i':0}
        fits.close()    

#        print self.objects
        " Loop on objects and select g,r,i images for each "
        fits = fitsio.FITS(fitsfile,'rw')
        ext = 0
        for hdu in fits:
            header = hdu.read_header()
            curObj = header["Object"]
            curBand = header["band"]
            typeE = header["TYPE"] # possible 'IMAGE' 'WEIGHT' 'PSF'  

            if curObj in self.objects:
                if typeE.find("IMAGE") >=0:
                    objind = self.objext.get(curObj)
                    if curBand.find('g') >= 0: objind['g'] = ext
                    if curBand.find('r') >= 0: objind['r'] = ext
                    if curBand.find('i') >= 0: objind['i'] = ext                   
                    self.objext[curObj] = objind
            ext+=1        
        fits.close()
        " loop on objects and create PNG images "
        fits = fitsio.FITS(fitsfile,'rw')
        for objectN in self.objects:
            objext = self.objext[objectN]
            outDir = "./gallery/"
            if not os.path.exists(outDir):
                os.makedirs(outDir)
            outFileName = "./gallery/"+objectN+'.png'
            imdict = {}
            extg = objext['g']
            header = fits[extg].read_header()
            imdict['g'] = fits[extg].read()
            extr = objext['r']
            imdict['r'] = fits[extr].read()
            exti = objext['i']
            imdict['i'] = fits[exti].read()
            makePNG = MakeHumViImage("")
            makePNG.setHead(header)
            makePNG.setParameters(self.Qval, self.alpha,imdict)
            makePNG.setSigmaCut(self.sCutMin)
            makePNG.setSaturationCut(self.satLevel)
            makePNG.setCorrFactors(self.Rcor, self.Gcor, self.Icor)
            makePNG.setEnhLevel(self.enhLevel)
            if self.MakeTag == 1:
                makePNG.setTag(objectN)
            makePNG.produce(outFileName)
        fits.close()
            
if __name__ == "__main__":
    print sys.argv
    outFile=""
    infile=""
    nbpar = len(sys.argv)
    if nbpar < 2:
        "Usage: ProducePNGForTile.py <required inputs>"
        print "  Required inputs:"
        print "  -i <input tile file> - like DES0118-5205_cutouts.fits"
        print "  -p <PNG parameters file> - like PNG.conf"
        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:p:",
                                  ["infile=","parfile="])
    except getopt.GetoptError:
        print "ProducePNGforTile.py -i <infile>" 
        sys.exit(2)
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "ProducePNGforTile.py -i <infile> -p <parameters file> " 
            sys.exit()
        elif opt in ("-i","--infile"):
            print "got -i arg= %s"%arg
            infile=arg
        elif opt in ("-p","--parfile"):
            print "got -p arg= %s"%arg
            parfile=arg
    prodPNG = ProducePNGforTile(parfile)
    prodPNG.getInput()
    prodPNG.produce(infile)
    sys.exit(0)