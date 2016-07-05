#!/usr/bin/env python
'''
Created on June 21, 2016
Modified on Jan 9 2014 to add wcs coordinates in the image metadata as suggested in
 http://blog.modp.com/2007/08/python-pil-and-png-metadata-take-2.html   
  
The program to create RGB image  from 3 fits files
representing  R,G,I or Z filters.
The input stamps should be from RGI or RGZ filters and will be
mapped to RGB colors of the output image.
The program provide HumVi type (Lupton's) scaling of the image 
 The image will be cut at given saturation level to enhance faint objects . 
 
Input parameters are provided in a single input file similar to following:
fileNameR=./data/out_r.fits
fileNameG=./data/out_g.fits
fileNameI=./data/out_i.fits
outFileName=testStamp
sigmaMin=1     # level of noise in the image
Rcor=1.0       # these correction factors can be considered as quantum efficiency
Gcor=1.0       # of the CCD to given band
Icor=1.0
enhLev=0.2     # additional enhancement of red and blue for contrast in search of strong lensing
               # one can play with value 0 - 0.2 and sign  of the correction. 0.2 will increase 
               # blue and decrease red.
satCut=0.8     # changes maximal value before saturation. The lower its value the brighter faint
               # objects will be and the more saturation.
Q=1.0          # Q factor in Lupton's scaling. Increase it to see fainter objects and more noise.
               # Recommended values 1.0 - 3.0
Alpha=0.06     # recommended values from 0.03  to 0.06. Also changes brightness and contrast.

Usage: ./MakeHumViImage.py -i inputParFile.txt


 
@author: kuropat
'''
import os
import sys
import pyfits
import exceptions
import string
import getopt
import math
import numpy as np
from PIL import Image
from PIL import PngImagePlugin
import ImageFont, ImageDraw




class MakeHumViImage(exceptions.Exception):
    
    def __init__(self, infile):
        self.infile = infile
        self.outfile = ""
        self.zp = 30.0
        self.enhLevel=0.0;
        self.Rcor=1.0
        self.Gcor=1.0
        self.Icor=1.0
        self.Q = 1.0
        self.Alpha = 0.06
        self.satCut = 0.8
        self.filters = ['R','G','I']
        self.crpix1 = 0
        self.crpix2 = 0
        self.crval1 = 0
        self.crval2 = 0
        self.cd1_1 = 0
        self.cd1_2 = 0
        self.cd2_1 = 0
        self.cd2_2 = 0
        self.scaleType = 0
        self.sigmaMin = 2
        self.sigmaMax = 6
        self.nrow = 0
        self.ncol = 0
        self.S = {}
        self.zpFf = {}
        self.zpF = {}
        self.inpar = {}
        self.stamps = {}
        self.impars = {}
        self.maxP = {}
        self.delta = {}
#        self.getInput()
        self.RAS = ''
        self.DECS = ''
        self.Tag = ''
        self.tagFlg = 0


    def printInputList(self):
        print "input file= "+self.infile
        
    ''' extract parameters from the file header '''
    def read(self,fitsfile):
        fitsfile = os.path.normpath(fitsfile)
        hdulist = pyfits.open(fitsfile)
        prihdr = hdulist[0].header
        self.ncol = prihdr["NAXIS1"]
        self.nrow = prihdr["NAXIS2"]
        expTime =  prihdr["EXPTIME"]
        self.crpix1 = prihdr["CRPIX1"]
        self.crpix2 = prihdr["CRPIX2"]
        self.crval1 = prihdr["CRVAL1"]
        self.crval2 = prihdr["CRVAL2"]
        self.cd1_1 = prihdr["CD1_1"]
        self.cd1_2 = prihdr["CD1_2"]
        self.cd2_1 = prihdr["CD2_1"]
        self.cd2_2 = prihdr["CD2_2"]
        try:
            self.zp = prihdr["SEXMGZPT"]
        except:
            self.zp = 30.

        image = hdulist[0].data
        hdulist.close()
        return image
    def setTag(self,tag):
        self.tagFlg = 1
        self.Tag = tag
        
    def makeTag(self,image):
        draw = ImageDraw.Draw(image)
        font = ImageFont.load_default()
        draw.text((1, 1),self.Tag,(255,255,255),font=font)
    
    def setHead(self,prihdr):
        self.ncol = prihdr["NAXIS1"]
        self.nrow = prihdr["NAXIS2"]
        expTime =  prihdr["EXPTIME"]
        self.crpix1 = prihdr["CRPIX1"]
        self.crpix2 = prihdr["CRPIX2"]
        self.crval1 = prihdr["CRVAL1"]
        self.crval2 = prihdr["CRVAL2"]
        self.cd1_1 = prihdr["CD1_1"]
        self.cd1_2 = prihdr["CD1_2"]
        self.cd2_1 = prihdr["CD2_1"]
        self.cd2_2 = prihdr["CD2_2"]
        try:
            self.zp = prihdr["SEXMGZPT"]
        except:
            self.zp = 30.
    '''
        processing input file
    '''   
    def getInput(self):
        f = open(self.infile,"r")
        while 1:
            line = f.readline()
            if not line:
                break
            line = line[0:len(line)-1]
            tokens = string.split(line, "=", 1)
            key = tokens[0]
            value = tokens[1]
            self.inpar[key] = value
        self.outfile = "./"+self.inpar["outFileName"]+".png"
        self.sigmaMin = string.atoi(self.inpar["sCutMin"])
#        self.sigmaMax = string.atoi(self.inpar["sCutMax"])
        self.enhLevel = string.atof(self.inpar["enhLevel"])
        self.Rcor = string.atof(self.inpar["Rcor"])
        self.Gcor = string.atof(self.inpar["Gcor"])
        self.Icor = string.atof(self.inpar["Icor"])
        self.satCut = string.atof(self.inpar["satCut"])
        self.Q = string.atof(self.inpar["Qval"])
        self.Alpha = string.atof(self.inpar["alpha"])
        filename=self.inpar["fileNameR"]
        self.stamps['R'] = self.read(filename)
        self.zpF["R"] = self.zp
        filename=self.inpar["fileNameG"]
        self.stamps['G'] = self.read(filename)
        self.zpF["G"] = self.zp
        filename=self.inpar["fileNameI"]
        self.stamps['I'] = self.read(filename)
        self.zpF["I"] = self.zp
        f.close()   
        return(0)
    
    def setSize(self,nrow,ncol):
        self.nrow = nrow
        self.ncol = ncol
        
    def setParameters(self,Q,alpha,imdict):
        self.imdict = imdict
        self.Q = string.atof(Q)
        self.Alpha = string.atof(alpha)
        self.stamps['R'] = self.imdict["r"]
        self.zpF["R"] = 30.
        self.stamps['G'] = self.imdict["g"]
        self.zpF["G"] = 30.
        self.stamps['I'] = self.imdict["i"]
        self.zpF["I"] = 30.
        
    def setSigmaCut(self,sigmaCut):
        self.sigmaMin = string.atof(sigmaCut)
    def setSaturationCut(self,sCut):
        self.sCut = string.atof(sCut)
    def setCorrFactors(self,Rcorr,Gcorr,Icorr):
        self.Rcor = string.atof(Rcorr)
        self.Gcor = string.atof(Gcorr)
        self.Icor = string.atof(Icorr)
    def setEnhLevel(self,enhLevel):
        self.enhLevel = string.atof(enhLevel)
               
    " calculate parameters of the image median, max value and sigma "
    def getStampPars(self,image):
        maxP = np.max(image)
        minP = np.min(image)
        subim = np.copy(image)
        subim = np.where(subim>=-19.9,subim,-19.9)
        subim = np.where(subim <= 20.,subim,20.)
        median = np.median(subim)
        sigma = np.std(subim)
#        print "im median=%f imSigma=%f maxP=%f"%(median,sigma,maxP)
        return [median,minP, maxP, sigma]


    '''
     Rescale image to make sigma = 1 and correct on the CCD efficiency
    '''  
    def filterIm(self, key, sigmaMin, sigmaMax, Fcor):
        median = self.impars[key][0]
        maxC = self.impars[key][2]
#        sigma = self.impars[key][3]
        image = self.stamps[key]
        a=(image + self.delta[key] -median)
        a = np.where(a>=0.,a,0.)
# 
        self.S[key] = Fcor 
        a = a*Fcor       
        self.maxP[key] = (maxC + self.delta[key] - median )*self.S[key]
#        print "Key= %s max= %f"%(key,np.max(a))
        return a 


    '''  Apply scaling factor and scale maximum to 255 '''    
    def rescale(self,key,X,scale):
        image = self.stamps[key]
        a = image*X
        a = a*scale
#        print " key= %s max= %f"%(key,np.max(a))
        return a 
    
    ''' apply artificial enhancement to R and B colors  '''
    def enhance(self):
        r = self.stamps["I"]  
        g = self.stamps["R"] 
        b = self.stamps["G"]

        r = r*(1 - self.enhLevel)
        b = b*(1 + self.enhLevel)
        r = np.where( r < 255,r, 255)
        g = np.where( g < 255,g, 255)
        b = np.where( b < 255,b, 255)

        
        self.stamps["I"] = r.astype(np.uint8)
        self.stamps["G"] = b.astype(np.uint8)
        self.stamps["R"] = g.astype(np.uint8)

    '''
         Perform processing of the input files creating png image stamp  
      '''  
    def processFile(self):

        "    "
#        print " zpFr= %f zpFg= %f zpFi= %f "%( self.zpF["R"],self.zpF["G"],self.zpF["I"]) 
        zpAv = (self.zpF["R"] + self.zpF["G"] + self.zpF["I"])/3
        self.Rcor = self.Rcor*math.pow(10, -0.4*(self.zpF["R"]-zpAv))
        self.Gcor = self.Gcor*math.pow(10, -0.4*(self.zpF["G"]-zpAv))
        self.Icor = self.Icor*math.pow(10, -0.4*(self.zpF["I"]-zpAv))
#        print " Rcor= %f Gcor= %f Icor= %f "%(self.Rcor,self.Gcor,self.Icor)
        
        self.impars['R'] = self.getStampPars(self.stamps['R'])
        self.delta['R'] = self.impars['R'][3]*self.sigmaMin
        self.impars['G'] = self.getStampPars(self.stamps['G'])
        self.delta['G'] = self.impars['G'][3]*self.sigmaMin
        self.impars['I'] = self.getStampPars(self.stamps['I'])
        self.delta['I'] = self.impars['I'][3]*self.sigmaMin
        self.stamps['R'] = self.filterIm('R',self.sigmaMin,self.sigmaMax,self.Rcor)
        self.stamps['G'] = self.filterIm('G',self.sigmaMin,self.sigmaMax,self.Gcor)
        self.stamps['I'] = self.filterIm('I',self.sigmaMin,self.sigmaMax,self.Icor)
        " now create scale image "
        Im = self.stamps["R"] + self.stamps["G"] + self.stamps["I"]
        d = Im*self.Q        
        d = np.where(d<=0.,0.0001*self.Q, d)
        c = Im*(self.Q*self.Alpha)
        F = np.arcsinh(c)
        X = F/d
        maxR = self.maxP['R']; maxG = self.maxP['G']; maxB = self.maxP['I']
        maxI = maxR + maxG + maxB
#        print " maxI=%f"%maxI
        scale = 765.0*self.Q/(math.asinh(maxI*self.Q*self.Alpha)*self.satCut)
#        print "Scale= %f"%scale

        self.stamps["R"] = self.rescale("R",X,scale)

        self.stamps["G"] = self.rescale("G",X,scale)

        self.stamps["I"] = self.rescale("I",X,scale)
        
        self.enhance()
        " now create PNG image " 
        rgbArray = np.zeros((self.nrow,self.ncol,3), 'uint8')
        rgbArray[..., 0] = self.stamps["I"]
        rgbArray[..., 1] = self.stamps["R"]
        rgbArray[..., 2] = self.stamps["G"]
        im = Image.fromarray(np.flipud(rgbArray))
#
# Now as we have the image let write metadata
#
        meta = PngImagePlugin.PngInfo()
        meta.add_text("crval1", str(self.crval1))
        meta.add_text("crval2", str(self.crval2))
        meta.add_text("crpix1", str(self.crpix1))
        meta.add_text("crpix2", str(self.crpix2))
        meta.add_text("cd1_1", str(self.cd1_1))
        meta.add_text("cd1_2", str(self.cd1_2))
        meta.add_text("cd2_1", str(self.cd2_1))
        meta.add_text("cd2_2", str(self.cd2_2))
        meta.add_text("copyright", "DES collaboration")
        reserved = ('interlace', 'gamma', 'dpi', 'transparency', 'aspect')
        f = open(self.outfile,"wc")
            # copy metadata into new object
        for k,v in im.info.iteritems():
            if k in reserved: continue
            meta.add_text(k, v, 0)
#
        im.save(f, "PNG", pnginfo=meta)
        f.close()
                    
    def produce(self,outfile):
        self.outfile = outfile
        zpAv = (self.zpF["R"] + self.zpF["G"] + self.zpF["I"])/3
        self.Rcor = self.Rcor*math.pow(10, -0.4*(self.zpF["R"]-zpAv))
        self.Gcor = self.Gcor*math.pow(10, -0.4*(self.zpF["G"]-zpAv))
        self.Icor = self.Icor*math.pow(10, -0.4*(self.zpF["I"]-zpAv))
#        print " Rcor= %f Gcor= %f Icor= %f "%(self.Rcor,self.Gcor,self.Icor)
        
        self.impars['R'] = self.getStampPars(self.stamps['R'])
        self.delta['R'] = self.impars['R'][3]*self.sigmaMin
        self.impars['G'] = self.getStampPars(self.stamps['G'])
        self.delta['G'] = self.impars['G'][3]*self.sigmaMin
        self.impars['I'] = self.getStampPars(self.stamps['I'])
        self.delta['I'] = self.impars['I'][3]*self.sigmaMin
        self.stamps['R'] = self.filterIm('R',self.sigmaMin,self.sigmaMax,self.Rcor)
        self.stamps['G'] = self.filterIm('G',self.sigmaMin,self.sigmaMax,self.Gcor)
        self.stamps['I'] = self.filterIm('I',self.sigmaMin,self.sigmaMax,self.Icor)
        " now create scale image "
        Im = self.stamps["R"] + self.stamps["G"] + self.stamps["I"]
        d = Im*self.Q        
        d = np.where(d<=0.,0.0001*self.Q, d)
        c = Im*(self.Q*self.Alpha)
        F = np.arcsinh(c)
        X = F/d
        maxR = self.maxP['R']; maxG = self.maxP['G']; maxB = self.maxP['I']
        maxI = maxR + maxG + maxB
#        print " maxI=%f"%maxI
        scale = 765.0*self.Q/(math.asinh(maxI*self.Q*self.Alpha)*self.satCut)
#        print "Scale= %f"%scale

        self.stamps["R"] = self.rescale("R",X,scale)

        self.stamps["G"] = self.rescale("G",X,scale)

        self.stamps["I"] = self.rescale("I",X,scale)
        
        self.enhance()
        " now create PNG image " 
        rgbArray = np.zeros((self.nrow,self.ncol,3), 'uint8')
        rgbArray[..., 0] = self.stamps["I"]
        rgbArray[..., 1] = self.stamps["R"]
        rgbArray[..., 2] = self.stamps["G"]
        im = Image.fromarray(np.flipud(rgbArray))
#
# Now as we have the image let write metadata
#
        meta = PngImagePlugin.PngInfo()
        meta.add_text("crval1", str(self.crval1))
        meta.add_text("crval2", str(self.crval2))
        meta.add_text("crpix1", str(self.crpix1))
        meta.add_text("crpix2", str(self.crpix2))
        meta.add_text("cd1_1", str(self.cd1_1))
        meta.add_text("cd1_2", str(self.cd1_2))
        meta.add_text("cd2_1", str(self.cd2_1))
        meta.add_text("cd2_2", str(self.cd2_2))
        meta.add_text("copyright", "DES collaboration")
        reserved = ('interlace', 'gamma', 'dpi', 'transparency', 'aspect')
        f = open(self.outfile,"wc")
            # copy metadata into new object
        for k,v in im.info.iteritems():
            if k in reserved: continue
            meta.add_text(k, v, 0)
#
        if self.tagFlg == 1:
#            print "Write tag"
            self.makeTag(im)
        im.save(f, "PNG", pnginfo=meta)

        f.close()            
                    
if __name__ == "__main__":
    print sys.argv
    outFile=""
    infile=""
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:",
                                  ["infile="])
    except getopt.GetoptError:
        print "MakeHumViImage.py -i <infile>" 
        sys.exit(2)
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "MakeHumViImage.py -i <infile> " 
            sys.exit()
        elif opt in ("-i","--infile"):
            print "got -i arg= %s"%arg
            infile=arg
    mrgbi = MakeHumViImage(infile)
    mrgbi.getInput()
    mrgbi.printInputList()
    mrgbi.processFile()
    sys.exit(0)
