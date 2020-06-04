import multiprocessing as mp
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import struct
import math

startp = 800
maxp = 1000
fild = 100
heightmin = 400
heightmax = 600

errTol = 0.001

convfac = 0.010584

fanl = open('bedInfo.txt','w')

fanl.write("Timefile\tdensity(kg/m^3)\tnumber of touching particle\toverlap(nm)\theight(nm)\n")

infoMat = [[0 for x in range(4)] for y in range(3)]
ctful = 0

for aa in range(startp,maxp+fild,fild):

        filenam = "particleOut%d.txt" %aa
        fild = open(filenam,'r')

        ctr = 0

        fileAr = fild.readlines()
        
        nop = len(fileAr)
        
        ctr = 0

        X = [0 for x in range(nop)]
        Y = [0 for x in range(nop)]
        Z = [0 for x in range(nop)]
        R = [0 for x in range(nop)]

        Xo = [0 for x in range(nop)]
        Yo = [0 for x in range(nop)]
        Zo = [0 for x in range(nop)]
        Ro = [0 for x in range(nop)]

        Xmin = [0 for x in range(3)]
      
        for line in fileAr:
                cols = line.split('\t')  

                Xc = float(cols[0])#[line.split('\t')[0] for line in fild]
                Yc = float(cols[1])#[line.split('\t')[1] for line in fild]
                Zc = float(cols[2])#[line.split('\t')[2] for line in fild]
                Rc = float(cols[3])

                Xo[ctr] = Xc
                Yo[ctr] = Yc
                Zo[ctr] = Zc
                Ro[ctr] = Rc
                
                X[ctr] = int(round(Xc/convfac))
                Y[ctr] = int(round(Yc/convfac))
                Z[ctr] = int(round(Zc/convfac))
                R[ctr] = int(round(Rc/convfac))
                
                ctr = ctr+1
                
        maxx = 0
        maxy = 0
        maxz = 0
        minx = 10
        miny = 10
        minz = 10

        nop = ctr
        
        for ai in range(nop):
                if maxx < X[ai]+R[ai]:              
                        maxx = X[ai]+R[ai]
                if maxy < Y[ai]+R[ai]: 
                        maxy = Y[ai]+R[ai]
                if maxz < Z[ai]+R[ai]: 
                        maxz = Z[ai]+R[ai]
                if minx > X[ai]-R[ai]:              
                        minx = X[ai]-R[ai]
                if miny > Y[ai]-R[ai]: 
                        miny = Y[ai]-R[ai]
                if minz > Z[ai]-R[ai]: 
                        minz = Z[ai]-R[ai]

        boxS = max(maxx,maxy,maxz)+10

        boxSizex = boxS
        boxSizey = boxS
        boxSizez = boxS

        collid = 0
        sumlid = 0.0

        for ai in range(ctr):
                for bi in range(ai+1,ctr):
                        Xmin[0] = Xo[ai]-Xo[bi]
                        Xmin[1] = Yo[ai]-Yo[bi]
                        Xmin[2] = Zo[ai]-Zo[bi]

                        if sqrt(Xmin[0]**2+Xmin[1]**2+Xmin[2]**2)<Ro[ai]+Ro[bi]:
                                #print "Touching particles %d + %d" %(ai,bi)
                                #print "overlap = %f nm" %(1000*(Ro[ai] + Ro[bi] - sqrt(Xmin[0]**2+Xmin[1]**2+Xmin[2]**2)))

                                #print "%f %f %f %f %f" %(Xmin[0], Xmin[1], Xmin[2], R[ai], R[bi])
                                sumlid = sumlid + 1000*(Ro[ai] + Ro[bi] - sqrt(Xmin[0]**2+Xmin[1]**2+Xmin[2]**2))
                                collid = collid+1

        #print "Number of touching particles:> %d" %collid
        #print "Average overlap = %f nm" %(sumlid/(collid))

        countpix = 0

        boxX = [minx, maxx]
        boxY = [miny, maxy]
        boxZ = [minz, maxz]

        Size = (boxX[1])*(boxY[1])*(boxZ[1])

        for ctr in range(nop):
                rad = R[ctr]
                xa = X[ctr]
                ya = Y[ctr]
                za = Z[ctr]
                
                for zz in range(boxZ[0],boxZ[1]):
                        for yy in range(boxY[0],boxY[1]):
                                for xx in range(boxX[0],boxX[1]):
                                        if ((zz-za)**2+(yy-ya)**2+(xx-xa)**2)<=rad**2:
                                                countpix = countpix+1

        countemp = Size-countpix

        infoMat[ctful][0] = float(countpix)*0.00896/Size
        infoMat[ctful][1] = sumlid/collid	
	infoMat[ctful][2] = max(Zo)
        infoMat[ctful][3] = (max(Zo) + Ro[Zo.index(max(Zo))])*1000
        
        fanl.write("%d\t%f\t%d\t%f\t%f\n" %(aa,infoMat[ctful][0],collid,infoMat[ctful][1],infoMat[ctful][3]))
        
        #print "done=> %d" %aa
        ctful = ctful + 1

fanl.close()

errorprop = 0

for ii in range(3):
        for jj in range(3):
                if (ii!=jj):
                        errorprop = errorprop+abs(infoMat[ii][0]-infoMat[jj][0])
                        errorprop = errorprop+abs(infoMat[ii][1]-infoMat[jj][1])
                        errorprop = errorprop+abs(infoMat[ii][2]-infoMat[jj][2])

if (errorprop < errTol and infoMat[-1][-1] >= heightmin and infoMat[-1][-1]<=heightmax):
	fanlb = open('bed.txt','w')
	fangu = open('bedguide.txt','w')

        aa = 1000
        filenam = "particleOut%d.txt" %aa
        fild = open(filenam,'r')

        ctr = 0
        for line in fild.readlines():
                cols = line.split('\t')  
                
                X[ctr] = int(round(float(cols[0])/convfac))+2
                Y[ctr] = int(round(float(cols[1])/convfac))+2
                Z[ctr] = int(round(float(cols[2])/convfac))+2
                R[ctr] = int(round(float(cols[3])/convfac))

                fanlb.write("%d\t%d\t%d\t%d\n" %(R[ctr],X[ctr],Y[ctr],Z[ctr]))
                ctr = ctr+1

        maxz = 0

        nop = ctr
        
        for ai in range(nop):
                if maxz < Z[ai]+R[ai]: 
                        maxz = Z[ai]+R[ai]

        fangu.write("%d\t%d" %(nop,maxz))

	fanlb.close()
	fangu.close()
