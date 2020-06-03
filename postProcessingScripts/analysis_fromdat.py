import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.animation as animation
from pylab import *
import struct
import math
import os
import random

maxfile = 10000
filed = 50

#SNOall = [50400002]
SNOall = [2,3,7,27,28,29,295]
TimeC = [150,150,300,400,450,550,550]

for SNoc in range(len(SNOall)):

        SNo = int('5040000'+str(SNOall[SNoc]))

        boxSize = 95

        xsize = 52
        ysize = 192
        zsize = 192

        ylb = 47
        zlb = 47
        
        yhb = ylb + boxSize
        zhb = zlb + boxSize

        limv = 0.025
        lim = 0 #0.15
        fanl = open("SN%ddatafullbox_x%d_y%d.txt" %(SNo,zlb,ylb),"w")
        fanl.write("filetimex100\ttotal count\tmaxz\tdensity\tshrinkage\tpixel count relative density\tempty pixels\tpore density\ttotal sum of rho\tsum in bounds\tminimum z\tvolume\tpixel size relative density\n")
        
        ctr = 0
        
        print "SNo: %d, Start:%d, %d, start Analysis" %(SNo,zlb,ylb)
        maxfile=TimeC[SNoc]
        for ai in range(0,maxfile+filed,filed):

                bnf = open("rho%dSN%d.dat" %(ai,SNo),"rb")
                fcont = bnf.read()
                #print len(fcont)
                a = len(fcont)/8
                Vd = struct.unpack('d'*a,fcont[:a*8])

                size = len(Vd)

                Zd = [0 for x in range(size)]
                # print "length of an array is: %d" % size

                cd = np.zeros(size)

                sumvt = 0
                counti = 0

                for z in range(zsize):
                        for y in range(ysize):
                                for x in range(xsize):
                                        sumvt = sumvt + Vd[counti]
                                        Zd[counti] = x
                                        counti = counti+1

                xlb = 0 
                xhb = xsize

                countl = 0
                totalb = 0
                sumv = 0

                for c in range(zlb,zhb):
                        for b in range(ylb,yhb):
                                for a in range(xlb,xhb):
                                        k = a + xsize*b + xsize*ysize*c
                                        if Vd[k] > lim:
                                                countl = countl + 1
                                                cd[totalb] = Zd[k]
                                                sumv = sumv + Vd[k]
                                        else:
                                                countl = countl + 0
                                                cd[totalb] = 1000.0
                                        totalb = totalb + 1

                val = np.array([1000.0])

                z = np.setdiff1d(cd,val,assume_unique = True)

                if ai == 0:
                        Init_maxz = max(z)
                        Init_density = countl
                        Init_empty = totalb - countl
                        maxz = max(z)
                        shrinkage = Init_maxz - maxz
                        density = countl
                        rel_density_pix = float(density)/float(totalb)
                        empty_pix = totalb - countl
                        pore_density = float(empty_pix)/float(totalb)
                        volume = maxz*(yhb-ylb)*(zhb-zlb)
                        rel_density_size = sumv/float(volume)
                        sumv0 = sumv
                else:
                        maxz = max(z)
                        shrinkage = Init_maxz - maxz
                        density = countl
                        rel_density_pix = float(density)/float(totalb)
                        empty_pix = totalb - countl
                        pore_density = float(empty_pix)/float(totalb)
                        volume = (maxz-min(z))*(yhb-ylb)*(zhb-zlb)
                        rel_density_size = sumv/float(volume)

                fanl.write("%d\t%d\t%d\t%d\t%d\t%f\t%d\t%f\t%f\t%f\t%d\t%d\t%f\t%f\n" %(ai,totalb,maxz,density,shrinkage,rel_density_pix,empty_pix,pore_density,sumvt,sumv,min(z),volume,rel_density_size,(sumv-sumv0)/sumv))
                
                #print "%d:%d done" %(num, ai)

        print "SNo: %d, Start:%d, %d, done Analysis, start Calibration" %(SNo,zlb,ylb)

        fanl.close()
