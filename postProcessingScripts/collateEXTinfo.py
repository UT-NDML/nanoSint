import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.animation as animation
from pylab import *
import shutil 
import struct
import os

startbeds = 0
nobeds = 18

dirstd = str(os.getcwd())

fill = open("CollateExtractedInfo.txt",'w')
fill.write("boxno\tminCalib\tmaxCalib\tAverageCalibration\tminError\tmaxError\tAveragePercentageError\n")

for bedno in range(startbeds,nobeds):

    calibfac = []
    percerr = []

    newdr = dirstd+'/bed'+str(bedno)
    os.chdir(newdr)

    print str(os.getcwd())
    
    if os.path.exists("Extractbed.txt"):

        xtf = open("Extractbed.txt",'r')

        for lin in xtf.readlines():
            rows = lin.split(" ")

            calibfac.append(float(rows[2]))

            fin = rows[len(rows)-1].split("\n")

            percerr.append(float(fin[0]))
                
        os.chdir(dirstd)

        fill.write("%d\t%f\t%f\t%f\t%f\t%f\t%f\n" %(bedno,min(calibfac),max(calibfac),mean(calibfac),min(percerr),max(percerr),mean(percerr)))

    os.chdir(dirstd)

fill.close()
