import multiprocessing as mp
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import struct
import math
import sys

#zs = raw_input('Height: ') #Z height
#nb = raw_input('Serial Number: ') #Serial number relating constants
#minnb = raw_input('Minimum file Number: ') #minimum file to read
#maxnb = raw_input('Maximum file Number: ') #maximum file to read
#distnb = raw_input('Distance betwen files: ') #distance between files
#procnb = raw_input('Number of processors: ') #number of processors
#runnb = raw_input('To run just whole type 1, to run just bounds type 2, to run whole and bounds type 3: ') #defining if the simulation should run the full bed or within set bounds used for analysis

#xsize = int(zs)
#SNo = int(nb)
#minfile = int(minnb)
#maxfile = int(maxnb)
#filed = int(distnb)
#sizep = int(procnb)
#simzoom = int(runnb)

zsize = int(sys.argv[1])
ysize = int(sys.argv[2])
xsize = int(sys.argv[3])
SNo = int(sys.argv[4])
minfile = int(sys.argv[5])
maxfile = int(sys.argv[6])
filed = int(sys.argv[7])
sizep = int(sys.argv[8])
simzoom = int(sys.argv[9])

wholeruns = int(math.floor((maxfile + filed - minfile)/filed))
left = int(math.floor(wholeruns/sizep))

print wholeruns
print left

def worker (num):
    global sizep
    global left
    global filed
    for zt in range(0,left+1):
        for ang in range(0,1):

	    if (zt < left):
            	ai = num*(filed) + zt*filed*sizep + minfile
	    else:
		if (num < (wholeruns - (sizep*left))):
		    ai = num*(filed) + zt*filed*sizep + minfile
	  	else:
		    exit()

            if ang%2==0:
                elev = 40 #70
            else:
                elev = 70

            bnf = open("rho%dSN%d.dat" %(ai,SNo),"rb")
            fcont = bnf.read()
            size = len(fcont)/8
            Vd = struct.unpack('d'*size,fcont[:size*8])

            Xd = [0 for x in range(size)]
            Yd = [0 for x in range(size)]
            Zd = [0 for x in range(size)]

            #xsize = 51
            #ysize = 192
            #zsize = 192

            ysizep = 100 #box size for plot if not ysize
            ysizep2 = 150 #box size for plot if not ysize		

            counti = 0
            for z in range(zsize):
                for y in range(ysize):
                    for x in range(xsize):
                        Zd[counti] = x
                        Yd[counti] = y
                        Xd[counti] = z
                        counti = counti+1

	    lim = 0.1

	    if (simzoom == 1): #whole

		countl = 0

		for i in range(0,size):
		    if Vd[i] > lim:
		        countl = countl+1
		    else:
		        countl = countl+0

		x = np.zeros(countl)
		y = np.zeros(countl)
		z = np.zeros(countl)
		v = np.zeros(countl)
		totalb = 0

		for c in range(zsize):
		    for b in range(ysize):
		        for a in range(xsize):
			    k = a + xsize*b + xsize*ysize*c
			    if Vd[k] > lim:
			        x[totalb] = Xd[k]
			        y[totalb] = Yd[k]
			        z[totalb] = Zd[k]
			        v[totalb] = Vd[k]
			        totalb = totalb + 1

		#fig = plt.figure(frameon=False)
		fig = plt.figure()
		ax = Axes3D(fig)
		#ax.axis('off')

		ax.set_xlim3d([50,ysizep2])
		ax.set_ylim3d([50,ysizep2])
		ax.set_zlim3d([0,ysizep])
		ax.scatter(x,y,z,s=100*v,c='white',marker='o',edgecolor='purple')
		ax.view_init(elev,azim=50)

		fig.savefig("rho%dSN%dfull.png" %(ai,SNo))

                plt.close('all')
		print "%d:%d whole done" %(num,ai)

	    elif (simzoom == 2): #bounds

		xlb = 0
		xhb = xsize
		ylb = 47
		yhb = 142
		zlb = 47
		zhb = 142

		totalb = 0

		rdsize = (zhb-zlb)*(yhb-ylb)*(xhb-xlb)
		xrd = np.zeros(rdsize)
		yrd = np.zeros(rdsize)
		zrd = np.zeros(rdsize)
		vrd = np.zeros(rdsize)

		for c in range(zlb,zhb):
		    for b in range(ylb,yhb):
		        for a in range(xlb,xhb):
		            k = a + xsize*b + xsize*ysize*c
		            if Vd[k] > lim:
		                xrd[totalb] = Xd[k]
		                yrd[totalb] = Yd[k]
		                zrd[totalb] = Zd[k]
		                vrd[totalb] = Vd[k]
		                totalb = totalb + 1

		#fig = plt.figure(frameon=False)
		fig = plt.figure()
		ax = Axes3D(fig)
		#ax.axis('off')

		ax.set_xlim3d([50,ysizep2])
		ax.set_ylim3d([50,ysizep2])
		ax.set_zlim3d([0,ysizep])
		ax.scatter(xrd,yrd,zrd,s=100*vrd,c='cyan',marker='o',edgecolor='black')
		ax.view_init(elev,azim=50)

		fig.savefig("rho%dSN%dbdd.png" %(ai,SNo))

                plt.close('all')
		print "%d:%d bounds done" %(num,ai)

	    else: #in and out

		countl = 0

		for i in range(0,size):
		    if Vd[i] > lim:
		        countl = countl+1
		    else:
		        countl = countl+0

		x = np.zeros(countl)
		y = np.zeros(countl)
		z = np.zeros(countl)
		v = np.zeros(countl)
		totalb = 0

		for c in range(zsize):
		    for b in range(ysize):
		        for a in range(xsize):
			    k = a + xsize*b + xsize*ysize*c
			    if Vd[k] > lim:
			        x[totalb] = Xd[k]
			        y[totalb] = Yd[k]
			        z[totalb] = Zd[k]
			        v[totalb] = Vd[k]
			        totalb = totalb + 1

		xlb = 0
		xhb = xsize
		ylb = 47
		yhb = 142
		zlb = 47
		zhb = 142

		totalb = 0

		rdsize = (zhb-zlb)*(yhb-ylb)*(xhb-xlb)
		xrd = np.zeros(rdsize)
		yrd = np.zeros(rdsize)
		zrd = np.zeros(rdsize)
		vrd = np.zeros(rdsize)

		for c in range(zlb,zhb):
		    for b in range(ylb,yhb):
		        for a in range(xlb,xhb):
		            k = a + xsize*b + xsize*ysize*c
		            if Vd[k] > lim:
		                xrd[totalb] = Xd[k]
		                yrd[totalb] = Yd[k]
		                zrd[totalb] = Zd[k]
		                vrd[totalb] = Vd[k]
		            	totalb = totalb + 1

		for zm in range(0,2):
		    if (zm==0): #whole	 
			#fig = plt.figure(frameon=False)
			fig = plt.figure()
			ax = Axes3D(fig)
			#ax.axis('off')

			ax.set_xlim3d([50,ysizep2])
			ax.set_ylim3d([50,ysizep2])
			ax.set_zlim3d([0,ysizep])
			ax.scatter(x,y,z,s=100*v,c='white',marker='o',edgecolor='purple')
			ax.view_init(elev,azim=50)

			fig.savefig("rho%dSN%dfull2.png" %(ai,SNo))
			print "%d:%d whole done" %(num,ai)

		    elif (zm == 1): #bounds 
			#fig = plt.figure(frameon=False)
			fig = plt.figure()
			ax = Axes3D(fig)
			#ax.axis('off')

			ax.set_xlim3d([50,ysizep2])
			ax.set_ylim3d([50,ysizep2])
			ax.set_zlim3d([0,ysizep])
			ax.scatter(xrd,yrd,zrd,s=100*vrd,c='cyan',marker='o',edgecolor='black')
			ax.view_init(elev,azim=50)

			fig.savefig("rho%dSN%dbdd2.png" %(ai,SNo))
			print "%d:%d bounds done" %(num,ai)

                    plt.close('all')
		    print "%d:%d done" %(num,ai)
    return

"""
if __name__== '__main__':
    jobs = []
    for i in range(size):
        p = mp.Process(target=worker, args=(i,))
        jobs.append(p)
        p.start()
"""

jobs = []
for i in range(sizep):
    p = mp.Process(target=worker, args=(i,))
    jobs.append(p)
    p.start()
