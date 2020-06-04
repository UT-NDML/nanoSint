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

nb = 50400002 #raw_input('Serial Number: ')
maxnb = 40800 #raw_input('Maximum file Number: ')
distnb = 50 #raw_input('Distance betwen files: ')
procnb = 48 #raw_input('Number of processors: ')

SNo = int(nb)
maxfile = int(maxnb)
filed = int(distnb)
sizep = int(procnb)

def actd(t,temp):			

	if temp == 450:
		T = [0.1572,1.2099,0.2854,0.8681,-0.0053]

	elif temp == 500:
		T = [0.2354,0.5326,0.8508,1.212,-0.4228]

	elif temp == 550:
		T = [0.2543,0.2889,0.5528,1.2493,-0.5253]
		
	elif temp == 600:
		T = [0.1651,0.4672,0.1331,0.8051,-0.011]

	valu = T[0]*exp((-T[1]/(t+T[2]))+T[3])+T[4]

	return valu

def minID(Col):
	valm = min(Col)
	i = 0
	while i < len(Col):
		if (Col[i] == valm):
			vali = i
			i = len(Col)
		else:
			i = i + 1

	return (valm,vali)

def worker (num):

        startbeds = 0
        endbeds = 18
        
	global maxfile
	global filed

        height = []

        with open('TACCguide.log','r') as lin:
                for dks in lin:
                        height.append(10*int(dks[3])+int(dks[4])+2)

        boxfl = open('BoxSize.txt','r')
        boxAll = []

        for line in boxfl.readlines():
            cols = line.split('\n')
            boxAll.append(int(cols[0]))

        Allbds = open('AllBounds.txt','r')
        dirtogo = []
	xAll = []
	yAll = []

	for line in Allbds.readlines():
		cols = line.split(' ')
	
		dirtogo.append(int(cols[0]))
		xAll.append(int(cols[1]))

		colsu = cols[2].split('\r')
		yAll.append(int(colsu[0]))
	
        dirstd = str(os.getcwd())

	#defining the necessary start and stop conditions for processor

	Nobds = len(dirtogo) #total number of bounds cobination to search
	allmax = Nobds/sizep #should floor
	
	left = Nobds - allmax*sizep #leftovers

	#calculate the maximum number of files each processor has to analyze
	if num < left:
		procmax = allmax + 1
	else:
		procmax = allmax

	pctr = 0 #initiate processor counter
	
	#print "proc: %d all: %d left: %d procmax: %d" %(num,allmax,left,procmax)

        for pno in range(procmax):

		pctr = pno*sizep + num

		bedno = dirtogo[pctr]

                os.chdir(dirstd)
                newdr = dirstd+'/bed'+str(bedno)
                os.chdir(newdr)

		boxSize = boxAll[bedno]

		xsize = height[bedno]
		ysize = 104
		zsize = 104

		ylb = yAll[pctr]
		zlb = xAll[pctr]
		
		yhb = ylb + boxSize
		zhb = zlb + boxSize

		limv = 0.025
		lim = 0.15
		fanl = open("SN%ddatafullbox_550CboundOUT_x%d_y%d.txt" %(SNo,zlb,ylb),"w")
		fanl.write("filetimex100\ttotal count\tmaxz\tdensity\tshrinkage\tpixel count relative density\tempty pixels\tpore density\ttotal sum of rho\tsum in bounds\tminimum z\tvolume\tpixel size relative density\n")

		rawsimt = []
		fullsimd = []
		
		ctr = 0

                print "Processor: %d, File line: %d, Start:%d, %d, dir: %s, start Analysis" %(num,pctr,zlb,ylb,str(os.getcwd()))

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
			else:
				maxz = max(z)
				shrinkage = Init_maxz - maxz
				density = countl
				rel_density_pix = float(density)/float(totalb)
				empty_pix = totalb - countl
				pore_density = float(empty_pix)/float(totalb)
				volume = (maxz-min(z))*(yhb-ylb)*(zhb-zlb)
				rel_density_size = sumv/float(volume)

			fanl.write("%d\t%d\t%d\t%d\t%d\t%f\t%d\t%f\t%f\t%f\t%d\t%d\t%f\n" %(ai,totalb,maxz,density,shrinkage,rel_density_pix,empty_pix,pore_density,sumvt,sumv,min(z),volume,rel_density_size))

			rawsimt.append(ai)
			fullsimd.append(sumv)
			
			#print "%d:%d done" %(num, ai)

                print "Processor: %d, File line: %d, Start:%d, %d done Analysis, start Calibration" %(num,pctr,zlb,ylb)

		startlow = 100
		starthigh = 10000
		maxstep = 0.005

		A1 = startlow
		A2 = starthigh
		Sstep = 100

		temp = 550

		initerr = [0 for x in range(len(rawsimt))]
		rawsimd = [0 for x in range(len(rawsimt))]

		noT = len(rawsimt)

		for i in range(0,noT):
			rawsimd[i] = float(fullsimd[i] - fullsimd[0])/float(fullsimd[i])

		for  i in range(0,noT):
			initerr[i] = abs(rawsimd[i]-actd(0,temp))

		(inY,inI) = minID(initerr)

		mintime = rawsimt[inI]

		#print "Closest fit to initial experiments value is: %f\nWith error value: %f\nAt simulation time: %d" %(rawsimd[inI],inY,mintime)

		simt = [0 for x in range(len(rawsimt))]
		simd = [0 for x in range(len(rawsimt))]

		for i in range(0,len(simt)):
			simt[i] = rawsimt[i]
			simd[i] = rawsimd[i]

		err1 = [0 for x in range(len(simt))]
		err2 = [0 for x in range(len(simt))]

		while (Sstep > maxstep):
			Ebr = {}
			ct = 0
			while ( A2-A1 > 0):

				for i in range(0,len(simt)):

					act1 = actd(float(simt[i]-simt[0])/A1,temp)
					act2 = actd(float(simt[i]-simt[0])/A2,temp)

					err1[i] = abs(act1 - simd[i])**2
					err2[i] = abs(act2 - simd[i])**2

				sumerr1 = np.mean(err1)
				sumerr2 = np.mean(err2)

				#for i in range(0,len(err1)):
					#sumerr1 = sumerr1 + (err1[i]**0.005)
					#sumerr2 = sumerr2 + (err2[i]**0.005)

				SSE1 = sumerr1
				SSE2 = sumerr2

				Ebr[ct,1] = A1
				Ebr[ct,2] = SSE1
				Ebr[ct,3] = A2
				Ebr[ct,4] = SSE2

				A1 = A1 + Sstep
				A2 = A2 - Sstep
				ct = ct + 1

			Ebr2 = [0 for x in range(ct)]
			Ebr4 = [0 for x in range(ct)]

			for i in range(0,ct):
				Ebr2[i] = Ebr[i,2]
				Ebr4[i] = Ebr[i,4]

			if (min(Ebr2) < min(Ebr4)):
				ad = 1
				helpS = Ebr2
			else:
				ad = 3
				helpS = Ebr4

			SI = np.argsort(helpS)
			SY = sorted(helpS)

			if (Ebr[SI[1],ad] < Ebr[SI[2],ad]):
				A1 = Ebr[SI[1],ad]
				A2 = Ebr[SI[2],ad]
			else:
				A1 = Ebr[SI[2],ad]
				A2 = Ebr[SI[1],ad]


			Sstep = float(A2 - A1)/10

		av1 = np.mean(err1);
		max1 = max(err1);
		av2 = np.mean(err2);
		max2 = max(err2);

		#Finding the highest maximum and average error between the simulation and the data

		if (av1 > av2):
			average = av1
		else: 
			average = av2

		if (max1 > max2):
			maximum = max1
		else: 
			maximum = max2

		(mY,mI) = minID(helpS)
		scalingfactor = Ebr[mI,ad]

		#print "Minimum Squared Sum of Errors: %f\nBest Scaling factor for time: %f\nAverage error: %f\nMaximum error: %f\n" %(mY,scalingfactor,average,maximum)

		actt = [0 for x in range(len(simt))]
		expd = [0 for x in range(len(simt))]
		perc = [0 for x in range(len(simt))]

		for i in range(0,len(simt)):
			actt[i] = float(simt[i]-simt[0])/Ebr[mI,ad]

		for i in range(0,len(simt)):
			expd[i] = actd(actt[i],temp)
			perc[i] = abs(expd[i]-simd[i])*100/expd[i]

		fanl.write("%d\t%d\t%f\t%f\t%f\t%f\n" %(xlb,ylb,scalingfactor,average,maximum,np.mean(perc)))

		fanl.close()
		
		print "Processor: %d, File line: %d, Start:%d, %d done Calibration" %(num,pctr,zlb,ylb)

	return

jobs = []
for i in range(sizep):
	p = mp.Process(target=worker, args=(i,))
	jobs.append(p)
	p.start()
