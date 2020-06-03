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

#hs = raw_input('Height: ')
#nb = raw_input('Serial Number: ')
#maxnb = raw_input('Maximum file Number: ')
#distnb = raw_input('Distance betwen files: ')
#procnb = raw_input('Number of processors: ')

#xsize = int(hs)
#SNo = int(nb)
#maxfile = int(maxnb)
#filed = int(distnb)
#sizep = int(procnb)

zsize = int(sys.argv[1])
ysize = int(sys.argv[2])
xsize = int(sys.argv[3])
SNo = int(sys.argv[4])
minfile = int(sys.argv[5])
maxfile = int(sys.argv[6])
filed = int(sys.argv[7])
sizep = int(sys.argv[8])

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

	dirstd = str(os.getcwd())
        
	global maxfile
	global filed

	boxSize = 95

	#xsize = 51
	#ysize = 192
	#zsize = 192
		

	xvec0 = [25,28,30,32,35,38,58,60,62,63,65,70,40,45,47,48,50,55]
	yvec0 = [25,28,30,32,35,38,58,60,62,63,65,70,40,45,47,48,50,55]

        xvec1 = [25,28,30,32,35,38,40,45,47,48,50,58,55,60,62,63,65,70]
        yvec1 = [25,28,30,32,35,38,40,45,47,48,50,58,55,60,62,63,65,70] 

	xvec2 = [23,24,25,24,25,23]
	yvec2 = [18,19,19,18,18,19]

	LIMS = ['bed5']#['Newbedsdp6/bed5Code4_z52/cutOffR/R00007'] #['Newbedsdp6/bed5SN2030', 'Newbedsdp6/bed5SN2038', 'Newbedsdp6/bed5SN2054', 'Newbedsdp6/bed5SN2081','Newbedsdp6/bed5SN2084']# ['Newbedsdp6/bed5SN2084'] #['1by1bedsUncertainty/2019Calibration/SNsto5000/SN2030','1by1bedsUncertainty/2019Calibration/SNsto5000/SN2038']  #['bed1','bed3','bed5', 'bed6', 'bed11','bed14','bed15','bed16','bed19', 'Newbedsdp6/bed0', 'Newbedsdp6/bed1', 'Newbedsdp6/bed4', 'Newbedsdp6/bed5']
	
	Xprop = [192]#, 192, 192, 192, 192]#[192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 190, 192]
	Yprop = [192]#, 192, 192, 192, 192]#[192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 196, 192]
	Zprop = [66]#, 52, 52, 52, 52]#[53, 54, 66, 50, 56, 57, 59, 45, 48, 52, 52, 49, 52]

	SNoA = [504]#[50400007]#[2030, 2038, 2054, 2081, 2084]#[50400002, 50400002, 50400002, 50400002, 50400002, 50400002, 50400002, 50400002, 50400002, 50400002, 50400002, 50400002, 50400002]  #[2084]

	procSize = len(xvec0)

	Ct1 = procSize*len(LIMS)/sizep
	Ct2 = procSize*len(LIMS) - Ct1*sizep	
	if num < Ct2: numCt = Ct1+1
	else: numCt = Ct1

	for numctr in range(numCt):

		numSH = num + numctr*sizep
		limc = numSH/procSize
		lims = numSH%procSize

                lim = 0

		limb = LIMS[limc]
		zsize = Xprop[limc]
		ysize = Yprop[limc]
		xsize = Zprop[limc]
		SNo = SNoA[limc]

        	newdr = dirstd+'/'+limb

        	os.chdir(newdr)

		if limb=='bed6':
			ylb = yvec1[lims]
			zlb = xvec1[lims]
                        yhb = ylb + boxSize
                        zhb = zlb + boxSize
			strOpen = "SN%ddatafull2by2square_xl%d_xh%d_lim_%s.txt" %(SNo,zlb,zhb,str(lim))
		elif limb[0]=='1':
                        ylb = yvec2[lims]
                        zlb = xvec2[lims]
			lim = 0 #0.15
			boxSize = 56
                        yhb = ylb + boxSize
                        zhb = zlb + boxSize
			strOpen = "SN%ddatafull_x%d_y%d_lim_%s.txt" %(SNo,zlb,ylb,str(lim))			
		else:
                	ylb = yvec0[lims]
                	zlb = xvec0[lims]
                	yhb = ylb + boxSize
                	zhb = zlb + boxSize		
			strOpen = "SN%ddatafull2by2square_xl%d_xh%d_lim_%s.txt" %(SNo,zlb,zhb,str(lim))

		#limv = 0.025
		#lim = 0.15
		fanl = open(strOpen,"w")
		fanl.write("filetimex100\ttotal count\tmaxz\tdensity\tshrinkage\tpixel count relative density\tempty pixels\tpore density\ttotal sum of rho\tsum in bounds\tminimum z\tvolume\tpixel size relative density\n")

		rawsimt = []
		fullsimd = []
		
		ctr = 0
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
				volume = (maxz-min(z))*(yhb-ylb)*(zhb-zlb)
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

		print "Start:%d, %d, %s done Analsis\nStart Calibration" %(zlb,ylb,str(limb))

		if False:
			startlow = 100
			starthigh = 100000
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
			print

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
						err1[i] = abs(actd(float(simt[i]-simt[0])/A1,temp) - simd[i])
						err2[i] = abs(actd(float(simt[i]-simt[0])/A2,temp) - simd[i])

					sumerr1 = 0
					sumerr2 = 0
					for i in range(0,len(err1)):
						sumerr1 = sumerr1 + (err1[i]**0.005)
						sumerr2 = sumerr2 + (err2[i]**0.005)

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
			
			print "Start:%d, %d done Calibration" %(zlb,ylb)
			
		os.chdir(dirstd)
	return

jobs = []
for i in range(sizep):
	p = mp.Process(target=worker, args=(i,))
	jobs.append(p)
	p.start()
