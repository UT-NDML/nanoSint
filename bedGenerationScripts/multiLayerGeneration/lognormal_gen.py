import numpy as np
import random as rand
import math
import sys

xmax = float(sys.argv[1])
ymax = float(sys.argv[2])
zmax = float(sys.argv[3])
avg = float(sys.argv[4])
devn = float(sys.argv[5])
nobeds = int(sys.argv[6])
noB = int(sys.argv[7])

outfile = "particle_input.txt"

if len(sys.argv)>8: outfile = "particle_input_z"+str(zmax)+"_noB"+str(noB)+".txt" 

upB = 10000

#noB = 25

maxv = 50


for ip in range(nobeds):

	particles = []
	moreR = 0

	while (moreR < noB): 

		moreP = 1	
		ctw = 0

		while (moreP==1 or len(particles)<=20):

			if (avg<0.0000001 and devn<0.0000001):
				radN = np.random.normal(5.367,0.3987)#default for intrinsiq glass ink
			else:
				radN = np.random.normal(avg,devn)

			rad = math.exp(radN)/2000

			for j in range(upB):

				x = rand.uniform(rad,xmax-rad)
				y = rand.uniform(rad,ymax-rad)
				z = rand.uniform(rad,zmax-rad)

				for part in particles:
					distApart = math.sqrt((x-part[0])**2+(y-part[1])**2+(z-part[2])**2)
					if (distApart<=rad+part[3]):
						break
				else:
					particles.append([x,y,z,rad])
					break

				if j==upB-1:
					moreP=0

		moreR += 1

	print "%d: %d particles" %(ip,len(particles))

	outf = open(outfile,'w')

	for part in particles:
		outf.write('%.9f\t%.9f\t%.9f\t%.6f\t0.00896\t%f\t%f\t%f\t0\t0\t0\n' %(part[0],part[1],part[2],part[3],rand.uniform(-maxv,maxv),rand.uniform(-maxv,maxv),rand.uniform(-maxv,maxv)))

	outf.close()
