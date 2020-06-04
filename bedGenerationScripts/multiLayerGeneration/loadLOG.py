import os
import sys
import logG

'''
xmax = float(sys.argv[1])
ymax = float(sys.argv[2])
zmax = float(sys.argv[3])
avg = float(sys.argv[4])
devn = float(sys.argv[5])
nobeds = int(sys.argv[6])
noB = int(sys.argv[7])
outfile = "particle_input.txt"

if len(sys.argv)>8: outfile = "particle_input_z"+str(zmax)+"_noB"+str(noB)+".txt" 
'''

outfile = "particle_input.txt"

totB = 10

mainDOut = os.getcwd()

pars = {'1.1':range(1,11), '1.2':range(1,11),'1.3':range(1,11), '1.4':range(1,11),'1.5':range(1,11), '1.6':range(1,11)}#{'0.9':range(1,21), '1':range(1,21)}

#pars = {'0.8':range(11,21)}#pars = {'0.8':range(1,11)}

zlow = 0.7 #should be left out for single layers, but is present for multiple layers to shift the particles up

for jk in pars.keys():
	zmax = eval(jk)

	#Make new directory for the z value and navigate to that directory
	dirstr = 'zmax_'+jk
	os.makedirs(dirstr)
	os.chdir(dirstr)
	mainDIn = os.getcwd()

	for noB in pars[jk]:
		dirstr = 'noB_'+str(noB)
		os.makedirs(dirstr)
		os.chdir(dirstr)
		mainDInE = os.getcwd()

		fileO = open('ParticleInfo.log','w')

		for bed in range(totB):
			dirstr = 'bedgen_'+str(bed)
			os.makedirs(dirstr)
			os.chdir(dirstr)

			noP = logG.logen(zmax+zlow,noB,outfile = outfile,zlow = zlow)

			os.chdir(mainDInE)

			fileO.write('%d\n' %noP)

			print ("zmax_{} noB_{} bed_{}".format(zmax,noB,bed))

		fileO.close()
		os.chdir(mainDIn)
	os.chdir(mainDOut)



