import numpy as np
import struct

#bulk value of density
bulkR = 8960 #kg/m^3

#simulation bed properties
xsize = 61
ysize = 104
zsize = 104
noparts = 21

#rho limit for analysis
lim = 0

#Serial number for constants
SNo = 50400002

#analysis box properties
xlb = 22
ylb = 22
boxSize = 55
xhb = xlb + boxSize
yhb = ylb + boxSize

#define list of timeSteps
AIlist = list(range(0,12,2)) + list(range(20,210,10)) + list(range(400,1200,200)) + [1400,2000,2400,3000,3400,4000]

#open file to write density information into
fileName = 'densityInfo.txt'
fill = open(fileName,'w')
fill.write('timeStep\tdensity(kg/m^3)\n')

for ai in AIlist:
    
    #read file and unpack the values into single Vdfull array
    fileName = "fullT"+str(ai)+"SN"+str(SNo)+".dat" 
    
    bnf = open(fileName,"rb")
    fcont = bnf.read()
    bnf.close()
    
    fullsize = int(len(fcont)/8)
    Vdfull = struct.unpack('d'*fullsize,fcont[:fullsize*8]) #contains float information for the data file
    
    size = xsize*ysize*zsize
    
    Xd = []
    Yd = []
    Zd = []
    Vd = []
    
    #unpack Vdfull buffer into rho array and create array for axes
    counti = 0
    totalP = 0
    for z in range(zsize):
        for y in range(ysize):
            for x in range(xsize):
                countv = counti*(noparts+1)            
                
                if Vdfull[countv] > lim:
                    Zd.append(x)
                    Yd.append(y)
                    Xd.append(z)
                    Vd.append(Vdfull[countv])
                    totalP += 1
                    
                counti += 1
    
    Xd = np.array(Xd)
    Yd = np.array(Yd)
    Zd = np.array(Zd)
    Vd = np.array(Vd)
    
    #find values in the analysis box
    idd = (Xd >= xlb) & (Xd <= xhb) & (Yd >= ylb) & (Yd <= yhb)
    xbox = Xd[idd]; ybox = Yd[idd]; zbox = Zd[idd]; vbox = Vd[idd]
    
    str2write = str(ai) + '\t'
    if ai == 0:
        init_volume = (max(Xd) - min(Xd))*(max(Yd) - min(Yd))*(max(Zd) - min(Zd))
        init_rho_ratio = sum(Vd>0.1)/init_volume
        sumrho_init = sum(vbox)
        init_density = bulkR*init_rho_ratio
        str2write += str(init_density) + '\n'
    else:
        sumrho = sum(vbox)
        rel_rho_change = (sumrho - sumrho_init)/sumrho
        density = init_density/(1-rel_rho_change)
        str2write += str(density) + '\n'
    fill.write('%s' %str2write)
    print(ai,'-- done')

fill.close()