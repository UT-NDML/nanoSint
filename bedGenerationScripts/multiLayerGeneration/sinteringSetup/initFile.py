import struct
from array import array

def CreateParticles(lis,xsize,ysize,zsize,circ):
    
    for z in range(zsize):
        for y in range(ysize):
            for x in range(xsize):
                
                if ((x-circ.midX)**2+(y-circ.midY)**2+(z-circ.midZ)**2)<=circ.radius**2:
                    i = xsize*ysize*z + xsize*y + x
                    lis[i]=circ.value                    
    

class circleProps:
    def __init__(self):
        self.midX = 0
        self.midY = 0
        self.midZ = 0
        self.radius = 0
        self.value = 0
    
    

rhoSolid = 0.9998
rhoVap = 0.000000089
etaSolid = 1.0
    
SNo = 50400002
ai = 10

bnf = open("fullT%dSN%d.dat" %(ai,SNo),"rb")
fcont = bnf.read()
fullsize = int(len(fcont)/8)
Vdfull = struct.unpack('d'*fullsize,fcont[:fullsize*8]) #contains float information for the data file
bnf.close()

print("length of an array is: %d" % fullsize)

#Old bed extaction
xsize = 61
ysize = 104
zsize = 104

size = xsize*ysize*zsize

pCut = 1
etaCut = 0.5

noparts = 21

xsizeNew = 94
nopartsNew = 58

sizeNew = xsizeNew*ysize*zsize

eta = [[0 for y in range(sizeNew)] for x in range(nopartsNew+1)]

trackx = 0
ctx = 0

for counti in range(size):
    countv = counti*(noparts+1)
    
    countl=counti + ctx*(xsizeNew-xsize)
    
    for p in range(noparts+1):
        eta[p][countl] = Vdfull[countv+p]
            
    if trackx < xsize-1:
        trackx += 1
    elif trackx == xsize-1:
        trackx = 0
        ctx += 1

'''
## testing only bottom layer
## comment out
        
#transfer data into single array
Vdfull2 = [0 for x in range(sizeNew*(nopartsNew+1))]
fill = open('newLayer_bottom.dat','wb')
ct = 0

for p in range(sizeNew):
    if eta[0][p] == 0:
        eta[0][p]=rhoVap
    for g in range(nopartsNew+1):
        Vdfull2[ct]=eta[g][p]
        ct += 1

s = array('d',Vdfull2)
s.tofile(fill)
fill.close()
'''

#new Layer extraction
circleProp = []
fileName='zmax0.8_res0.99_noB9_bed6_bed.txt'
with open(fileName,'r') as finl:
    for ii in finl.readlines():
        [r,z,y,x] = ii.split()  
        circleProp.append([float(r),float(x),float(y),float(z)])
        
        
#Make new particles
circ = circleProps()
ctr = 0
for g in range(noparts+1,nopartsNew+1):
    circ.radius = circleProp[ctr][0]
    circ.midX = circleProp[ctr][1]
    circ.midY = circleProp[ctr][2]
    circ.midZ = circleProp[ctr][3]
    circ.value = rhoSolid
    CreateParticles(eta[0], xsizeNew, ysize, zsize, circ)
    circ.value = etaSolid
    CreateParticles(eta[g], xsizeNew, ysize, zsize, circ)    
    ctr += 1
    
    
#transfer data into single array
Vdfull2 = [0 for x in range(sizeNew*(nopartsNew+1))]
fill = open('newLayer.dat','wb')
ct = 0

for p in range(sizeNew):
    if eta[0][p] == 0:
        eta[0][p]=rhoVap
    for g in range(nopartsNew+1):
        Vdfull2[ct]=eta[g][p]
        ct += 1

s = array('d',Vdfull2)
s.tofile(fill)
fill.close()
