import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math as m
import struct
import os

def resize(List,ch=1): 
    if ch:
        newList = []
        for ix,iy,iz in zip(List[0],List[1],List[2]):
            newList.append([ix,iy,iz])
        newList.sort(key = lambda x:(x[0],x[1],x[2])) #(x[2],x[1],x[0])
    else:
        newList = [[],[],[]]
        for ii in List:
            newList[0].append(ii[0])
            newList[1].append(ii[1])
            newList[2].append(ii[2])
    return newList

def makeParticlesAppend(pinA):
    
    sphereList = [[],[],[]]
    for xC,yC,zC,R in zip(pinA[0],pinA[1],pinA[2],pinA[3]):
        
        lowb = [m.floor(xC-R),m.floor(yC-R),m.floor(zC-R)]
        diam = m.ceil(2*R)
        
        #print ('Circle properties-->',xC,yC,R)
        for ix in range(lowb[0],lowb[0]+diam+3):
            for iy in range(lowb[1],lowb[1]+diam+3):
                for iz in range(lowb[2],lowb[2]+diam+3):
                    if (ix-xC)**2+(iy-yC)**2+(iz-zC)**2<=R**2:
                        sphereList[0].append(ix)
                        sphereList[1].append(iy)
                        sphereList[2].append(iz)
    
    sphereList = resize(resize(sphereList),0)
    
    return sphereList


def makeParticles(pinA,xsize,ysize,zsize,ch=0):
    
    sphereList = [[],[],[]]
    size = xsize*ysize*zsize
    sphereList[0] = [1000 for ii in range(size)]
    sphereList[1] = [1000 for ii in range(size)]
    sphereList[2] = [1000 for ii in range(size)]
    
    for xC,yC,zC,R in zip(pinA[0],pinA[1],pinA[2],pinA[3]):
        
        if ch:
            xC = round(xC); yC = round(yC); zC = round(zC); R = round(R)
        
        lowb = [m.floor(xC-R),m.floor(yC-R),m.floor(zC-R)]
        diam = m.ceil(2*R)
        highb = [lowb[0]+diam+3,lowb[1]+diam+3,lowb[2]+diam+3]
        
        #print ('Circle properties-->',xC,yC,R)
        for ix in range(lowb[0],highb[0]):
            for iy in range(lowb[1],highb[1]):
                for iz in range(lowb[2],highb[2]):
                    if ((ix-xC)**2+(iy-yC)**2+(iz-zC)**2)<=R**2:
                        i = xsize*ysize*iz + xsize*iy + ix
                        sphereList[0][i] = ix
                        sphereList[1][i] = iy
                        sphereList[2][i] = iz
    
    for ct in range(3):
        sphereList[ct] = np.setdiff1d(sphereList[ct],1000,assume_unique = True)
    
    sphereList = resize(resize(sphereList),0)
    
    return sphereList

xsize = 94
ysize = 104
zsize = 104
zBox = 110 

lim = 0.1

elev = 30
azim = 45

colbd = 'black'
colfc = 'white'

convfac = 0.010584 #pixel to um conversion factor

curdr = str(os.getcwd())

dirList = ['generation1','generation2'] #['generation2']#['sintering1','sintering2']#,

for dirName in dirList:
    
    os.chdir(dirName)
    
    if dirName[:-1] == 'sintering':
        if dirName[-1] == '0':
            SNo = 50400002
            xsize = 61
            fStart = 0
            fEnd = 1
        elif dirName[-1] == '1':
            SNo = 50400002
            xsize = 61
            fStart = 0
            fEnd = 102            
        else:
            SNo = 504
            xsize = 94
            fStart = 0
            fEnd = 102
            
        for ai in range(fStart,fEnd,2):
        
            bnf = open("rho%dSN%s.dat" %(ai,SNo),"rb")
            fcont = bnf.read()
            size = int(len(fcont)/8)
            Vd = struct.unpack('d'*size,fcont[:size*8])
            bnf.close()   
            
            x = []; y = []; z = []; v = [];
        
            counti = 0
            for c in range(zsize):
                for b in range(ysize):
                    for a in range(xsize):
                        if Vd[counti] > lim:
                            z.append(a)
                            y.append(b)
                            x.append(c)
                            v.append(Vd[counti])
                        counti = counti+1    
            
            v = np.array(v)
            fig = plt.figure()
            ax = Axes3D(fig)
        
            ax.set_xlim([0,ysize])
            ax.set_ylim([0,ysize])
            ax.set_zlim([0,zBox])
        
            ax.scatter(x,y,z,s=100*v,c= colfc,marker='o',edgecolor = colbd)
            ax.view_init(elev,azim)
            
            fig.savefig("plotrho%dSN%s.png" %(ai,SNo))
            
            plt.close()
            
            print(dirName,'file',ai,'Done')
    
    else:
        if dirName=='generation2':
   
            inRead = 'ASurfaceNorms.txt' 
            
            surf = [[],[],[]]
             
            with open(inRead,'r') as finl:
                for ii in finl.readlines():
                    cols = ii.split('\t')
                    surf[0].append(int(cols[0]))
                    surf[1].append(int(cols[1]))
                    surf[2].append(int(cols[2]))  
            
            xsize = 94
            fStart = 0
            fstep = 1
            fEnd = 160
            ch = 0
        else:
            xsize = 61
            fStart = 0
            fstep = 1
            fEnd = 67
            ch = 1
            
        for ai in range(fStart,fEnd,fstep):
            pinA = [[],[],[],[]]
            finalName = 'particleOut'+ str(ai)+'.txt'
            with open(finalName,'r') as finl:
                for ii in finl.readlines():
                    cols = ii.split('\t')
                    pinA[0].append(float(cols[0])/convfac)
                    pinA[1].append(float(cols[1])/convfac)
                    pinA[2].append(float(cols[2])/convfac)
                    pinA[3].append(float(cols[3])/convfac)
                    
            maxZ = m.ceil(max([z+r for z,r in zip(pinA[2],pinA[3])])+1)
            
            sphereList = makeParticlesAppend(pinA,zsize,ysize,maxZ,ch)
            
            fig = plt.figure()
            ax = Axes3D(fig)
                            
            ax.set_xlim([0,ysize])
            ax.set_ylim([0,ysize])
            ax.set_zlim([0,zBox])
            
            ax.scatter(sphereList[0],sphereList[1],sphereList[2],s = 9.8,c= colfc,marker='o',edgecolor = colbd)
            
            if dirName=='generation2':
                ax.scatter(surf[0],surf[1],surf[2],s = 9.8,c= colfc,marker='o',edgecolor = colbd)
                
            ax.view_init(elev,azim)
                      
            fig.savefig("Append2plotParts%d.png" %ai)
            plt.close() 
            
            print(dirName,'file',ai,'Done')
    
    os.chdir(curdr)