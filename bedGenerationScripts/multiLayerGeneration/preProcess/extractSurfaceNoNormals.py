import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.animation as animation
from pylab import *
import struct
import random
import time

#Print ints in array to file. Prints list with [] if there is a list in the array
def toFile(Pd,fname='outPosition.txt',filout='int'):

    outf = open(fname,'w')
    
    for ij in Pd:
        for kk in ij:
            if type(kk)==list:
                outf.write('[')
                for ji in kk:
                    if type(ji)==int:
                        outf.write('%d,' %ji)
                    else:
                        if filout=='float':
                           outf.write('%f,' %ji) 
                outf.write(']')
            else:
                if type(kk)==int:
                    outf.write('%d' %kk)
            outf.write('\t')
        outf.write('\n')

    outf.close()

start_t = time.time()

SNoAll = [50400002]

aiAll = [100] #range(0,500,50)+range(500,1000,100)+[1000,5000]+range(10000,40000,10000)+[40800]

colbd = 'green'
colfc = 'white'

lim = 0.1
limH = 0.25
                
ctt = 1

SurfaceExtract = 0 #0 if reading from file, 1 if calculating normal
    
for SNo in SNoAll:
    
    xsize = 61
    ysize = 104
    zsize = 104
            
    for ai in aiAll:

        if SurfaceExtract==1:
            bnf = open("rho%dSN%d.dat" %(ai,SNo),"rb")
            fcont = bnf.read()
            size = len(fcont)/8
            Vd = struct.unpack('d'*size,fcont[:size*8])

            print "length of an array is: %d" % size

            #surfout = open('PointsSurface2.txt', 'w')

            Xd = [0 for x in range(size)]
            Yd = [0 for x in range(size)]
            Zd = [0 for x in range(size)]
            Px = [0 for x in range(size)]
            ad = np.zeros(size)
            bd = np.zeros(size)
            cd = np.zeros(size)
            vd = np.zeros(size)

            counti = 0
            for z in range(zsize):
                for y in range(ysize):
                    for x in range(xsize):
                        Zd[counti] = x
                        Yd[counti] = y
                        Xd[counti] = z
                        Px[counti] = counti
                        counti = counti+1

            x = []
            y = []
            z = []
            v = []
            p = []

            for ii in range(size):
                if Vd[ii]>lim:
                    x.append(Xd[ii])
                    y.append(Yd[ii])
                    z.append(Zd[ii])
                    v.append(Vd[ii])
                    p.append(Px[ii])

            ssize = len(x)
            pos = []

            pSet = set(p)
            
            print ssize
            print "Isolating surface"
            
            for xx in range(ssize):
                pl = p[xx]-1
                pr = p[xx]+1
                pt = x[xx]*xsize*ysize + (y[xx]+1)*xsize + z[xx]
                pd = x[xx]*xsize*ysize + (y[xx]-1)*xsize + z[xx]
                pzh = (x[xx]+1)*xsize*ysize + y[xx]*xsize + z[xx]
                pzl = (x[xx]-1)*xsize*ysize + y[xx]*xsize + z[xx]

                if Vd[pl]>lim and Vd[pr]>lim and Vd[pt]>lim and Vd[pd]>lim and Vd[pzh]>lim and Vd[pzl]>lim:
                    pass
                else:
                    #print x[xx],y[xx],z[xx],v[xx]
                    #surfout.write('%d\t%d\t%d\n' %(x[xx],y[xx],z[xx]))
                    pos.append([x[xx],y[xx],z[xx],v[xx]])
                    
                '''
                if pl in pSet and pr in pSet and pt in pSet and pd in pSet and pzh in pSet and pzl in pSet:
                    pass
                else:
                    print x[xx],y[xx],z[xx],v[xx]
                    surfout.write('%d\t%d\t%d\n' %(x[xx],y[xx],z[xx]))
                    pos.append([x[xx],y[xx],z[xx],v[xx]])
                '''

            #surfout.close()

            toFile(pos,'PointSurfaceFunc.txt')

            print "Surface found\nElapsed Time", (time.time()-start_t)/60, "mins\nWorking on sparser cloud"

            #closout = open('ClosestIndicesMin.txt', 'w')
            
            noc = 5 #number of sparse points to leave
            spax = 5 #spacing to allow in grid search
            vac = 2000 #amount of points to consider for search

            
            lenpos = len(pos)
            ctpp = 0
            ctr = 0
            posSend = []
            
            for pp in pos:
                #print pp
                
                if ctpp==0:
                    posSend.append(pp[:3])                
                elif ctpp == noc:
                    ctpp = 0
                    posSend.append(pp[:3])
                elif pp == pos[-1]:
                    posSend.append(pp[:3])
                
                ctpp+=1

                
        
            print "Done with closest points\nElapsed Time", (time.time()-start_t)/60, "mins"

            print "Saving Information to File"
            
            toFile(posSend,'ASurfacePoints%d.txt' %noc)

            print "Done!\nElapsed Time", (time.time()-start_t)/60, "mins"

            print "All ==> ", len(pos),";", noc, " ==> ", len(posSend)
            
        else:
            pos = []

            noc = 1
            
            with open('ASurfacePoints%d.txt' %noc,'r') as rd:
                for lin in rd:
                    temp = []
                    lin = lin.replace('[','\t')
                    lin = lin.replace(']','\t')
                    lin = lin.replace(',','\t')
                    lin = lin.split('\t')
                    for ik in lin:
                        try:
                            temp.append(int(ik))
                        except ValueError:
                            try:
                                temp.append(float(ik))
                            except ValueError: pass
                    pos.append(temp)
       
            print "Done reading in file"
            print "Plotting Surface points\n"
            
            shdplot = 1
                    
            if shdplot==1:
                gdd = 0
                #grid = [[0,20],[78,98],[35,55]]
                #grid = [[0,30],[25,55],[30,60]]
                #grid = [[30,60],[40,70],[30,60]]
                grid = [[0,zsize],[0,ysize],[0,xsize]]      
                
                print "Plotting Surface"
                temp = [[],[],[],[]]
                
                for pp in pos:
                    if pp[0]>=grid[0][0] and pp[0] <= grid[0][1] and pp[1]>=grid[1][0] and pp[1] <= grid[1][1] and pp[2]>=grid[2][0] and pp[2] <= grid[2][1]:
                        temp[0].append(pp[0])
                        temp[1].append(pp[1])
                        temp[2].append(pp[2])
                        temp[3].append(9.8)
                        
                fig = plt.figure()
                ax = Axes3D(fig)

                ax.set_xlim([0,ysize])
                ax.set_ylim([0,ysize])
                ax.set_zlim([0,ysize])
                
                ax.scatter(temp[0],temp[1],temp[2],s=temp[3],c= colfc,marker='o',edgecolor = colbd)
                #ax.view_init(elev=0,azim=0)
                ax.view_init(elev=40,azim=50)
                #ax.view_init(elev=90,azim=50)

                #print "Plotting Normals"
                
                #grid = [[15,20],[78,85],[35,55]]
                #grid = [[0,20],[78,98],[35,55]]

                #print "saving"
                fig.savefig("A2CheckSurfaceSparse%d.png" %noc)
                close()
                #print gdd,"=> x: ",grid[0],"y: ",grid[1],"z: ",grid[2]
                gdd += 1

    ctt = ctt + 1
    
    print "%d===>Done" %SNo
