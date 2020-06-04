import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.animation as animation
from pylab import *
import struct
                        
SNo = 50400002

colbd = 'black'
colfc = 'red'

old=True
plotEta = False
plotRho = False

for ai in [10]:
    for ang in range(0,1):

        if ang%2==0:
            elev = 40
        else:
            elev = 70

        if old:
            bnf = open("fullT%dSN%d.dat" %(ai,SNo),"rb")
            xsize = 61
            noparts = 21
        else:
            bnf = open("newLayer_bottom.dat","rb")#open("newLayer.dat","rb")
            xsize = 94
            noparts = 58
            
        fcont = bnf.read()
        fullsize = len(fcont)/8
        Vdfull = struct.unpack('d'*fullsize,fcont[:fullsize*8]) #contains float information for the data file
        
        print "length of an array is: %d" % fullsize

        ysize = 104
        zsize = 104

        pCut = 1
        etaCut = 0.1
        
        size = xsize*ysize*zsize

        Xd = [0 for x in range(size)]
        Yd = [0 for x in range(size)]
        Zd = [0 for x in range(size)]
        Vd = [0 for x in range(size)]
        Ed = [[0 for y in range(noparts)] for x in range(size)]
        
        ad = np.zeros(size)
        bd = np.zeros(size)
        cd = np.zeros(size)
        vd = np.zeros(size)

        ysizep = ysize;

        counti = 0
        for z in range(zsize):
            for y in range(ysize):
                for x in range(xsize):
                    Zd[counti] = x
                    Yd[counti] = y
                    Xd[counti] = z
                    
                    countv = counti*(noparts+1)
                    Vd[counti] = Vdfull[countv]

                    for p in range(noparts):
                        Ed[counti][p] = Vdfull[countv+p+1]
                        
                    counti = counti+1
            
        lim = 0.1
        if plotRho:
            colfc = 'green'
            fig2 = plt.figure()
            ax2 = Axes3D(fig2)

            ax2.set_xlim([0,ysizep])
            ax2.set_ylim([0,ysizep])
            ax2.set_zlim([0,ysizep])
                
            ax2.view_init(elev,azim=50)

            countl = 0
            for i in range(0,size):
                if Vd[i] > lim:
                    ad[i] = Xd[i]
                    bd[i] = Yd[i]
                    cd[i] = Zd[i]
                    vd[i] = Vd[i]
                    countl = countl+1
                else:
                    ad[i] = 1000.0
                    bd[i] = 1000.0
                    cd[i] = 1000.0
                    vd[i] = 1000.0

            countth = np.count_nonzero(vd==1000.0)

            val = np.array([1000.0])

            x = np.setdiff1d(ad,val,assume_unique = True)
            y = np.setdiff1d(bd,val,assume_unique = True)
            z = np.setdiff1d(cd,val,assume_unique = True)
            v = np.setdiff1d(vd,val,assume_unique = True)

            ax2.scatter(x,y,z,s=100*v,c= colfc,marker='o',edgecolor = colbd)
            
            fig2.savefig("newOldjustrho%dSN%d_All.png" %(ai, SNo))
            
        if plotEta:
            fig2 = plt.figure()
            ax2 = Axes3D(fig2)

            ax2.set_xlim([0,ysizep])
            ax2.set_ylim([0,ysizep])
            ax2.set_zlim([0,ysizep])
                
            ax2.view_init(elev,azim=50)
            
            for pCut in range(noparts):
                countl = 0
                for i in range(0,size):
                    if Vd[i] > lim and Ed[i][pCut] >= etaCut:
                        ad[i] = Xd[i]
                        bd[i] = Yd[i]
                        cd[i] = Zd[i]
                        vd[i] = Vd[i]
                        countl = countl+1
                    else:
                        ad[i] = 1000.0
                        bd[i] = 1000.0
                        cd[i] = 1000.0
                        vd[i] = 1000.0

                countth = np.count_nonzero(vd==1000.0)

                val = np.array([1000.0])

                x = np.setdiff1d(ad,val,assume_unique = True)
                y = np.setdiff1d(bd,val,assume_unique = True)
                z = np.setdiff1d(cd,val,assume_unique = True)
                v = np.setdiff1d(vd,val,assume_unique = True)

                fig = plt.figure()
                ax = Axes3D(fig)

                ax.set_xlim([0,ysizep])
                ax.set_ylim([0,ysizep])
                ax.set_zlim([0,ysizep])

                ax.scatter(x,y,z,s=100*v,c= colfc,marker='o',edgecolor = colbd)
                ax.view_init(elev,azim=50)

                ax2.scatter(x,y,z,s=100*v,c= colfc,marker='o',edgecolor = colbd)
                
                #plt.show()
                print pCut
                fig.savefig("newOldrho%dSN%d_fromFull_eta%d_mt%f.png" %(ai, SNo,pCut, etaCut))
                close()

            fig2.savefig("newOldrho%dSN%d_All.png" %(ai, SNo))
        bnf.close()
        print "%d_%d done" %(ai, elev)
