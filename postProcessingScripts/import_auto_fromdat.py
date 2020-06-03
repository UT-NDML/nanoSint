import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.animation as animation
from pylab import *
import struct

SNoAll = ['504']

colbd = 'black'
colfc = 'white'

bdn = 5

ctt = 40
    
for SNo in SNoAll:
    
    xsize = 94
    ysize = 104
    zsize = 104
    lowb = 25
    highb = 75
    
    for ai in range(0,5050,50):
        for ang in range(0,1):

            if ang%2==0:
                elev = 40
            else:
                elev = 70

            #bnf = open("rho%dSN%d_bed%dp6.dat" %(ai,SNo,bdn),"rb")
            bnf = open("rho%dSN%s.dat" %(ai,SNo),"rb")
            fcont = bnf.read()
            size = len(fcont)/8
            Vd = struct.unpack('d'*size,fcont[:size*8])
            bnf.close()
            print "length of an array is: %d" % size

            Xd = [0 for x in range(size)]
            Yd = [0 for x in range(size)]
            Zd = [0 for x in range(size)]
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
                        counti = counti+1

            lim = 0.1
            countl = 0

            temp = []
            for i in range(0,size):
                if Vd[i] > lim and lowb<=Yd[i]<=highb and lowb<=Xd[i]<=highb:#Yd[i] < ysize/2:
                    x = Xd[i]
                    y = Yd[i]
                    z = Zd[i]
                    rad = sqrt((x-50)**2+(y-50)**2+(z-100)**2)
                    if rad<=30: temp.append('r')
                    elif rad<=40: temp.append('y')
                    else: temp.append(colbd)
                    ad[i] = x
                    bd[i] = y
                    cd[i] = z
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

            ax.set_xlim([0,ysize])
            ax.set_ylim([0,ysize])
            ax.set_zlim([0,ysize])

            ax.scatter(x,y,z,s=100*v,c= colfc,marker='o',edgecolor = temp)
            ax.view_init(elev,azim=50)
            #ax.view_init(elev=10,azim=80)

            #plt.show()
            print ai
            print elev
            fig.savefig("Boxrho%dSN%s.png" %(ai,SNo))
            close()
            print "%d_%d done" %(ai, elev)

    ctt = ctt + 1
    
    print "%s===>Done" %SNo
