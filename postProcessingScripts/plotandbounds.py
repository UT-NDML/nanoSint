import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.animation as animation
from pylab import *
import shutil 
import struct
import os

SNo = 50400002
ai = 40800

startbeds = 0
nobeds = 18

colbd = 'purple'
colfc = 'white'

height = []

dirstd = str(os.getcwd())
'''
with open('TACCguide.log','r') as lin:
    for dks in lin:
        height.append(10*int(dks[3])+int(dks[4])+2)

boxS = open('BoxSize.txt','w')

for bedno in range(startbeds,nobeds):

    boxSize = 1
    
    newdr = dirstd+'/bed'+str(bedno)
    os.chdir(newdr)

    print str(os.getcwd())
    
    if os.path.exists("rho%dSN%d.dat" %(ai,SNo)):
        
        xsize = height[bedno]
        ysize = 104
        zsize = 104
        
        for ai in range(0,ai+50,ai):
            for ang in range(0,1):

                if ang%2==0:
                    elev = 40
                else:
                    elev = 70

                bnf = open("rho%dSN%d.dat" %(ai,SNo),"rb")
                fcont = bnf.read()
                size = len(fcont)/8
                Vd = struct.unpack('d'*size,fcont[:size*8])

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

                fig = plt.figure()
                ax = Axes3D(fig)

                ysizep = ysize

                ax.set_xlim([0,ysizep])
                ax.set_ylim([0,ysizep])
                ax.set_zlim([0,ysizep])

                ax.scatter(x,y,z,s=100*v,c= colfc,marker='o',edgecolor = colbd)
                ax.view_init(elev=40,azim=50)

                #plt.show()
                print ai
                print elev
                fig.savefig("rho%dSN%d_elev%d.png" %(ai, SNo, elev))
                close()
                print "%d_%d done" %(ai, elev)

                bnf.close()
                

        bnf = open("rho%dSN%d.dat" %(ai,SNo),"rb")
        fcont = bnf.read()
        size = len(fcont)/8
        Vd = struct.unpack('d'*size,fcont[:size*8])

        Xd = [0 for x in range(size)]
        Yd = [0 for x in range(size)]
        Zd = [0 for x in range(size)]

        boxSize = 40

        errfac = 0.3
        #errfacz = 0.6

        counti = 0
        for z in range(zsize):
                for y in range(ysize):
                        for x in range(xsize):
                                Zd[counti] = x
                                Yd[counti] = y
                                Xd[counti] = z
                                counti = counti+1

        lim = 0.1

        ZstripInfo = [[0 for y in range(6)] for x in range(ysize*zsize)]

        #devnfile = open('BottomDeviation2.txt', 'w')
        #devnfile.write("zvalue(x)\tyvalue\tminz\tnorm_min\n")
                      
        ctr = 0;

        for zi in range(0,zsize):
                for b in range(ysize):
                        countz = 0
                        for a in range(xsize):
                                c = zi
                                k = a + xsize*b + xsize*ysize*c
                                if Vd[k] > lim:
                                        countz = countz+1
                                else:
                                        countz = countz+0

                        if (countz!=0):
                                zz = np.zeros(countz)
                                totalz = 0
                                for a in range(xsize):
                                        c = zi
                                        k = a + xsize*b + xsize*ysize*c
                                        if Vd[k] > lim:
                                                zz[totalz] = Zd[k]
                                                totalz = totalz + 1
                
                                minzz = min(zz)
                                maxzz = max(zz)
                        else:
                                minzz = xsize
                                maxzz = 0
                                
                        ZstripInfo[ctr][0] = zi
                        ZstripInfo[ctr][1] = b
                        ZstripInfo[ctr][2] = minzz
                        ZstripInfo[ctr][3] = maxzz

                        ctr = ctr+1;
                        
                        #devnfile.write("%d\t%d\t%d\t%d\n" %(zi,b,minzz,maxzz))

        minallZ = min([x[2] for x in ZstripInfo]) #find the minimum non zero z
        maxallZ = max([x[3] for x in ZstripInfo]) #find the maximum z

        for sp in range(zsize*ysize):
                        
                ZstripInfo[sp][4] = (ZstripInfo[sp][2] - minallZ)/(xsize-minallZ)
                ZstripInfo[sp][5] = (maxallZ - ZstripInfo[sp][3])/maxallZ

                #devnfile.write("%d\t%d\t%d\t%f\n" %(ZstripInfo[sp][0],ZstripInfo[sp][1],ZstripInfo[sp][2],ZstripInfo[sp][4]))
                
        #devnfile.close()

        pixlimit = 5
        boxlimit = 57
                
        outlog = open("logFile_bed%d.txt" %bedno,'w')
        for pixdev in np.arange(2,pixlimit+1):

                errfac = pixdev/(xsize-minallZ)

                outfile = open("ValidStrips_deviation_%dpixels.txt" %pixdev,'w')
                
                #using only the deviation from the bottom to determine acceptable beds
                #errfac above tells how much leeway we're willing to give in this deviation from the base

                rdsize = 0
                rdsizexy = 0
                
                for sp in range(zsize*ysize):
                        if (ZstripInfo[sp][4] <= errfac):
                                rdsizexy = rdsizexy+1
                                for zc in range(int(ZstripInfo[sp][2]),int(ZstripInfo[sp][3])+1):
                                        rdsize = rdsize+1

                xrd = np.zeros(rdsize)
                yrd = np.zeros(rdsize)
                zrd = np.zeros(rdsize)
                vrd = np.zeros(rdsize)

                xrdxy = np.zeros(rdsizexy)
                yrdxy = np.zeros(rdsizexy)
                
                totc = 0
                totcxy = 0
                
                for spr in range(zsize*ysize):
                        if (ZstripInfo[spr][4] <= errfac):
                                
                                xrdxy[totcxy] = ZstripInfo[spr][0]
                                yrdxy[totcxy] = ZstripInfo[spr][1]

                                totcxy = totcxy + 1 

                                for zc in range(int(ZstripInfo[spr][2]),int(ZstripInfo[spr][3])+1):
                                        k = zc + xsize*ZstripInfo[spr][1] + xsize*ysize*ZstripInfo[spr][0]
                                        xrd[totc] = ZstripInfo[spr][0]
                                        yrd[totc] = ZstripInfo[spr][1]
                                        zrd[totc] = zc
                                        vrd[totc] = Vd[k]

                                        totc = totc + 1
                                        
                                #print "%d\t%d" %(ZstripInfo[spr][0],ZstripInfo[spr][1])          
                                outfile.write("%d\t%d\n" %(xrdxy[totcxy-1],yrdxy[totcxy-1]))


                outfile.close()
             
                #fig = plt.figure(frameon=False)
                fig = plt.figure()
                ax = Axes3D(fig)
                #ax.axis('off')

                ax.set_xlim3d([0,100])
                ax.set_ylim3d([0,100])
                ax.set_zlim3d([0,100])
                plt.xlabel('x')
                plt.ylabel('y')
                ax.scatter(xrd,yrd,zrd,s=100*vrd,c='m',marker='o',edgecolor='black')
                ax.view_init(elev=40,azim=50)

                fig.savefig("Bottom_rho%dSN%d_pixelDeviation_%d.png" %(ai,SNo,pixdev))

                print "Plot deviation: %d__errorfac: %f done,start boxfind" %(pixdev,errfac)

                #Finding appropriate boxes

                outlog.write("Pixel deviation: %d\n" %pixdev)

                BSize = boxSize*boxSize

                xboxG = []
                yboxG = []

                if pixdev==pixlimit:
                    boxSize = 45
                else:
                    boxSize = boxlimit

                while boxSize < boxlimit+1:
                        Xall = []
                        Yall = []
                        gdb = 0
                        
                        for xpos in range(int(min(xrdxy)),int(max(xrdxy))-boxSize+1):
                                for ypos in range(int(min(yrdxy)),int(max(yrdxy))-boxSize+1):
                                        ctbox = 0

                                        
                                        ##first start by testing for edges of the box, proceed inwards if edges are present
                                        for ig in range(4):
                                                xbc = xpos+(boxSize-1)*(ig/2)
                                                ybc = ypos+(boxSize-1)*mod(ig,2)
                     
                                                xind=[i for i, x in enumerate(xrdxy) if x == xbc] 
                                                yind=[i for i, y in enumerate(yrdxy) if y == ybc]
                                                tmparr=set(xind).intersection(yind)

                                                if len(tmparr)==1:
                                                        ctbox=ctbox+1
                                                        #print "good %d" %ctbox
                                                else:
                                                        break
                
                                        if (ctbox==4):
                                                Xall.append(xpos)
                                                Yall.append(ypos)
                                                gdb = gdb+1

                        if gdb!=0:
                                outfilebx = open("BiggestBeds_ValidBoxStarts_deviation_%dpixs_%dSize.txt" %(pixdev,boxSize),'w')
                                for ii in range(gdb):
                                        outfilebx.write("%d\t%d\n" %(Xall[ii],Yall[ii]))
                                outfilebx.close()
                                
                                
                        if gdb==0:
                                break

                        outlog.write("%d complete boxes found in good area_Size %d\n" %(gdb,boxSize))                                  
                        boxSize = boxSize+1
                    
                print "BoxFind deviation: %d__errorfac: %f done" %(pixdev,errfac)

                if pixdev==pixlimit or boxSize == boxlimit+1:
                    shutil.copyfile("BiggestBeds_ValidBoxStarts_deviation_%dpixs_%dSize.txt" %(pixdev,boxSize-1),"BoundstoUse.txt")
                    break

        outlog.close()
                
    os.chdir(dirstd)
    print "BoxSize: %d" %(boxSize-1)
    boxS.write("%d\n" %(boxSize-1))
    
boxS.close()
'''
os.chdir(dirstd)

dirtomake = dirstd+'/BoundsInfo'

os.mkdir(dirtomake)

shutil.copy2("BoxSize.txt",dirtomake)

for bedno in range(startbeds,nobeds):

    newdr = dirstd+'/bed'+str(bedno)
    mdir = dirtomake+'/bed'+str(bedno)

    os.mkdir(mdir)
    os.chdir(newdr)

    print str(os.getcwd())
    
    if os.path.exists("rho%dSN%d.dat" %(ai,SNo)):
        shutil.copy2("BoundstoUse.txt",mdir)

    os.chdir(dirstd)
