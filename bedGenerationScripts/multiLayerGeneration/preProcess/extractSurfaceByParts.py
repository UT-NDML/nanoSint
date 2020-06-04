import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import matplotlib.animation as animation
from pylab import *
import struct
import random

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

#Calculates the normal vector given three points <only x,y and z>
def calcNormal(midP,p1,p2,dbg=0):
    vec1 = [0 for ii in range(3)]
    vec2 = vec1[:]

    for ii in range(3):
        vec1[ii] = p1[ii] - midP[ii]
        vec2[ii] = p2[ii] - midP[ii]

    normF = np.cross(vec1,vec2)

    if dbg==1:
        print vec1
        print vec2
        print normF
    
    return [normF[0].item(),normF[1].item(),normF[2].item()]

def findNormal(pp,pos):
    
    fdNorm = calcNormal(pp[0:3],pos[pp[6][0]][0:3],pos[pp[6][1]][0:3])

    id2 = 1
    
    while sum(np.abs(fdNorm))==0:
        #print pp[0:3]
        #print pos[pp[6][0]][0:3]
        #print pos[pp[6][id2]][0:3]
        id2 += 1
        #print "id ==> %d, norm = 0, try %d" %(pp[4],id2)
        fdNorm = calcNormal(pp[0:3],pos[pp[6][0]][0:3],pos[pp[6][id2]][0:3])

    return fdNorm

def signChange(arr):
    ct = 0
    for ii in arr:
        if ct == 0:
            sn = np.sign(ii)           
        else:
            if np.sign(ii)==sn or ii==0:
                pass
            else:
                return True
        ct +=1
        
    return False

SNoAll = [50400002]

aiAll = [100] #range(0,500,50)+range(500,1000,100)+[1000,5000]+range(10000,40000,10000)+[40800]

colbd = 'green'
colfc = 'white'

lim = 0.1
limH = 0.25
                
ctt = 1

NormalCalc = 0 #0 if reading from file, 1 if calculating normal
    
for SNo in SNoAll:
    
    xsize = 61
    ysize = 104
    zsize = 104
            
    for ai in aiAll:

        if NormalCalc==1:
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
            
            calcNormal(pos[0][0:3],pos[1][0:3],pos[3][0:3])

            print "Surface found\nWorking on closest points"

            #closout = open('ClosestIndicesMin.txt', 'w')
            
            noc = 10 #number of closest points to get
            spax = 5 #spacing to allow in grid search
            vac = 2000 #amount of points to consider for search

            
            lenpos = len(pos)
            ctpp = 0
            
            for pp in pos:
                #print pp
                closest = []
                inspax = []
                for ip in range(max(ctpp-vac,0),min(ctpp+vac,lenpos)):
                    if abs(pos[ip][0]-pp[0])<=spax and abs(pos[ip][1]-pp[1])<=spax and abs(pos[ip][2]-pp[2])<=spax:
                        dist = (pos[ip][0]-pp[0])**2+(pos[ip][1]-pp[1])**2+(pos[ip][2]-pp[2])**2
                        if dist==0:
                            pp.append(ip)
                        else:
                            inspax.append([ip,dist])

                inspax = sorted(inspax,key = lambda x:x[1])
                inspax = inspax[:noc]

                pp.append(inspax[0][1])
                for ii in inspax:
                    closest.append(ii[0])
                    
                pp.append(closest)
                #closout.write('%d\t%d\t%d\t%d\t%d\t[' %(pp[0],pp[1],pp[2],pp[4],pp[5]))
                #for ii in pp[6]:
                    #closout.write('%d,' %ii)
                #closout.write(']\n')
                #print pp
                ctpp+=1
            #closout.close()         
            print "Done with closest points\nCalculating Normals"

            for pp in pos:
                pp.append(findNormal(pp,pos))

            toFile(pos,'AllSurfaceInfo.txt')

        else:
            pos = []

            opf = 1

            if opf == 0:
                with open('AllSurfaceInfo.txt','r') as rd:
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
                        pos.append(temp[0:4]+temp[-3:])

            posN = []
            if opf == 1:
                with open('AllSurfaceInfoCenNew_tt5.txt','r') as rd:
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
                        posN.append(temp)

            opf = 2
            
            if opf == 2:
                with open('AllSurfaceInfoCen_indices.txt','r') as rd:
                    for lineL in rd:
                        temp = []
                        lineL = lineL.split('[')
                        for lin in lineL:
                            tt = []
                            lin = lin.replace(']','\t')
                            lin = lin.replace(',','\t')
                            lin = lin.split('\t')
                            for ik in lin:
                                try:
                                    tt.append(int(ik))
                                except ValueError:
                                    try:
                                        tt.append(float(ik))
                                    except ValueError: pass
                            temp.append(tt)
                        pos.append(temp)
       
            print "Done reading in file"
            print "Plotting Surface points\nCalculating Centroids"

            doCen = 0

            tt = 1

            if doCen==1:
                plusMin = 0
                
                if plusMin == 0:
                    for tt in [3,4]:
                        print "New",tt
                        
                        for iip in pos:
                            
                            #if iip[3]in[1997,2487,2489,2493,2495,2936,2940,2942,3348,3713,5655]:

                            #print iip[-1]
                            if tt > 3:
                                if type(iip[-1])==list:
                                    print 'list'
                                    del(iip[-1])
                            
                            Nnom = iip[-3:]

                            x2p = np.add(np.dot(tt,Nnom),iip[0:3])
                            x2m = np.add(np.dot(-1*tt,Nnom),iip[0:3])
                            #print Nnom
                            #print x2p
                            #print x2m
                            
                            ctp = 0; ctm = 0;
                            mx = []; my = []; mz = [];
                            px = []; py = []; pz = [];
                            
                            tmx = 0; tmy = 0; tmz = 0;
                            tpx = 0; tpy = 0; tpz = 0;                            
                            for pc in pos:
                                if (x2p[0]==pc[0] and x2p[1]==pc[1]) and (pc[2]<x2p[2] or pc[2]>x2p[2]):
                                    if tpz==0:
                                        pz.append(pc[2]-x2p[2])
                                        kp = 'keep'
                                    else:
                                        if pc[2]==tpz+1:                                                
                                            kp = 'pass'
                                        else:
                                            pz.append(pc[2]-x2p[2])
                                            kp = 'keep'
                                    tpz = pc[2]
                                    #print "pz: ", pc[0:3],tpz,pz,' ',kp
                                if (x2p[1]==pc[1] and x2p[2]==pc[2]) and (pc[0]<x2p[0] or pc[0]>x2p[0]):
                                    if tpx == 0:
                                        px.append(pc[0]-x2p[0])
                                        kp = 'keep'
                                    else:
                                        if pc[0]==tpx+1:                                    
                                            kp = 'pass'
                                        else:
                                            px.append(pc[0]-x2p[0])
                                            kp = 'keep'
                                    tpx = pc[0]
                                    #print "px: ", pc[0:3],tpx,px,' ',kp
                                if (x2p[0]==pc[0] and x2p[2]==pc[2]) and (pc[1]<x2p[1] or pc[1]>x2p[1]):
                                    if py == 0:
                                        py.append(pc[1]-x2p[1])
                                        kp = 'keep'
                                    else:
                                        if pc[1]==tpy+1:
                                            kp = 'pass'
                                        else:
                                            py.append(pc[1]-x2p[1])
                                            kp = 'keep'
                                    tpy = pc[1]
                                    #print "py: ", pc[0:3],tpy,py,' ',kp
                                if (x2m[0]==pc[0] and x2m[1]==pc[1]) and (pc[2]<x2m[2] or pc[2]>x2m[2]):
                                    if tmz == 0:
                                        mz.append(pc[2]-x2m[2])
                                        kp = 'keep'
                                    else:
                                        if pc[2]==tmz+1:
                                            kp = 'pass'
                                        else:
                                            mz.append(pc[2]-x2m[2])
                                            kp = 'keep'
                                    tmz = pc[2]
                                    #print "mz: ", pc[0:3],tmz,mz,' ',kp
                                if (x2m[1]==pc[1] and x2m[2]==pc[2]) and (pc[0]<x2m[0] or pc[0]>x2m[0]):
                                    if tmx == 0:
                                        mx.append(pc[0]-x2m[0])
                                        kp = 'keep'
                                    else:
                                        if pc[0]==tmx+1:                            
                                            kp = 'pass'
                                        else:
                                            mx.append(pc[0]-x2m[0])
                                            kp = 'keep'
                                    tmx = pc[0]
                                    #print "mx: ", pc[0:3],tmx,mx,' ',kp
                                if (x2m[0]==pc[0] and x2m[2]==pc[2]) and (pc[1]<x2m[1] or pc[1]>x2m[1]):
                                    if tmy == 0:
                                        my.append(pc[1]-x2m[1])
                                        kp = 'keep'
                                    else:
                                        if pc[1]==tmy+1:
                                            kp = 'pass'
                                        else:
                                            my.append(pc[1]-x2m[1])
                                            kp = 'keep'
                                    tmy = pc[1]
                                    #print "my: ", pc[0:3],tmy,my,' ',kp

                            if signChange(pz):ctp += 1
                            if signChange(px):ctp += 1
                            if signChange(py):ctp += 1

                            if signChange(mz): ctm += 1
                            if signChange(mx): ctm += 1
                            if signChange(my): ctm += 1
                            
                            if ctp < ctm:
                                mult = 1
                            elif ctm < ctp:
                                mult = -1
                            else:
                                mult = 0

                            Nnorm = [mult*Nnom[0],mult*Nnom[1],mult*Nnom[2]]

                            #print iip[3],iip[0],"multiply by ",ctp,ctm,mult,Nnorm
                            #iip.append(ctp)
                            #iip.append(ctm)
                            iip.append(Nnorm)
                            
                        toFile(pos,'AllSurfaceInfoCenNew%d.txt' %tt,'float')
                        
                    '''
                    xy = []
                    xz = []
                    yz = []
                    for pc in pos:
                        if iip[0]==pc[0] and iip[1]==pc[1]:
                            xy.append(pc[3])
                        if iip[0]==pc[0] and iip[2]==pc[2]:
                            xz.append(pc[3])
                        if iip[0]==pc[1] and iip[2]==pc[2]:
                            yz.append(pc[3])

                    iip.append(xy)
                    iip.append(xz)
                    iip.append(yz)
                    '''
                else:           
                    step = 30
                    #step = 10

                    Xx = [[],[],[]]

                    for pp in pos:
                        Xx[0].append(pp[0])
                        Xx[1].append(pp[1])
                        Xx[2].append(pp[2])
                        
                    print len(Xx[0]),len(pos)

                    #print Xx.index(4)
                    oldx = min(Xx[0])

                    maxx = max(Xx[0])
                    maxy = max(Xx[1])
                    maxz = max(Xx[2])
                    
                    for iip in pos:
                        startx = iip[0]; starty = iip[1]; startz = iip[2]
                        temp = [[],[],[]]

                        gh = startx-step/2
                        ghm = startx+1+step/2

                        #print gh,gm
                        if iip[3]==0 or startx!=oldx:
                            oldx = startx
                            print "finding index: ",oldx
                            if gh in Xx[0]:
                                sti = Xx[0].index(gh)
                            else:
                                sti = 0

                            if ghm in Xx[0]:
                                ste = Xx[0].index(ghm)
                            else:
                                ste = len(Xx[0])
                            stc = iip[3]

                        #print iip[3],startx,oldx," => searching [",sti,",",ste,"]"
                        
                        grid = [[startx-step/2,startx+step/2],[starty-step/2,starty+step/2],[startz-step/2,startz+step/2]]
                        #grid = [[startx,startx+step],[starty,starty+step],[startz,startz+step]]
                        '''
                        stfi = stc
                        stfe = ste
                        if startx > maxx-step:
                            stfi = sti
                            stfe = stc
                            grid[0]=[startx-step,startx]
                        if starty > maxy-step:
                            grid[1]=[starty-step,starty]
                        if startz > maxz-step:
                            grid[2]=[startz-step,startz]
                        '''
                        print iip[3],grid," => searching [",sti,",",ste,"]"
                        for ic in range(sti,ste):
                            ppc = pos[ic]
                            if ppc[0]>=grid[0][0] and ppc[0] <= grid[0][1] and ppc[1]>=grid[1][0] and ppc[1] <= grid[1][1] and ppc[2]>=grid[2][0] and ppc[2] <= grid[2][1]:
                                temp[0].append(ppc[0])
                                temp[1].append(ppc[1])
                                temp[2].append(ppc[2])

                        centroid = [np.mean(temp[0]).item(),np.mean(temp[1]).item(),np.mean(temp[2]).item()]
                        #print centroid

                        iip.append(centroid)
                    
                #toFile(pos,'AllSurfaceInfoCen.txt','float')
                
            print "Preparing points to plot"
            
            shdplot = 1
                    
            if shdplot==1:
                gdd = 0
                #grid = [[0,20],[78,98],[35,55]]
                #grid = [[0,30],[25,55],[30,60]]
                #grid = [[30,60],[40,70],[30,60]]
                #grid = [[0,zsize],[0,ysize],[0,xsize]]
                grid = [[65,95],[40,70],[30,60]]
                
                #for xsh in range(0,zsize,30):
                #    grid[0] = [xsh,xsh+30]
                #    for ysh in range(0,ysize,30):
                #        grid[1] = [ysh,ysh+30]
                #        for zsh in range(0,xsize,xsize):
                #            grid[2] = [zsh,xsize]
                #            

                '''
                ctr = 0
                temp = [[],[],[],[]]
                for pp in pos:  
                    if pp[0]>=grid[0][0] and pp[0] <= grid[0][1] and pp[1]>=grid[1][0] and pp[1] <= grid[1][1] and pp[2]>=grid[2][0] and pp[2] <= grid[2][1]:
                        #print pp
                        temp[0].append(pp[0])
                        temp[1].append(pp[1])
                        temp[2].append(pp[2])
                        temp[3].append(0.98)
                        ctr+=1
                '''
                
                print "Finding tops"

                tops = []

                '''
                for ii in pos:
                    if ii[0][3] not in tops:
                        iz = ii[0][2]
                        xy = ii[1]
                        group = []
                        grp = 0
                        grplist = [0]
                        lastz = pos[ii[1][0]][0][2]
                        
                        for hc in xy:
                            newz = pos[hc][0][2]
                            if abs(newz-lastz)>1:
                                grp += 1
                                grplist.append(grp)
                            group.append(grp)
                            lastz = newz


                        gdgrp = grplist[-1]
                        idgd = group.index(gdgrp)

                        for jj in range(idgd,len(group)):
                            tops.append(xy[jj])                        
                '''

                for pp in posN:
                    Nnom = pp[-3:]
                    if Nnom[2]<0: pass
                    else:
                        tops.append(pp[3])
                        
                print "Plotting Surface"

                temp = [[],[],[],[]]
                for tp in tops:#range(len(pos)):
                    pp = pos[tp][0]
                    if pp[0]>=grid[0][0] and pp[0] <= grid[0][1] and pp[1]>=grid[1][0] and pp[1] <= grid[1][1] and pp[2]>=grid[2][0] and pp[2] <= grid[2][1]:
                        temp[0].append(pp[0])
                        temp[1].append(pp[1])
                        temp[2].append(pp[2])
                        temp[3].append(0.98)
                        
                fig = plt.figure()
                ax = Axes3D(fig)

                ax.set_xlim([0,ysize])
                ax.set_ylim([0,ysize])
                ax.set_zlim([0,ysize])
                
                ax.scatter(temp[0],temp[1],temp[2],s=100*temp[3],c= colfc,marker='o',edgecolor = colbd)
                #ax.view_init(elev=0,azim=0)
                ax.view_init(elev=40,azim=50)
                #ax.view_init(elev=90,azim=50)

                #print "Plotting Normals"
                
                #grid = [[15,20],[78,85],[35,55]]
                #grid = [[0,20],[78,98],[35,55]]
                grid = [[85,95],[50,70],[30,60]]
                
                #print "saving"
                fig.savefig("A2CheckSurfaceGrids%d.png" %gdd)
                close()
                #print gdd,"=> x: ",grid[0],"y: ",grid[1],"z: ",grid[2]
                gdd += 1

                plotnorm = True
                
                if plotnorm:

                    print "Plotting normals"
                    
                    tt = 5

                    ax2 = Axes3D(fig)

                    ax2.set_xlim([0,ysize])
                    ax2.set_ylim([0,ysize])
                    ax2.set_zlim([0,ysize])

                    #grid = [[0,30],[25,55],[30,40]]

                    ct = 0
                    for tp in tops:
                        tempA = [[],[],[]]
                        pp = posN[tp]
                        #print pp
                        Nnom = pp[-3:]
                        if True:#abs(Nnom[2])==1:#abs(Nnom[0])==1 and abs(Nnom[1])==0 and abs(Nnom[2])==0:#sum(np.abs(Nnom))!=1 and sum(np.abs(Nnom))!=0:# pp[3]in[2420]:#
                            #print (pp[-3:])
                            if pp[0]>=grid[0][0] and pp[0] <= grid[0][1] and pp[1]>=grid[1][0] and pp[1] <= grid[1][1] and pp[2]>=grid[2][0] and pp[2] <= grid[2][1]:
                                tempA[0].append(pp[0])
                                tempA[1].append(pp[1])
                                tempA[2].append(pp[2])
                                
                                #if abs(Nnom[2])==0:
                                    #print pp[3], Nnom
                                    
                                #centroid = pp[-3:]
                                #print pp
                                #print Nnom
                                
                                nom = np.dot(1/np.linalg.norm(Nnom),Nnom)
                                #nom = Nnom[:]
                                #print nom

                                #print np.linalg.norm(nom)

                                #for tt in range(Ttend):
                                #x2p = np.add(np.dot(tt,Nnom),pp[0:3])
                                #x2m = np.add(np.dot(-1*tt,Nnom),pp[0:3])

                                #print "Point:",pp[0:3]
                                #print "Positive vector:",x2p
                                #print "Negative vector:",x2m

                                #x2pd = np.subtract(x2p,centroid)
                                #x2md = np.subtract(x2m,centroid)

                                #print x2pd
                                #print x2md

                                #pdnorm = x2pd[0]**2+x2pd[1]**2+x2pd[2]**2
                                #mdnorm = x2md[0]**2+x2md[1]**2+x2md[2]**2

                                #print centroid
                                #print "distance +ve = ",pdnorm
                                #print "distance -ve = ",mdnorm
                                
                                #if pdnorm>mdnorm:
                                x2p = np.add(np.dot(tt,nom),pp[0:3])
                                x2=[x2p[0].item(),x2p[1].item(),x2p[2].item()]
                                #else:
                                #    x2m = np.add(np.dot(-1*tt,nom),pp[0:3])
                                #    x2=[x2m[0].item(),x2m[1].item(),x2m[2].item()]
                                
                                tempA[0].append(x2[0])
                                tempA[1].append(x2[1])
                                tempA[2].append(x2[2])

                                #ax2 = Axes3D(fig)

                                #ax2.set_xlim([0,ysize])
                                #ax2.set_ylim([0,ysize])
                                #ax2.set_zlim([0,ysize])
                                #ax2.view_init(elev=30,azim=225)
                                #ax2.view_init(elev=40,azim=50)
                                
                                ax2.plot(tempA[0],tempA[1],tempA[2],c= 'black')
                                #ax2.scatter(temp[0],temp[1],temp[2],s=100*temp[3],c= colfc,marker='o',edgecolor = colbd)
                                #ax2.scatter(centroid[0],centroid[1],centroid[2],s=98,c= 'red',marker='*',edgecolor = colbd) 
                                
                                #fig.savefig("A2CheckSurfaceGrid2ArrID_%d.png" %pp[3])
                                close()
                                ct += 1
                                #print pp[3]#, ct, '/', ctr
                                
                    #ax2.view_init(elev=40,azim=50)
                    #ax2.view_init(elev=0,azim=250)
                    #ax2.view_init(elev=0,azim=0)
                    #ax2.view_init(elev=20,azim=150)
                    ax2.view_init(elev=20,azim=50)
                    
                    ax2.scatter(temp[0],temp[1],temp[2],s=100*temp[3],c= colfc,marker='o',edgecolor = colbd)                        
                    fig.savefig("A2CheckSurfaceGrid.png")

                    print "saving"
                    #fig.savefig("CheckSurfaceGridArr.png")
                    
                    close()
        
        '''
        randIDS = []
        for ii in range(100):
            randIDS.append(random.randint(0,lenpos-1))

        for idd in randIDS:

            temp = [[],[],[],[]]

            temp[0].append(pos[idd][0])
            temp[1].append(pos[idd][1])
            temp[2].append(pos[idd][2])
            temp[3].append(pos[idd][3])

            #print pos[idd][6]
            
            for pp in pos[idd][6]:
                temp[0].append(pos[pp][0])
                temp[1].append(pos[pp][1])
                temp[2].append(pos[pp][2])
                temp[3].append(pos[pp][3])
            
            fig = plt.figure()
            ax = Axes3D(fig)

            ax.set_xlim([0,ysize])
            ax.set_ylim([0,ysize])
            ax.set_zlim([0,ysize])

            ax.scatter(temp[0],temp[1],temp[2],s=100*temp[3],c= colfc,marker='o',edgecolor = colbd)
            ax.view_init(elev=50,azim=50)

            #plt.show()
            #print ai
            print idd
            fig.savefig("ID_%d_closestpts.png" %idd)
            close()
        '''
        
        '''  
        for zLIM in [10,15,20,30,45,61]:

            temp = [[],[],[],[]]

            for pp in pos:
                if pp[2] <= zLIM:
                    temp[0].append(pp[0])
                    temp[1].append(pp[1])
                    temp[2].append(pp[2])
                    temp[3].append(pp[3])
            
            fig = plt.figure()
            ax = Axes3D(fig)

            ax.set_xlim([0,ysize])
            ax.set_ylim([0,ysize])
            ax.set_zlim([0,ysize])

            ax.scatter(temp[0],temp[1],temp[2],s=100*temp[3],c= colfc,marker='o',edgecolor = colbd)
            ax.view_init(elev=50,azim=50)

            #plt.show()
            print ai
            print zLIM
            fig.savefig("zLIM_%d_justSurface.png" %zLIM)
            close()
        '''

    ctt = ctt + 1
    
    print "%d===>Done" %SNo
