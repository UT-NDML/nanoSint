import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import struct
import time
import math as m
import sys
from os.path import isfile
import multiprocessing as mp

sizep = int(sys.argv[1]) #number of processors

def worker(num):

    def makeSphere(xC,yC,zC,R):
        
        lowb = [m.floor(xC-R),m.floor(yC-R),m.floor(zC-R)]
        diam = m.ceil(2*R)
        
        spherePoints = [[],[],[]]
        
        #print ('Circle properties-->',xC,yC,R)
        for ix in range(lowb[0],lowb[0]+diam+3):
            for iy in range(lowb[1],lowb[1]+diam+3):
                for iz in range(lowb[2],lowb[2]+diam+3):
                    if (ix-xC)**2+(iy-yC)**2+(iz-zC)**2<=R**2:
                        spherePoints[0].append(ix)
                        spherePoints[1].append(iy)
                        spherePoints[2].append(iz)
        
        return spherePoints
    
    def sortPlot(toPlot):
        newList = []
        for ix,iy,iz,c in zip(toPlot[0],toPlot[1],toPlot[2],toPlot[3]):
            newList.append([ix,iy,iz,c])
        newList.sort(key = lambda x:(x[0],x[1],x[2]))
        
        newL = [[],[],[],[]]
        for xx in newList:
            for nn in range(4): newL[nn].append(xx[nn])
        
        return newL
    
    def addArrays(*x, listlen=2):
        allOut = []
        for xx in x:
            newOut = []
            for ii in range(listlen):
                newOut += list(xx[ii])
            allOut.append(np.array(newOut))
            
        if len(allOut)==1: return allOut[0]
        else: return allOut
        
    def props(x,y,z,debug=False,Center = [0]):
        cenIn = len(Center)
        if cenIn != 3: Center = np.array([(max(x)+min(x))/2, (max(y)+min(y))/2, (max(z)+min(z))/2])
        mid = (x-Center[0])**2 + (y-Center[1])**2 + (z-Center[2])**2
        #rad = max(max(x)-min(x),max(y)-min(y),max(z)-min(z))/2 #sqrt(2)*maxdistance/2
        #print('low ', rad, end = ' ')
        #rad = 1.415*max(max(x)-min(x),max(y)-min(y),max(z)-min(z))/2 #sqrt(2)*maxdistance/2
        #print('high', rad, end = ' ')
        rad = m.ceil(np.sqrt(max(mid)))
        #print('current', rad, end = ' ')
        
        if debug:
            radSq = rad*rad
            check = mid >= radSq
            
            if check.any():
                print("Not all pixels contained in sphere")
            else:
                print("All pixels contained in sphere")
            
        if cenIn != 3: return rad,Center
        else: return rad
    
    def lineEQ(t,p0,p1):
        return p0 - t*(p0 - p1)
        
    def planeEQ(N,p,pi):
        return np.dot(N,pi-p)
    
    #%%
    def plotNeck(ovlp,figName,elevN=40,azimA=50):
        fig = plt.figure()
        ax = Axes3D(fig)
    
        ax.set_xlim([0,ysize])
        ax.set_ylim([0,ysize])
        ax.set_zlim([0,ysize])
        
        ovlp = sortPlot(ovlp)
        #ax.scatter(ovlp[0],ovlp[1],ovlp[2],s=100*np.array(ovlp[3]),c= colfc,marker='o',edgecolor = colbd)
        ax.scatter(ovlp[0],ovlp[1],ovlp[2],s=9.8,c= colfc,marker='o',edgecolor = colbd)
        #ax.scatter(cenN[0],cenN[1],cenN[2],s=100*radN,c= 'green',marker='o',edgecolor = colbd)
        ax.view_init(elev=elevN,azim=azimA)
    
        fig.savefig(figName)
        plt.close() 
        
    def testParticles(p1,p2):
        x = [[],[]]; y = [[],[]]; z = [[],[]]; v = [[],[]]
        t1 = Ed[p1] >= etaCut
        idd = Vt&t1
        ip = 0; x[ip] = Xd[idd]; y[ip] = Yd[idd]; z[ip] = Zd[idd]; v[ip] = Vd[idd] 
        rad0, cen0 = props(x[0],y[0],z[0])
        t1 = Ed[p2] >= etaCut
        idd = Vt&t1
        ip = 1; x[ip] = Xd[idd]; y[ip] = Yd[idd]; z[ip] = Zd[idd]; v[ip] = Vd[idd]
        rad1, cen1 = props(x[1],y[1],z[1])
        titleName, ovlpA = inContact(x,y,z,v,[rad0,rad1],[cen0,cen1],False,True,p1=p1,p2=p2)
        if titleName:
            #unpack ovlp
            #ovlp,radN,cenN,finPos,overlap01,lenN,neckAll = ovlpA
            ovlp, radN, overlap01, lenN, neckAll = ovlpA
            print('Eta cut -->',etaCut,'Touching\nNeck radius =','{:0.4f}'.format(radN))
            #plotNeck(overlap01,'cen2cen'+str(p1)+'_and_'+str(p2)+'.png')
            plotNeck(ovlp,'ATEST_b4Narrowed'+str(p1)+'_and_'+str(p2)+'.png')
            #plotNeck(finPos,'ATEST_lim'+str(etaCut)+'_NarrowedNeck'+str(p1)+'_and_'+str(p2)+'.png')
            plotNeck(neckAll,'ATEST_lim'+str(etaCut)+'_NarrowedNeckRegion'+str(p1)+'_and_'+str(p2)+'.png')
        else:
            print('Eta cut -->',etaCut,'Not touching')
        return titleName, ovlpA
    
    def inContact(x,y,z,v,RAD,CEN,ovlp01 = False,debug=False,**info):
        global pli
    
        def mapPlane(x,y,z):
            global pli
            p = np.array([x,y,z]) - pli
            return np.dot(N,p)
        
        def narrowNeck1plane(left = True, debug = False):
            global pli
            if left: 
                ii = ilft
                tt = tl
                step = delT
                if debug: print('left index is ', ii)
            else: 
                ii = irgt
                step = -delT
                tt = tr
                if debug: print('right index is ', ii)
            
            xl,yl,zl,vl = x[ii][togP[ii]],y[ii][togP[ii]],z[ii][togP[ii]],v[ii][togP[ii]]
            
            while tl<=tt<=tr:
                if debug: print(tt,"size ->", len(xl))
                tt += step
                pli = lineEQ(tt,cen0,cen1)
                arr0 = np.array(list(map(mapPlane,xl,yl,zl)))
                if left: pog = arr0 >= 0
                else: pog = arr0 <= 0
                if not pog.any():
                    return tt,xl,yl,zl,vl
                xl,yl,zl,vl = xl[pog],yl[pog],zl[pog],vl[pog]
            return tt,xl,yl,zl,vl
    
        def narrowNeck(debug = True):
            global pli
            tlt = tl
            trt = tr
            
            leftMove = True
            rightMove = True
            
            if debug: print('left index is ', ilft,'right index is ', irgt)
            
            ii = ilft; xl,yl,zl,vl = x[ii][togP[ii]],y[ii][togP[ii]],z[ii][togP[ii]],v[ii][togP[ii]]
            ii = irgt; xr,yr,zr,vr = x[ii][togP[ii]],y[ii][togP[ii]],z[ii][togP[ii]],v[ii][togP[ii]]
            
            while tlt < trt:
                arrL = [[],[]]; arrR = [[],[]]; pog = [0,0];
                if debug: print('left:','{:0.4f}'.format(tlt),"size ->", len(xl),'|| right:','{:0.4f}'.format(trt),"size ->", len(xr))
                
                if leftMove:
                    tlt += delT; 
                    pli = lineEQ(tlt,cen0,cen1); 
                    ii = ilft; arrL[ii] = np.array(list(map(mapPlane,xl,yl,zl)))
                    ii = irgt; arrL[ii] = np.array(list(map(mapPlane,xr,yr,zr)))
                    if not rightMove:
                        pog[ilft] = arrL[ilft]>=0; pog[irgt] = arrL[irgt]>=0;
                        if not pog[ilft].any() or not pog[irgt].any():
                            tlt -= delT; 
                            if debug: print('(trt - tlt)*(rad0+rad1) = ',(trt - tlt)*(rad0+rad1))
                            if (trt - tlt)*(rad0+rad1) > 1.5: 
                                #true if planes are more than one pixel apart in which case
                                #there is empty space between planes and they aren't touching
                                return False, [xl,xr],[yl,yr],[zl,zr],[vl,vr]
                            else: break
                        else: 
                            ii = ilft; xl,yl,zl,vl = xl[pog[ii]],yl[pog[ii]],zl[pog[ii]],vl[pog[ii]]
                            ii = irgt; xr,yr,zr,vr = xr[pog[ii]],yr[pog[ii]],zr[pog[ii]],vr[pog[ii]]
                            continue
                        
                if rightMove:
                    trt -= delT
                    pli = lineEQ(trt,cen0,cen1); 
                    ii = ilft; arrR[ii] = np.array(list(map(mapPlane,xl,yl,zl)))
                    ii = irgt; arrR[ii] = np.array(list(map(mapPlane,xr,yr,zr)))
                    if not leftMove:
                        pog[irgt] = arrR[irgt]<=0; pog[ilft] = arrR[ilft]<=0;
                        if not pog[ilft].any() or not pog[irgt].any():
                            trt += delT
                            if debug: print('(trt - tlt)*(rad0+rad1) = ',(trt - tlt)*(rad0+rad1))
                            if (trt - tlt)*(rad0+rad1) > 1.5:
                                return False, [xl,xr],[yl,yr],[zl,zr],[vl,vr]
                            else: break
                        else:
                            ii = ilft; xl,yl,zl,vl = xl[pog[ii]],yl[pog[ii]],zl[pog[ii]],vl[pog[ii]]
                            ii = irgt; xr,yr,zr,vr = xr[pog[ii]],yr[pog[ii]],zr[pog[ii]],vr[pog[ii]]
                            continue
                
                lft = [0,0]; rgt = [0,0]
                lft[ilft] = arrL[ilft]>=0; lft[irgt] = arrL[irgt]>=0
                rgt[ilft] = arrR[ilft]<=0; rgt[irgt] = arrR[irgt]<=0
                pog[ilft] = lft[ilft]&rgt[ilft]; pog[irgt] = lft[irgt]&rgt[irgt]
                
                if not pog[ilft].any() and pog[irgt].any():
                    leftMove = False
                    tlt -= delT; 
                    continue
                elif not pog[irgt].any() and pog[ilft].any():
                    rightMove = False
                    trt += delT
                    continue
                elif not pog[ilft].any() and not pog[irgt].any():
                    tlt -= delT; trt += delT
                    if debug: print('(trt - tlt)*(rad0+rad1) = ',(trt - tlt)*(rad0+rad1))
                    if (trt - tlt)*(rad0+rad1) > 1.5:
                        return False, [xl,xr],[yl,yr],[zl,zr],[vl,vr]
                    else: break                
                    #empty space between particles and so particles aren't touching
                    return False, [xl,xr],[yl,yr],[zl,zr],[vl,vr]
    
                ii = ilft; xl,yl,zl,vl = xl[pog[ii]],yl[pog[ii]],zl[pog[ii]],vl[pog[ii]]
                ii = irgt; xr,yr,zr,vr = xr[pog[ii]],yr[pog[ii]],zr[pog[ii]],vr[pog[ii]]
            
            return True, [xl,xr],[yl,yr],[zl,zr],[vl,vr]
        
        def neckSweep(debug = False):
            def findIDS(radA):
                temp = list(set(radA))
                temp.remove(min(temp))
                radt = min(temp)
                radm = min(radA)
                
                tempA = np.array(radA)
                
                temp = np.where(tempA == radt)[0]
                id0 = min(temp); id1 = max(temp)
                temp = np.where(tempA == radm)[0]
                id0 = min(id0,min(temp))
                id1 = max(id1,max(temp))
                
                return radm, radt, id0, id1
            
            def pointsBetPlanes(tsw,tswr):
                global pli
                pli = lineEQ(tswr,cen0,cen1)
                ii = ilft; arrR[ii] = np.array(list(map(mapPlane,xl,yl,zl)))
                ii = irgt; arrR[ii] = np.array(list(map(mapPlane,xr,yr,zr)))
                
                pli = lineEQ(tsw,cen0,cen1)
                ii = ilft; arrL[ii] = np.array(list(map(mapPlane,xl,yl,zl)))
                ii = irgt; arrL[ii] = np.array(list(map(mapPlane,xr,yr,zr)))
                
                lft = [0,0]; rgt = [0,0]
                lft[ilft] = arrL[ilft]>=0; lft[irgt] = arrL[irgt]>=0
                rgt[ilft] = arrR[ilft]<=0; rgt[irgt] = arrR[irgt]<=0
                pog[ilft] = lft[ilft]&rgt[ilft]; pog[irgt] = lft[irgt]&rgt[irgt]            
                
                ii = ilft; xlt,ylt,zlt,vlt = xl[pog[ii]],yl[pog[ii]],zl[pog[ii]],vl[pog[ii]]
                ii = irgt; xrt,yrt,zrt,vrt = xr[pog[ii]],yr[pog[ii]],zr[pog[ii]],vr[pog[ii]]
                
                return [xlt,xrt],[ylt,yrt],[zlt,zrt],[vlt,vrt]
            
            global pli
    
            #radA = [[],[],[]]
            apos = [[[],[]],[[],[]],[[],[]],[[],[]]]
            ctN = 0
            rt0 = 0; 
            radAs = []
            
            ii = ilft; xl,yl,zl,vl = x[ii][togP[ii]],y[ii][togP[ii]],z[ii][togP[ii]],v[ii][togP[ii]]
            ii = irgt; xr,yr,zr,vr = x[ii][togP[ii]],y[ii][togP[ii]],z[ii][togP[ii]],v[ii][togP[ii]]
            
            tsw = tl + delT
            
            while tsw <= tr - delT:
                arrL = [[],[]]; arrR = [[],[]]; pog = [0,0];
                tswr = tsw + delT
                pos = pointsBetPlanes(tsw,tswr)
                
                if len(pos[0][1]) != 0: rt0 += 1
                if (rt0 == 1 and len(pos[0][0]) == 0) or (len(pos[0][1]) == 0 and len(pos[0][0]) == 0): return False,0,0,[],[],[],[]
                
                if debug: print("planes: <",'{:0.4f}'.format(tsw),'{:0.4f}'.format(tswr),"> #left in:",len(pos[0][0]),"#right in:",len(pos[0][1]),end = '')           
                
                if rt0 > 0 and len(pos[0][0])==0:
                    if debug: print('\n')
                    break
                if rt0 > 0 and len(pos[0][0])!=0:
                    for jj in range(4):
                        for ij in range(2):
                            apos[jj][ij] = addArrays([apos[jj][ij],pos[jj][ij]])
                    sweepPos = addArrays(pos[0],pos[1],pos[2],pos[3])
                    radS = props(sweepPos[0],sweepPos[1],sweepPos[2],Center = pli)
                    radAs.append(radS)
                    ctN += 1
                    if debug: print(" radius of slice:", radS)
                else: 
                    if debug: print()
                #sweepPos = addArrays(pos[0],pos[1],pos[2],pos[3])
                #radS = props(sweepPos[0],sweepPos[1],sweepPos[2],Center = pli)
                #if debug: print("radius of slice:", radS)
                #radA[0].append(tsw); radA[1].append(tswr); radA[2].append(radS)
                tsw = tswr
    
            #radm, radt, id0, id1 = findIDS(radA[2])
            #radN = (radm+radt)/2
            #ctN = id1- id0 + 1
            #tsw = radA[0][id0]
            #tswr = radA[1][id1]
            #apos = pointsBetPlanes(tsw,tswr) 
            
            radN = np.mean(radAs)
            
            return True,radN,ctN,apos[0],apos[1],apos[2],apos[3]
        
        overlap = []; overlap01 = []
        rad0, rad1 = RAD 
        cen0, cen1 = CEN   
        ilft,irgt = 0,1
        
        if info:
            if 'p1' in info.keys(): p1=info['p1']
            if 'p2' in info.keys(): p2=info['p2']
        
        N = cen0 - cen1
        cenD = np.sqrt(sum(N**2))
        
        if cenD <= rad0 + rad1:
            if planeEQ(N,cen1,cen0)>0:
                cen1,cen0 = cen0,cen1
                rad1,rad0 = rad0,rad1
                ilft,irgt = 1,0
            
            if ovlp01:
                arr1 = [[],[]]; arr2 = [[],[]]
                pli = cen0
                ii = ilft; arr1[ii] = np.array(list(map(mapPlane,x[ii],y[ii],z[ii])))
                ii = irgt; arr1[ii] = np.array(list(map(mapPlane,x[ii],y[ii],z[ii])))
                pli = cen1
                ii = ilft; arr2[ii] = np.array(list(map(mapPlane,x[ii],y[ii],z[ii])))
                ii = irgt; arr2[ii] = np.array(list(map(mapPlane,x[ii],y[ii],z[ii])))
                lft = [0,0]; rgt = [0,0]
                lft[ilft] = arr1[ilft]>=0; lft[irgt] = arr1[irgt]>=0
                rgt[ilft] = arr2[ilft]<=0; rgt[irgt] = arr2[irgt]<=0
                tog = addArrays(lft)&addArrays(rgt)
                newx, newy, newz, newv = addArrays(x,y,z,v)        
                overlap01 = [newx[tog], newy[tog], newz[tog], newv[tog]]
                
            ovl = rad0 + rad1 - cenD
            tr = rad0/(rad0 + rad1)
            ovl = ovl/(rad0 + rad1)
            tl = tr - ovl
            tl = 0; tr = 1
            pl1 = lineEQ(tl,cen0,cen1)
            pl2 = lineEQ(tr,cen0,cen1)
            
            arr1 = [[],[]]; arr2 = [[],[]]
            pli = pl1
            ii = ilft; arr1[ii] = np.array(list(map(mapPlane,x[ii],y[ii],z[ii])))
            ii = irgt; arr1[ii] = np.array(list(map(mapPlane,x[ii],y[ii],z[ii])))
            
            pli = pl2
            ii = ilft; arr2[ii] = np.array(list(map(mapPlane,x[ii],y[ii],z[ii])))
            ii = irgt; arr2[ii] = np.array(list(map(mapPlane,x[ii],y[ii],z[ii])))
            
            lft = [0,0]; rgt = [0,0]
            lft[ilft] = arr1[ilft]>=0; lft[irgt] = arr1[irgt]>=0
            rgt[ilft] = arr2[ilft]<=0; rgt[irgt] = arr2[irgt]<=0
    
            tog = addArrays(lft)&addArrays(rgt)
    
            newx, newy, newz, newv = addArrays(x,y,z,v)
            
            overlap = [newx[tog], newy[tog], newz[tog], newv[tog]]
            
            delT = 1/(rad0+rad1)
                
            togP = [0,0]; togP[ilft] = lft[ilft]&rgt[ilft]; togP[irgt] = lft[irgt]&rgt[irgt]
            
            if togP[0].any() and togP[1].any():
                
                #testing with moving 1 plane at a time
                #xNar = [[],[]]; yNar = [[],[]]; zNar = [[],[]]; vNar = [[],[]];
                #ii = ilft; xNar[ii], yNar[ii], zNar[ii], vNar[ii] = narrowNeck1plane()
                #radL, cenL = props(xNar[ii], yNar[ii], zNar[ii])
                #ii = irgt; xNar[ii], yNar[ii], zNar[ii], vNar[ii] = narrowNeck1plane(False)
                #radR, cenR = props(xNar[ii], yNar[ii], zNar[ii])
                
                touchCheck,radN,lenN,xA,yA,zA,vA = neckSweep(debug)
                #touchCheck, xA, yA, zA, vA = narrowNeck(debug)
                if touchCheck:
                    neckAll = addArrays(xA,yA,zA,vA)
                    #finPos = addArrays(xA,yA,zA,vA)
                    #radN, cenN = props(finPos[0],finPos[1],finPos[2])
    
                    #if debug: print("\n")
                    #lenN,neckAll = neckSweep(debug)
                    #if debug: print("\n")
                    
                    print('p',p1,' rad = ',rad0,' and p',p2,' rad = ',rad1, sep = '',end = ' | ')
                    print('lengths:=',len(xA[0]),'and',len(xA[1]),end = ' | ')
                    print('cenD =',cenD,end = ' | ')
                    print('neck radius =','{:0.4f}'.format(radN),'neck length =',lenN)
                    #return True, [overlap, radN, cenN, finPos, overlap01, lenN, neckAll]
                    return True, [overlap, radN, overlap01, lenN, neckAll]
                else:
                    return False, overlap
            else:
                return False, overlap
        else:
            return False, overlap
    
    def calcRes(r0,r1,rN,L):
        reff = (r0+r1)/2
        frac = reff/rN + m.log(2*reff/rN)/m.pi
        res = frac*Rbulk
        A = m.pi*rN*rN
        return res*L*Sconv/A 
    
    #%%
    def writeVolList(maxV,fileName,Vdc = 15):
        filv = open(fileName,'w')
        for ix in range(maxV+1):
            for iy in range(ix+1,maxV+1):
                filv.write('V1 %d %d dc %d\n' %(ix,iy,Vdc))
        filv.close()
        
    def writeNodeList(nodeList,fileName,Vn1,Vn2,Vdc = 15,vfileName = 'VoltageFile.txt'):
        
        #fill.write('Effective Resistance Circuit\n')
    
        filv = open(vfileName,'w')
        for ix in Vn1[0]:
            for iy in Vn1[1]:
                filv.write('V1 %d %d dc %d\n' %(ix,iy,Vdc))
        for ix in Vn2[0]:
            for iy in Vn2[1]:
                filv.write('V1 %d %d dc %d\n' %(ix,iy,Vdc))
        filv.close()
           
        #fill.write('V1 %d %d dc %d\n' %(Vn1,Vn2,Vdc))
        fill = open(fileName,'w')
        for ct, res in enumerate(nodeList[2]):
            rt = 'R'+str(ct+1)
            fill.write('%s %d %d %f\n' %(rt,nodeList[0][ct],nodeList[1][ct],res))
        fill.write('.dc V1 %d %d 1\n' %(Vdc,Vdc))
        fill.write('.print dc i(V1)\n')
        fill.write('.end\n')
        fill.close()
    #%%
        
    def getRead(fileCheck):
        infoDic = {}
        with open(fileCheck,'r') as fild:
            for ii in fild.readlines():
                if ii[0] == 'p':
                    lin = ii.split('|')
                    fl1 = lin[0]; fl2 = lin[-1]
                    fl1 = fl1.replace('rad = ','')
                    fl1 = fl1.split()
                    keyN = fl1[0]+fl1[3]
                    fl2 = fl2.replace('neck radius','').replace('neck length','').replace('=','').split()
                    props = [eval(fl2[0]),eval(fl2[1])]
                    infoDic[keyN] = props
        return infoDic
    
    def dicExtract(xminA):
        outt = []
        for ii in xminA:
            outt.append(nodeDic[ii])
        return outt
    
    def getOpt(listO,op,giv=4):
        opL = []
        if op==0:
            optV = min(listO)
            for ii in range(giv):
                optV = optV + ii
                opL += list(np.array(iilist)[np.where(np.array(listO)==optV)])
        elif op==1:
            optV = max(listO)
            for ii in range(giv):
                optV = optV - ii
                opL += list(np.array(iilist)[np.where(np.array(listO)==optV)])
        return opL 
    
    def rhoAnalysis(ai):
    
        xbox = Xd[Vt]; ybox = Yd[Vt]; zbox = Zd[Vt]; vbox = Vd[Vt]
        idd = (xbox >= xlb) & (xbox <= xhb) & (ybox >= ylb) & (ybox <= yhb)
        vbox = vbox[idd]      
        sumrho = sum(vbox)
        
        if ai == 0:
            init_volume = (max(xbox) - min(xbox))*(max(ybox) - min(ybox))*(max(zbox) - min(zbox))
            init_rho_ratio = sum(Vt)/init_volume
            init_density = bulkR*init_rho_ratio     
            print("init_density =", init_density)
        
        print("sumrho =", sumrho)
        
        return
    #%%
    
    Rbulk = 1.72e-6 #ohms-cm, resistivity of bulk copper
    Sconv = 944822 #pixels/cm, size conversion from pixels to cm
    bulkR = 8960 #kg/m^3, density of bulk copper
    
    colbd = 'black'
    colfc = 'red'
    
    xsize = 52
    ysize = 192
    zsize = 192
    noparts = 97
    
    lim = 0.00001#0.1
    etaCut = 0.00001#0.0000001
    
    Alist = list(range(0,12,2)) + list(range(20,210,10)) + list(range(400,1200,200)) + list(range(1250,3250,250)) + [3150,3300]
    
    AIlist = []
    for ik in range(m.ceil(len(Alist)/sizep)):
        aid = num+ik*sizep
        if aid < len(Alist):
            AIlist.append(Alist[aid])
    
    for ai in AIlist:
    
        #ai = AIlist[num] #int(argv[1]) #10
        SNo = 504
    
        #analysis box properties
        xlb = 47
        ylb = 47
        boxSize = 95
        xhb = xlb + boxSize
        yhb = ylb + boxSize
    
        propsDic = {}
        nodectr = 0
        nodeDic = {}
        nodeList = [[],[],[]]
    
        starttime = time.time()
    
        fileName = "fullT"+str(ai)+"SN"+str(SNo)+".dat" 
    
        fileCheck = "fullT"+str(ai)+"Out.log"
        fileBool = isfile('./'+fileCheck)
    
        if fileBool:
            neckInfoDic = getRead(fileCheck)
        else:
            #direct output to file
            print_file = open(fileCheck, 'w')
            sys.stdout = print_file
    
        bnf = open(fileName,"rb")
        fcont = bnf.read()
        bnf.close()
    
        fullsize = int(len(fcont)/8)
        Vdfull = struct.unpack('d'*fullsize,fcont[:fullsize*8]) #contains float information for the data file
    
        size = xsize*ysize*zsize
    
        Xd = np.zeros(size)
        Yd = np.zeros(size)
        Zd = np.zeros(size)
        Vd = np.zeros(size)
        Ed = np.zeros([noparts,size])
    
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
                        Ed[p][counti] = Vdfull[countv+p+1]
                        
                    counti = counti+1
    
        x = [[],[]]; y = [[],[]]; z = [[],[]]; v = [[],[]]
    
        Vt = Vd>lim
    
        midtime = time.time()
        print('Time to initialize:')
        print(midtime-starttime,'s')
    
        rhoAnalysis(ai) 
        
        print('Time to do Density Analysis:')
        print(time.time()-midtime,'s')
    
        #%%
        midtime = time.time()
    
        plotEta = False
        plotCheckSpheres = False
        overlapAnalysis = True
        plotOvlp = False
        neckPlot = False
    
        for p1 in range(noparts):
            t1 = Ed[p1] >= etaCut
            idd = Vt&t1
            x[0] = Xd[idd]
            y[0] = Yd[idd]
            z[0] = Zd[idd]
            v[0] = Vd[idd]
            
            if v[0].any():
                sp = str(p1)
                if sp in propsDic:
                    rad0, cen0 = propsDic[sp][0]
                else:
                    #print("particle",p1,end=" ==> ")
                    rad0, cen0 = props(x[0],y[0],z[0])#,True) 
                    propsDic[sp]=[[rad0,cen0],[min(x[0]),max(x[0]),min(y[0]),max(y[0])]]
                
                if plotCheckSpheres:
                    pSphere = makeSphere(cen0[0],cen0[1],cen0[2],rad0)
                    fig = plt.figure()
                    ax = Axes3D(fig)
            
                    ax.set_xlim([0,ysize])
                    ax.set_ylim([0,ysize])
                    ax.set_zlim([0,ysize])
            
                    for ii in range(1):
                        ax.scatter(x[ii],y[ii],z[ii],s=9.8,c= colfc,marker='o',edgecolor = colbd)
                    
                    ax.scatter(pSphere[0],pSphere[1],pSphere[2],s=9.8,c= 'y',marker='o',edgecolor = colbd)
                    
                    ax.view_init(elev=40,azim=50)
            
                    titleName = "Sphere And Particle " + sp
                    ax.set_title(titleName)
            
                    print(p1,"done")
                    fig.savefig("furthestPixRadius_Sphere_cutOff%s_part%d.png" %(str(etaCut),p1))
                    plt.close() 
                                
                if overlapAnalysis:
                    
                    x[1] = []; y[1] = []; z[1] = []; v[1] = []
                    
                    for p2 in range(p1+1,noparts):
                        t1 = Ed[p2] >= etaCut
                        idd = Vt&t1
                        x[1] = Xd[idd]
                        y[1] = Yd[idd]
                        z[1] = Zd[idd]
                        v[1] = Vd[idd]
                
                        if v[1].any():
                            sp = str(p2)
                            if sp in propsDic:
                                rad1, cen1 = propsDic[sp][0]
                            else:
                                #print("particle",p2,end=" ==> ")
                                rad1, cen1 = props(x[1],y[1],z[1])#,True)
                                propsDic[sp]=[[rad1,cen1],[min(x[1]),max(x[1]),min(y[1]),max(y[1])]]
                            
                            if fileBool:
                                sp1 = str(p1); sp2 = str(p2)
                                keyN = 'p'+sp1+'p'+sp2
                                if keyN in neckInfoDic:
                                    radN,lenN = neckInfoDic[keyN]
                                    cenD = np.sqrt(sum((cen0 - cen1)**2))
                                    if sp1 not in nodeDic:
                                        nodeDic[sp1] = nodectr
                                        nodectr += 1
                                    if sp2 not in nodeDic:
                                        nodeDic[sp2] = nodectr
                                        nodectr += 1
                                    nodeList[0].append(nodeDic[sp1]) 
                                    nodeList[1].append(nodeDic[sp2])
                                    nodeList[2].append(calcRes(rad0,rad1,radN,lenN))
                            else:
                                titleName, ovlpA = inContact(x,y,z,v,[rad0,rad1],[cen0,cen1],False,False,p1=p1,p2=p2)
                                if titleName:
                                    cenD = np.sqrt(sum((cen0 - cen1)**2))
                
                                    #unpack ovlp
                                    #ovlp,radN,cenN,finPos,overlap01 = ovlpA
                                    ovlp, radN, overlap01, lenN, neckAll = ovlpA
                                    
                                    sp1 = str(p1); sp2 = str(p2)
                                    if sp1 not in nodeDic:
                                        nodeDic[sp1] = nodectr
                                        nodectr += 1
                                    if sp2 not in nodeDic:
                                        nodeDic[sp2] = nodectr
                                        nodectr += 1
                                    nodeList[0].append(nodeDic[sp1]) 
                                    nodeList[1].append(nodeDic[sp2])
                                    nodeList[2].append(calcRes(rad0,rad1,radN,lenN))
                                    
                                    if neckPlot:
                                        plotNeck(ovlp,'b4Narrowed'+str(p1)+'_and_'+str(p2)+'.png')
                                        #plotNeck(finPos,'lim'+str(etaCut)+'_NarrowedNeck'+str(p1)+'_and_'+str(p2)+'.png')
                                        plotNeck(neckAll,'lim'+str(etaCut)+'_NarrowedNeck'+str(p1)+'_and_'+str(p2)+'.png')
                                    
                                if plotEta:
                                    fig = plt.figure()
                                    ax = Axes3D(fig)
                            
                                    ax.set_xlim([0,ysize])
                                    ax.set_ylim([0,ysize])
                                    ax.set_zlim([0,ysize])
                            
                                    for ii in range(2):
                                        ax.scatter(x[ii],y[ii],z[ii],s=100*v[ii],c= colfc,marker='o',edgecolor = colbd)
                                    
                                    ax.view_init(elev=40,azim=50)
                            
                                    ax.set_title(titleName)
                        
                                    print(p1,"+",p2)
                                    fig.savefig("Part%d_and_%d.png" %(p1,p2))
                                    plt.close()  
                                
                                if plotOvlp:
                                    fig = plt.figure()
                                    ax = Axes3D(fig)
                            
                                    ax.set_xlim([0,ysize])
                                    ax.set_ylim([0,ysize])
                                    ax.set_zlim([0,ysize])
                                    
                                    ovlp = sortPlot(ovlp)
                                    ax.scatter(ovlp[0],ovlp[1],ovlp[2],s=100*np.array(ovlp[3]),c= colfc,marker='o',edgecolor = colbd)
                                    
                                    ax.view_init(elev=40,azim=50)
                                    titleName = "Overlap "+ titleName
                                    ax.set_title(titleName)
                        
                                    print(p1,"+",p2,"Overlap")
                                    fig.savefig("2LROverlapPart%d_and_%d.png" %(p1,p2))
                                    plt.close()                     
    
        minX = [0,100]; minY = [0,100]
        maxX = [0,0]; maxY = [0,0]
        listX = [[],[]]
        listY = [[],[]]
        iilist = []
    
        for ii,ij in propsDic.items():
            if ii in nodeDic:
                #using mins and maxs
                xmin,xmax,ymin,ymax = ij[1]
                listX[0].append(xmin)
                listX[1].append(xmax)
                listY[0].append(ymin)
                listY[1].append(ymax)
                iilist.append(ii)
                '''
                xval,yval = ij[0][1][0],ij[0][1][1]
                xmin,xmax,ymin,ymax = xval,xval,yval,yval
                if xmin < minX[1]: minX[0] = ii; minX[1] = xmin
                if ymin < minY[1]: minY[0] = ii; minY[1] = ymin    
                if xmax > maxX[1]: maxX[0] = ii; maxX[1] = xmax
                if ymax > maxY[1]: maxY[0] = ii; maxY[1] = ymax 
                '''
                
        xminA = getOpt(listX[0],0)
        xmaxA = getOpt(listX[1],1)
        yminA = getOpt(listY[0],0)
        ymaxA = getOpt(listY[1],1) 
    
        xbdsN = [dicExtract(xminA),dicExtract(xmaxA)]#[nodeDic[minX[0]],nodeDic[maxX[0]]]
        ybdsN = [dicExtract(yminA),dicExtract(ymaxA)]#[nodeDic[minY[0]],nodeDic[maxY[0]]]
    
        if nodeList[0]:
            writeNodeList(nodeList,'nodeFileT'+str(ai)+'.txt',xbdsN,ybdsN,vfileName = 'VoltageFileXY_T'+str(ai)+'.txt')
            maxV = max(max(nodeList[0]),max(nodeList[1]))
            writeVolList(maxV,'VoltageFileT'+str(ai)+'.txt')
            
        print('Time to end:')
        print(time.time()-midtime,'s')     
    
        if not fileBool:
            sys.stdout = sys.__stdout__
            print_file.close()
    
    return

jobs = []
for i in range(sizep):
        p = mp.Process(target=worker, args=(i,))
        jobs.append(p)
        p.start()
