import numpy as np
import math
from pylab import *
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def actd(t,temp):                       

        if temp == 450:
                T = [0.1572,1.2099,0.2854,0.8681,-0.0053]

        elif temp == 500:
                T = [0.2354,0.5326,0.8508,1.212,-0.4228]

        elif temp == 550:
                T = [0.2543,0.2889,0.5528,1.2493,-0.5253]
                
        elif temp == 600:
                T = [0.1651,0.4672,0.1331,0.8051,-0.011]

        valu = T[0]*exp((-T[1]/(t+T[2]))+T[3])+T[4]

        return valu

def actdTest(t,T):                      

        valu = T[0]*exp((-T[1]/(t+T[2]))+T[3])+T[4]

        return valu

def actdFit(t,T0,T1,T2,T3,T4):                  

        valu = T0*exp((-T1/(t+T2))+T3)+T4

        return valu

def minID(Col):
        valm = min(Col)
        i = 0
        while i < len(Col):
                if (Col[i] == valm):
                        vali = i
                        i = len(Col)
                else:
                        i = i + 1

        return (valm,vali)

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
                                        outf.write('%f' %kk)
                                        outf.write('\t')
                                        outf.write('\n')

        outf.close()

def inbetT(simT,simR,lowR,highR,calB):
                
        ctIn = 0
        ttA = []
        
        for tt,rr in zip(simT,simR):
                tc = float(tt)/calB
                lowB = actdTest(tc,lowR)
                highB = actdTest(tc,highR)
                t = (rr - lowB)/(highB - lowB)
                ttA.append(t)
                if t>=0 and t<=1:
                        ctIn+=1
                        
        ctOut = len(simR)-ctIn

        return ctOut,mean(ttA)

def COV(Xmin,Xmax,step):

        covv = 0
        for mn1,mx1 in zip(Xmin,Xmax):
                av1 = (mn1+mx1)/2
                for mn2,mx2 in zip(Xmin,Xmax):
                        av2 = (mn2+mx2)/2
                        summ = 0
                        ctr = 0
                        for ij in range(mn1,mx1+step,step):  
                                for iv in range(mn2,mx2+step,step):
                                        summ = summ+iv*ij
                                        ctr += 1
                                      
                        covv = covv+(summ/ctr-av1*av2)
        return covv

dirstd = str(os.getcwd())

TEMP = 550

if TEMP == 450:
        Ttime = [1,2,4,5,7.5,10,30,45]
        Td = [0.1572,1.2099,0.2854,0.8681,-0.0053]#[0.296129599858139, 0.573749703904609, 1.48273030967169, 1.42529100705374, -0.835644039167063]#
        ACof = [[0.1,0.5],[0.8,1.5],[0.1,0.5],[0.4,1.1],[-1,1]]
        Tvals = [0.141913443,0.217869154,0.267871757,0.27611845,0.316230947,0.358823953,0.363547613,0.384753104]
        TvalsE = [0.051631389,0.065867917,0.05344458,0.069630836,0.05056041,0.05242889,0.049289925,0.040739769]
        TvalsE_2 = [0.025815695,0.032933958,0.02672229,0.034815418,0.025280205,0.026214445,0.024644962,0.020369885]
elif TEMP == 500:
        Td = [0.2354, 0.5326, 0.8508, 1.212, -0.4228]#[0.348110047752255, 0.341680160776424, 1.2495102539606, 1.53282511992199, -1.22314729511812]#
        ACof = [[0.1,0.5],[0.3,1],[0.6,1.1],[0.8,1.5],[-1,1]]
        Ttime = [1,2,4,5,7.5,10,30,45]
        Tvals = [0.171896686705849,0.236480116951352,0.277864355419302,0.280633538625766,0.325941559448155,0.351375314372826,0.378378776903819,0.378169617402741]
        TvalsE = [0.0514792784973066,0.0475241693567705,0.052206541720242,0.0637139964013396,0.0417327581506577,0.0504100783455164,0.0426071738457872,0.0507123355883403]
        TvalsE_2 = [0.0257396392486533,0.0237620846783852,0.026103270860121,0.0318569982006698,0.0208663790753289,0.0252050391727582,0.0213035869228936,0.0253561677941702]
elif TEMP == 550:
        Td = [0.2543, 0.2889, 0.5528, 1.2493, -0.5253]#[0.260674054525433, 0.199837362299377, 0.726292975704178, 1.78007782323173, -1.17221051030361]#
        ACof = [[0.1,0.5],[0.05,0.7],[0.4,0.9],[0.8,1.5],[-1,1]]
        Ttime = [1,2,4,5,7.5,10,30,45]  
        Tvals = [0.212421170165928,0.279798362337399,0.283124724601518,0.297578737506606,0.352488688167634,0.341003082523426,0.370484097989162,0.376609403612535]       
        TvalsE = [0.0502117901605388,0.0526405201443995,0.0512600203744162,0.0499200468124537,0.0459796952315542,0.0405759465768843,0.0489677145014584,0.0353044414735159]
        TvalsE_2 = [0.0251058950802694,0.0263202600721998,0.0256300101872081,0.0249600234062269,0.0229898476157771,0.0202879732884421,0.0244838572507292,0.017652220736758]
elif TEMP == 600:
        Td = [0.1651, 0.4672, 0.1331, 0.8051, -0.011]#[0.268332715883733, 0.235242854887199, 0.479348492096105, 1.28109278342303, -0.59129594400154]#
        ACof = [[0.1,0.5],[0.05,0.75],[0.1,0.5],[0.4,1.1],[-1,1]]
        Ttime = [1,2,5,7.5,15,45]
        Tvals = [0.233405784023887,0.289472196525101,0.332755918847125,0.339088034575323,0.363536836588568,0.373289944549384]
        TvalsE = [0.0398590170555035,0.0421532471807141,0.0453847860437572,0.0455396341885881,0.038651636783154,0.0428095671582943]
        TvalsE_2 = [0.0199295085277518,0.0210766235903571,0.0226923930218786,0.0227698170942941,0.019325818391577,0.0214047835791472]


finalT = 330 #final time that simulation ran to
timestep = 50 #time steps used in file i.e. distance between files (filed in python codes)
noT = int((finalT*100/timestep) + 1)
boxSize = 95


ff0 = []
fileList = ["AllCoeficients_0_neg1.txt","AllCoeficients_0_pos1.txt","AllCoeficients_0_neg1.5.txt","AllCoeficients_0_pos1.5.txt"]
ctr = 0
for fileN in fileList:
        with open(fileN) as fil:
                fll = fil.readlines()

        
        for ff in fll:
                ff = ff.split('\t')
                ff0.append([])
                for ii in ff:
                        ff0[ctr].append(float(ii))
                ctr+=1        


ffno0 = []

fileList = ["AllCoeficients_neg1.txt","AllCoeficients_pos1.txt","AllCoeficients_neg1.5.txt","AllCoeficients_pos1.5.txt"]
ctr = 0
for fileN in fileList:
        with open(fileN) as fil:
                fll = fil.readlines()

        for ff in fll:
                ff = ff.split('\t')
                ffno0.append([])
                for ii in ff:
                        ffno0[ctr].append(float(ii))
                ctr+=1        

INS="OnlyZero"

dataEndtime = 3.1
Arrl = 2

bandS = 1

if abs(bandS-1)<1e-15:
        SF = 0
if abs(bandS-1.5)<1e-15:
        SF = 8
        
ffUse = ff0[:]
SNo = 50400002

idd = (TEMP-450)/50 + SF

startlow = 10000
starthigh = 100000
step = 100

filetord = "SN%d_48.txt" %SNo

bedList = ['bed1','bed3','bed5','bed6','bed11','bed14','bed15','bed16','bed19','Newbed0','Newbed1','Newbed4','Newbed5']

filAll = open("CollatedData_Band%sData_%s.txt" %(str(bandS),INS),'w')

for bedn in bedList:
        
        newdr = dirstd+'/'+bedn

        if bedn == 'bed6':
                xbd = [55,60,62,63,65,70]
                ybd = [55,60,62,63,65,70]              
        else:
                xbd = [40,45,47,48,50,55]
                ybd = [40,45,47,48,50,55]

        '''
        if bedn[-3]=='5':
                SNo=int(bedn[-3:]+'00002')
        else:
                SNo=int(bedn[-4:]) #newer beds
        '''
        
        os.chdir(newdr)

        filAll.write("%s\n" %bedn)

        filBed = open("%s_Band%sData_SN%d_T%dC_%s.txt" %(INS,str(bandS),SNo,TEMP,bedn), 'w')
        filBed.write("xlb\tylb\taverage\tstdev\t#pointsOut\t#CalibFactors\tminCalib\tmaxCalib\n")
        filAll.write("xlb\tylb\taverage\tstdev\t#pointsOut\t#CalibFactors\tminCalib\tmaxCalib\n")
        
        for SNoctr in range(len(xbd)):
                
                xlb = xbd[SNoctr]
                ylb = ybd[SNoctr]
                
                filetord = "SN%ddatafull2by2square_xl%d_xh%d_lim_0.txt" %(SNo,xlb,xlb+boxSize)
                fb = open(filetord,"r")
                lines = fb.readlines()
                fb.close()

                fres = []
                frst = []

                simT = [0 for x in range(noT)]
                simR = [0 for x in range(noT)]

                cttr = 0

                for i in lines:
                        if cttr <= noT:
                                fres.append(i.split('\t')[9])
                                frst.append(i.split('\t')[0])
                                cttr = cttr + 1
                        
                inR = float(fres[1])

                for i in range(1,noT+1):
                        fullsimd = float(fres[i])
                        simR[i-1] = float(fullsimd - inR)/float(fullsimd)
                        simT[i-1] = int(frst[i])

                clib = []
                smallT = 1000
                ctsS = 0

                for calB in range(startlow,starthigh,step):
                        lowR = ffUse[idd][:]
                        highR = ffUse[idd+4][:]
                        CtOut = inbetT(simT,simR,lowR,highR,calB)
                        clib.append(CtOut[0])
                        
                        if abs(CtOut[0]-smallT)<1e-15:
                                ctsS+=1
                        else:
                                if CtOut[0]<smallT:
                                        ctsS=0
                                
                        smallT = min(CtOut[0],smallT)
                        #print calB,' ==> ',float(CtOut[0])/len(simR),CtOut[0],len(simR)

                if smallT==0:
                        alow = startlow+minID(clib)[1]*step
                        bhigh = startlow+(minID(clib)[1]+ctsS)*step

                        filBed.write("%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\n" %(xlb,ylb,(alow+bhigh)/2,(bhigh-alow)/sqrt(12),smallT,ctsS,alow,bhigh))
                        filAll.write("%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\n" %(xlb,ylb,(alow+bhigh)/2,(bhigh-alow)/sqrt(12),smallT,ctsS,alow,bhigh))
                        
                        fig = plt.figure()

                        array = []
                        tarr = []
                        for ii in np.arange(0,dataEndtime,0.1):
                                array.append(actdTest(ii,Td))
                                tarr.append(ii)
                        plt.plot(tarr,array,'b')
                                
                        for io in range(2):
                                Tdtemp = ffUse[idd+4*io][:]
                                array = []
                                tarr = []
                                for ii in np.arange(0,dataEndtime,0.1):
                                        array.append(actdTest(ii,Tdtemp))
                                        tarr.append(ii)
                                plt.plot(tarr,array,'r')
                                
                        plt.errorbar(Ttime[0:Arrl], Tvals[0:Arrl], yerr=TvalsE[0:Arrl], fmt='.k')
                        plt.errorbar(Ttime[0:Arrl], Tvals[0:Arrl], yerr=TvalsE_2[0:Arrl], fmt='.y')

                        

                        for cal in range(alow,bhigh+step,step):
                                tarr = []
                                for it in simT:
                                        tarr.append(float(it)/cal)
                                plt.plot(tarr,simR,'g--')
                                
                        fig.savefig("%s_Band%sData_SN%d_T%dC_%s_box%d.png" %(INS,str(bandS),SNo,TEMP,bedn,xlb))
                        plt.close()

                print bedn," ---Done--- ",xlb
                
        filBed.close()
        
        os.chdir(dirstd)
        
filAll.close()
