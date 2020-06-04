from statistics import stdev,mean
import struct
import math as m
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from numpy import random
import os
from shutil import copyfile

def toFile(Pd,fname='outPosition.txt',filout='int',tag = []):

    outf = open(fname,'w')
    
    for ij in Pd:
        ct = 0
        for kk in ij:
            if ct < len(tag):
                outf.write(tag[ct])
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
                else:
                    if filout=='float':
                       outf.write('%f' %kk) 
            outf.write('\t')
            ct+=1
        outf.write('\n')

    outf.close()

def plotParticle(Particle,figName='BeautPartplot.png',colFace='white',colEdge_part='purple',ysize=104,elevN=30,azimA=45):
    
    sphereList = makeSphere(Particle[0],Particle[1],Particle[2],Particle[3])
    
    fig = plt.figure()
    ax = Axes3D(fig)
                    
    ax.set_xlim([0,ysize])
    ax.set_ylim([0,ysize])
    ax.set_zlim([0,ysize])
    
    ax.scatter(sphereList[0],sphereList[1],sphereList[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge_part)
    
    ax.view_init(elev=elevN,azim=azimA)
              
    fig.savefig(figName)
    plt.close() 

def plotSurf(Surf,figName='Beautplot.png',colFace='white',colEdge='green',ysize=104,elevN=30,azimA=45):
    
    fig = plt.figure()
    ax = Axes3D(fig)
                    
    ax.set_xlim([0,ysize])
    ax.set_ylim([0,ysize])
    ax.set_zlim([0,ysize])
    
    ax.scatter(Surf[0],Surf[1],Surf[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge)

    ax.view_init(elev=elevN,azim=azimA)
              
    fig.savefig(figName)
    plt.close()

def plotParticleSurf(Particle,Surf,figName='BeautPartSurfplot.png',colFace='white',colEdge_surf='green',colEdge_part='purple',ysize=104,elevN=30,azimA=45):
    
    sphereList = makeSphere(Particle[0],Particle[1],Particle[2],Particle[3])
    
    fig = plt.figure()
    ax = Axes3D(fig)
                    
    ax.set_xlim([0,ysize])
    ax.set_ylim([0,ysize])
    ax.set_zlim([0,ysize])
    
    ax.scatter(Surf[0],Surf[1],Surf[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge_surf)
    ax.scatter(sphereList[0],sphereList[1],sphereList[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge_part)
    
    ax.view_init(elev=elevN,azim=azimA)
              
    fig.savefig(figName)
    plt.close()    
    
def plotParticlesSurf(pinA,Surf,figName='BeautPartsSurfplot.png',colFace='white',colEdge_surf='green',colEdge_part='purple',ysize=104,elevN=30,azimA=45):
    
    sphereList = [[],[],[],[]]
    for x,y,z,r in zip(pinA[0],pinA[1],pinA[2],pinA[3]):
        partS = makeSphere(x,y,z,r)
        sphereList[0] += partS[0]
        sphereList[1] += partS[1] 
        sphereList[2] += partS[2] 
    
    fig = plt.figure()
    ax = Axes3D(fig)
                    
    ax.set_xlim([0,ysize])
    ax.set_ylim([0,ysize])
    ax.set_zlim([0,ysize])
    
    ax.scatter(Surf[0],Surf[1],Surf[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge_surf)
    ax.scatter(sphereList[0],sphereList[1],sphereList[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge_part)
    
    ax.view_init(elev=elevN,azim=azimA)
              
    fig.savefig(figName)
    plt.close() 

def plotParticles(pinA,figName='BeautPartsplot.png',colFace='white',colEdge_part='purple',ysize=104,elevN=30,azimA=45):
    
    sphereList = [[],[],[],[]]
    for x,y,z,r in zip(pinA[0],pinA[1],pinA[2],pinA[3]):
        partS = makeSphere(x,y,z,r)
        sphereList[0] += partS[0]
        sphereList[1] += partS[1] 
        sphereList[2] += partS[2] 
    
    fig = plt.figure()
    ax = Axes3D(fig)
                    
    ax.set_xlim([0,ysize])
    ax.set_ylim([0,ysize])
    ax.set_zlim([0,ysize])
    
    ax.scatter(sphereList[0],sphereList[1],sphereList[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge_part)
    
    ax.view_init(elev=elevN,azim=azimA)
              
    fig.savefig(figName)
    plt.close() 

def plotSpheres(sphereList,figName='BeautPartsplot.png',colFace='white',colEdge_part='purple',ysize=104,elevN=30,azimA=45):
        
    fig = plt.figure()
    ax = Axes3D(fig)
                    
    ax.set_xlim([0,ysize])
    ax.set_ylim([0,ysize])
    ax.set_zlim([0,ysize])
    
    ax.scatter(sphereList[0],sphereList[1],sphereList[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge_part)
    
    ax.view_init(elev=elevN,azim=azimA)
              
    fig.savefig(figName)
    plt.close() 

def plotSpheresSurf(sphereList,Surf,figName='BeautPartsSurfplot.png',colFace='white',colEdge_surf='green',colEdge_part='purple',ysize=104,elevN=30,azimA=45):

    fig = plt.figure()
    ax = Axes3D(fig)
                    
    ax.set_xlim([0,ysize])
    ax.set_ylim([0,ysize])
    ax.set_zlim([0,ysize])
    
    ax.scatter(Surf[0],Surf[1],Surf[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge_surf)
    ax.scatter(sphereList[0],sphereList[1],sphereList[2],s = 9.8,c= colFace,marker='o',edgecolor = colEdge_part)
    
    ax.view_init(elev=elevN,azim=azimA)
              
    fig.savefig(figName)
    plt.close()
    
def animaUpdate(num):
    
    num-=1
    
    if num==-1:
        surfPlot = resortSurf(surfTemp,1) 
        scat._offsets3d = (surfPlot[0],surfPlot[1],surfPlot[2])  
    else:
        xc = pinA[0][num]; yc = pinA[1][num]; zc = pinA[2][num]; rr = pinA[3][num]
        castShadow(surfTemp,[xc,yc,zc,rr])
        print('Done particle:', num, 'size of list now -->', len(surfTemp),file=logOut) 
        surfPlot = resortSurf(surfTemp,1) 
        scat._offsets3d = (surfPlot[0],surfPlot[1],surfPlot[2])  
        
def positionVariation(*fileno):
    """
    This function calculates the amount of variation between the particles
    in file numbers provided. This provides numerical quantification to determine 
    if steady state has been reached
    
    Parameters
    ----------
    *fileno : variadic list
        input the ID number for the list of files to iterate over to determine 
        if steady state has been reached

    Returns
    -------
    tuple of standard deviations for the x, y and z centers over the range of 
    files passed        
    """

    Xo = []; Yo = []; Zo = [];
    aaN = 0
    
    for aa in fileno:
        filenam = 'particleOut%d.txt' %aa       
        
        nop = 0
        
        with open(filenam,'r') as fild:
            for line in fild.readlines():
                cols = line.split('\t')

                if aaN==0:
                    Xo.append([])
                    Yo.append([])
                    Zo.append([])
                
                Xo[nop].append(float(cols[0]))
                Yo[nop].append(float(cols[1]))
                Zo[nop].append(float(cols[2]))
                nop+=1
       
        aaN+=1

    XoV = [stdev(ix) for ix in Xo]
    YoV = [stdev(ix) for ix in Yo]
    ZoV = [stdev(ix) for ix in Zo]

    return mean(XoV),mean(YoV),mean(ZoV)

def findOldHeight(surf):
    """
    Parameters
    ----------
    surf: list of the bottom layer surface cloud points

    Returns
    -------
    the maximum height of old layer
    the average height of old layer
    """

    heightDict = {}
    
    for x,y,z in zip(surf[0],surf[1],surf[2]):
        key = str(x)+','+str(y)
        if key in heightDict.keys():
            heightDict[key]=max(z,heightDict[key])
        else:
            heightDict[key] = z
    
    oldZmax = max(heightDict.values())
    oldZhighAvg = mean(heightDict.values())
    
    return oldZmax, oldZhighAvg

def makeParticles(pinA):
    
    sphereList = [[],[],[]]
    for x,y,z,r in zip(pinA[0],pinA[1],pinA[2],pinA[3]):
        partS = makeSphere(x,y,z,r)
        sphereList[0] += partS[0]
        sphereList[1] += partS[1] 
        sphereList[2] += partS[2] 
    
    return sphereList

def Overlap(pinA,xx = [],partpart=True,partSurf=True):
    """
    This function determines the amount of overlap in pixels between the particles
    of new layers as well as particles of new layer and the bottom layer    

    Parameters
    ----------
    pinA : nested list of new layer particle information
    xx : nested list of full bottom layer information
    partpart : flag to determine if particle-particle overlap should be calculated
        DESCRIPTION. The default is True.
    partSurf : flag to determine if particle-bottom layer overlap should be calculated
        DESCRIPTION. The default is True.

    Returns
    -------
    PartsOverlap : list
        contains:    
            particle #, list--> particle-particleOverlappings #, sum particle - particleoverlap
            particle surface overlap # of pixels overlapping, %of overlapping(wtr to full particle)

    """
    noParts = len(pinA[0])
    PartsOverlap = [] 
    xxlist = [str(ii) for ii in xx]
    
    for pii in range(noParts):
        pin = [pinA[0][pii], pinA[1][pii], pinA[2][pii], pinA[3][pii]]
        
        if partpart:    
            
            ovlpn = []
            overlap = 0
            
            #checking for overlap between particle pii and the other particles in 
            #the new layer
            for bi in range(noParts):
                if bi != pii:
                    xdif = pin[0] - pinA[0][bi]
                    ydif = pin[1] - pinA[1][bi]
                    zdif = pin[2] - pinA[2][bi]
                    
                    ovl = pin[3]+pinA[3][bi] - m.sqrt(xdif**2+ydif**2+zdif**2)
                    if ovl>=0:
                        ovlpn.append(bi)
                        overlap+=ovl
            
            PartsOverlap.append([pii,ovlpn,overlap])
   
        if partSurf:
            newxx = []
            
            pinB = [[m.floor(pin[0]-pin[3]),m.ceil(pin[0]+pin[3]+1)],[m.floor(pin[1]-pin[3]),m.ceil(pin[1]+pin[3]+1)],[m.floor(pin[2]-pin[3]),m.ceil(pin[2]+pin[3]+1)]]
            for x in range(*pinB[0]):
                for y in range(*pinB[1]):
                   for z in range(*pinB[2]):
                       if ((x-pin[0])**2+(y-pin[1])**2+(z-pin[2])**2)<=pin[3]*pin[3]:
                           newxx.append(str([x,y,z]))

            xxcomb = set(xxlist)&set(newxx)
            PartsOverlap[pii].append(len(xxcomb))
            PartsOverlap[pii].append(len(xxcomb)/len(newxx))
        
    return PartsOverlap

def analyzeOverlapInfo(overlapList,overlapLimit = 0.2):
    '''
    Calculates the number of particles that are in lower layer, where in lower
    layer is defined by having more overlapping % of pixels than overlapLimit

    Parameters
    ----------
    overlapList : list containing overlap info
    overlapLimit : double, optional
        cutoff used to determine if particle is in bottom layer. The default is 0.2.

    Returns
    -------
    ctr : int
        Number of particles from top layer that are in bottom layer
    '''
    
    ctr = 0
    for ii in overlapList:
        if ii[-1]>=overlapLimit:
            ctr+=1
            
    return ctr

def resortSurf(oldList,sortType = 0):
    
    if sortType:
        emptyList = [[],[],[]]
        for ii in oldList:
            emptyList[0].append(ii[0][0])
            emptyList[1].append(ii[0][1])
            emptyList[2].append(ii[1])
    else:
        emptyList = []
        for x,y,z in zip(oldList[0],oldList[1],oldList[2]):
            emptyList.append([[x,y],z])   
        
    return emptyList

def makeCircle(xC,yC,R):
    
    lowb = [m.floor(xC-R),m.floor(yC-R)]
    diam = m.ceil(2*R)
    
    circPoints = []
    
    #print ('Circle properties-->',xC,yC,R)
    for ix in range(lowb[0],lowb[0]+diam+3):
        for iy in range(lowb[1],lowb[1]+diam+3):
            if (ix-xC)**2+(iy-yC)**2<=R**2:
                circPoints.append([ix,iy])
    
    return circPoints

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

def castShadow(fullList,particle):
    '''
    FullList changed in this function with the points in the particles
    shadow removed

    Parameters
    ----------
    fullList : List of points in the surface point cloud with the list in the form
                [[x,y],z,v] so the modified/copied list with the points to be removed
    particle : List of particle [x_center,y_center,radius]

    Returns
    -------
    None
    '''
    
    #make list with pixels in the circle
    circleList = makeCircle(particle[0],particle[1],particle[3])
    
    zlow = particle[2]
    #print('Particle Info: ',*particle)
    #print('Cicle List:')
    #print(*circleList,sep = '\n')
    
    #Make copy for iterations
    listCopy = fullList.copy()
    
    for jk in circleList:
        for ii in listCopy:
            if jk in ii and ii[1] <= zlow:
                #print(jk,'--in-->',ii)
                fullList.remove(ii)
    
def extractTop(fullList):
    '''
    Function searches for the top layer of the left over surface and deletes all
    points below
    Parameters
    ----------
    fullList : list of layer

    Returns
    -------
    None.
    '''

    listCopy = fullList.copy()
    heightDict = {}
    fullList.clear()
    
    for [x,y],z in listCopy:

        key = str(x)+','+str(y)
        if key in heightDict.keys():
            heightDict[key].append(z)
        else:
            heightDict[key] = [z]
     
    for key,val in heightDict.items():
        ikey = key.split(',')
        fullList.append([[int(ikey[0]),int(ikey[1])],max(val)])
    
    '''
    for jk in listCopy:
        ij = jk[0]
        inPart = [kk for kk in fullList if ij in kk]
        #print(inPart,end='--->')
        minR = max(inPart,key = lambda x: x[0][1]+x[0][0]+x[1])
        inPart.remove(minR)
        #print(minR,'--->',inPart)
        for kk in inPart:
            fullList.remove(kk) 
    '''
    
def cluster(fullList,plotClusters = True):
    '''
    This function takes a point cloud and clusters it according to the pixel
    closeness values set below:
        pixDx -> proximity of x values to consider a cluster (1 pixel default)
        pixDy -> proximity of y values to consider a cluster (1 pixel default)
        pixDz -> proximity of z values to consider a cluster (3 pixels default)

    Parameters
    ----------
    fullList : List containing point cloud 
        The function removes the clustered points found from the list passed to it
    plotClusters : boolean, optional
        True or False boolean to state if the Clusters should be plotted. 
        The default is True.

    Returns
    -------
    clusters : list of class clussAtr
        The list returned contains a dummy class called clussAtr which has the 
        following cluster attributes:
            xSpread -> list of min and max x values in that cluster
            ySpread -> list of min and max y values in that cluster
            zSpread -> list of min and max z values in that cluster
            xyStart -> list of the x and y starting pixel for that cluster 
                        (basically gives one identifying pixel location for that cluster)
            noPixels -> number of pixels captured in the cluster
    '''
    
    pixDx = 1; pixDy = 1; pixDz = 3
    clusters = []
    ctp = 0
    ctt = 0
    
    class clussAttr:
        #Dummy class to make all the cluster attributes easier to unpack
        def __init__(self,ii,ctp,plotS):
            self.xSpread = [min(plotS[0]),max(plotS[0])]
            self.ySpread = [min(plotS[1]),max(plotS[1])]
            self.zSpread = [min(plotS[2]),max(plotS[2])]
            self.xyStart = [ii[0][0],ii[0][1]]
            self.noPixels = ctp
            
        def __str__(self):
            string = 'xyStart: ['+str(self.xyStart[0])+','+str(self.xyStart[1])+']'
            string += '\tnoPixels: '+str(self.noPixels)
            string += '\txSpread: ['+str(self.xSpread[0])+','+str(self.xSpread[1])+']'
            string += '\tySpread: ['+str(self.ySpread[0])+','+str(self.ySpread[1])+']'
            string += '\tzSpread: ['+str(self.zSpread[0])+','+str(self.zSpread[1])+']\n'
            return string
    
    if plotClusters:
        fig = plt.figure()
        ax = Axes3D(fig)
                        
        ax.set_xlim([0,ysize])
        ax.set_ylim([0,ysize])
        ax.set_zlim([0,ysize])
        
        colEdge_surf = tuple(random.rand(3)) #col[ctt%len(col)]
        
        ax.scatter([],[],[],s = 9.8,c= 'white',marker='o',edgecolor = colEdge_surf)

        ax.view_init(elev=elevN,azim=azimA)
        
        ctt+= 1          
            
    while len(fullList):
        
        ii = fullList[0]
        
        ctp = 0
        Exit = True
        
        while Exit:
            ij = fullList[0]
            if ctp == 0: 
                plotS = [[ij[0][0]],[ij[0][1]],[ij[1]]] 
                fullList.remove(ij)
            
            for enx,eny,enz in zip(plotS[0],plotS[1],plotS[2]):
                listCopy = [kk for kk in fullList if (abs(enx - kk[0][0]) <= pixDx and abs(eny - kk[0][1])<=pixDy and abs(enz - kk[1])<=pixDz )]
                for kk in listCopy:
                    #if len(fullList)>1 and abs(enx - fullList[1][0][0])+abs(eny - fullList[1][0][1])<=pixD:
                    plotS[0].append(kk[0][0])
                    plotS[1].append(kk[0][1])
                    plotS[2].append(kk[1])
                    fullList.remove(kk)
                    ctp+=1

            else:
                clussA = clussAttr(ii,ctp,plotS)
                clusters.append(clussA)
                Exit = False

            ctp+=1
                        
        if plotClusters and fullList:
            colEdge_surf = tuple(random.rand(3)) #col[ctt%len(col)]
            ax.scatter(plotS[0],plotS[1],plotS[2],s = 9.8,c= 'white',marker='o',edgecolor = colEdge_surf)
            
        
        ctt += 1
        
    if plotClusters:
        fig.savefig('clusterPlot.png')
        plt.close()  
            
    return clusters

def analyzeClusters(clusterList,pixelLimit = 10,arealSpreadLimit = 49):
    '''
    Parameters
    ----------
    clusterList : list containing cluster classes
    pixelLimit : int, optional
        max number of pixels in cluster to flag. The default is 10.
    arealSpreadLimit : int, optional
        largest x-y area of clusters to flag. The default is 49.

    Returns
    -------
    badC : int -- number of clusters that fail requirements (bad clusters)
    badCluss : list -- containt cluster classes of bad clusters
    '''
    
    badC = 0
    badCluss = []
    
    for cluss in clusterList:
        if cluss.noPixels>=pixelLimit:
            arealSpread = (cluss.xSpread[1]-cluss.xSpread[0])*(cluss.ySpread[1]-cluss.ySpread[0])
            if arealSpread>=arealSpreadLimit:
                badC += 1
                badCluss.append(cluss)
    
    return badC,badCluss



if __name__ == '__main__':
    
    #Define constants used
    
    heightMin = 0.2
    heightMax = 0.61
    
    convfac = 0.010584 #pixel to um conversion factor
    
    ai = 100 #rho file number to read in for full analysis
    SNo = 50400002 #tag for rho file
    xsize = 61 #zheight for bed to read in
    ysize = 104 #ysize for bed to read in
    zsize = 104 #zsize for bed to read in
    lim = 0.1 #limit to use for valid pixels
    
    colFace = 'white' #colour of pixel faces
    colEdge = 'green' #colour of edge for surface plots animation
    elevN = 30 #elevation angle
    azimA = 45 #azimuth angle
    
    #Define file with surface cloud information
    inRead = 'ASurfaceNorms.txt' 
    #Define file with particle information for top layer
    finalName = 'particleOut1000.txt'
    
    #set bools that govern the work to be done in the script
    partpartOvl = True #use to decide if particle to particle overlap should be calculated (within new layer)
    partSurfOvl = True #use to decide if particle to surface overlap should be calculated (new layer and full bottom)
    plotNew = False #plot only new particles
    plotFull = True #plot both old surface and new particles    
    remCirc = True #used to decide whether or not to go through process of casting shadows
    makeAnima = False #make animation or movie while casting shadow
    plotShad = False #plot separate shadowed plots from each particle
    
    logFile = 'log.txt'
    fullLogfile = 'AllDirsLog.txt'
    
    dirstd = str(os.getcwd())

    #Open and read the standard rho file to calculate overlap 
    #between the particles and the bottom layer
    bnf = open("rho%dSN%d.dat" %(ai,SNo),"rb")
    fcont = bnf.read()
    size = int(len(fcont)/8)
    Vd = struct.unpack('d'*size,fcont[:size*8])
    
    xx = []
    ii = 0
    
    for z in range(zsize):
        for y in range(ysize):
            for x in range(xsize):
                if Vd[ii]>lim:
                    xx.append([z,y,x])
                ii+=1
                 
    #open and generate list that contains the bottom surface information
    
    surf = [[],[],[]]
    
    with open(inRead,'r') as finl:
        for ii in finl.readlines():
            cols = ii.split('\t')
            surf[0].append(int(cols[0]))
            surf[1].append(int(cols[1]))
            surf[2].append(int(cols[2]))          

    oldZmax, oldZhighAvg = findOldHeight(surf)
    
    #Navigate to desired directories and perform analysis:
    rootdir = '.'
    
    logFull = open(fullLogfile,'w')
    logFull.write('Bottom Layer maximum height:\t%f\t%f\nBottom Layer average height:\t%f\t%f\n' %(oldZmax,oldZmax*convfac,oldZhighAvg,oldZhighAvg*convfac))
    logFull.write('dirName\tSteadyStateChangeTest\t\tHeightTest\t\t#partsInBottomTest\t\t#badClustersTest\t\n')
    
    for dirN, dirs, files in os.walk(rootdir):
        for file in files:
            if (file=='particleOut10.txt'):

                os.chdir(dirN)
                logFull.write('%s\t' %dirN)
                print('Started analysis on', dirN)
                
                #Start analysis:
                #---------------#
                
                #Open print log file
                logOut = open(logFile,'w')

                #calculate the maximum position variation
                pathd = str(os.getcwd())
                fileList = os.listdir(pathd)
                blank = []
                for ii in fileList:
                    if 'particleOut' in ii:
                        blank.append(int(ii.replace('particleOut','').replace('.txt','')))            
                highVal = max(blank)
                if highVal == 1000:
                    a = max(positionVariation(800,900,1000))
                else:
                    a = 0
                    fileN = 'particleOut'+str(highVal)+'.txt'
                    copyfile(fileN,finalName)
                    print("Analysis file",fileN,file=logOut)

                #open and generate list that contains the full bed 
                #information for the final file to be 
                #considered as the steady state position of the new layer
                pinA = [[],[],[],[]]
                
                with open(finalName,'r') as finl:
                    for ii in finl.readlines():
                        cols = ii.split('\t')
                        pinA[0].append(float(cols[0])/convfac)
                        pinA[1].append(float(cols[1])/convfac)
                        pinA[2].append(float(cols[2])/convfac)
                        pinA[3].append(float(cols[3])/convfac)
                
                nop = len(pinA[0])
                
                print("Max x,y,z variation in position --->",a,'um',file=logOut)
                
                if a<= 1e-5:
                    print('Passed steady state Test. Moving on to next test')
                    testSteadyState = True
                else:
                    print('Failed steady state Test.')
                    testSteadyState = False       
                
                logFull.write('%s\t%f\t' %(testSteadyState,a))
                
                #Calculate height
                sphereList = makeParticles(pinA)
                newZmax, newZhighAvg = findOldHeight(sphereList)
                
                topheight = newZhighAvg - oldZhighAvg
                print("Height of new layer Info\nMaximum: --->",newZmax,'pix |',newZmax*convfac,'um',file=logOut)
                print("Average: --->",newZhighAvg,'pix |',newZhighAvg*convfac,'um',file=logOut)
                
                if heightMin<=topheight*convfac<=heightMax:
                    print('Passed Height Test. Moving on to next test')
                    testHeight = True
                else:
                    print('Failed Height Test.')
                    testHeight = False
                
                logFull.write('%s\t%f\t' %(testHeight,topheight*convfac))
                
                if plotFull:
                    #plotParticleSurf([pinA[0][0],pinA[1][0],pinA[2][0],pinA[3][0]], surf) #plots only a single particle with surface
                    print('Plotting bottom layer with new layer',file=logOut)
                    plotSpheresSurf(sphereList,surf)
                
                if plotNew:
                    print('Plotting new layer',file=logOut)
                    plotSpheres(sphereList)
            
                #Calculate overlap information
                if partpartOvl:
                    print('OverlapInfo\n','-'*len('OverlapInfo'),sep='',file=logOut)
                    OverlapInfo = Overlap(pinA,xx,partpartOvl,partSurfOvl) #Turning it off
                    print(*OverlapInfo,sep='\n',file=logOut)
                
                    #Write Overlap information to file
                    print('Writing to file',file=logOut)
                    TAGS = ['Particle: ','Particle-particle overlap: ', 'sum ', 'Particle-surface overlap: # ', '% ']
                    toFile(OverlapInfo,fname='OverlapInfoNew2.txt',filout='float', tag=TAGS)
                    
                    #analyze overlap info
                    noIN = analyzeOverlapInfo(OverlapInfo)
                    if noIN:
                        print('Failed overlay test --> There are',noIN,'particles in the bottom layer.')
                        testOnLayer = False  
                    else:
                        print('Passed overlay test --> There are no particles in the bottom layer. Moving on to next test')
                        testOnLayer = True            
                    
                    logFull.write('%s\t%d\t' %(testOnLayer,noIN))
                    
                if remCirc:
                    
                    #Get rid of pixels in the particle's shadow
                    surfTemp = resortSurf(surf)
                    print('Casting shadow\nsize of fullList before -->',len(surfTemp),file=logOut)
                    
                    #Casting shadow --- basically removing all the points in largest circle
                    #cast by the center 2D circle of the sphere
            
                    ##Testing Casting with single particle
                    #castShadow(surfTemp, [pinA[0][0],pinA[1][0],pinA[2][0],pinA[3][0]])
                    
                    #make animation...
                    if makeAnima:
                        
                        Writer = animation.writers['ffmpeg']
                        writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
                
                        fig = plt.figure()
                        ax = Axes3D(fig)
                                        
                        ax.set_xlim([0,ysize])
                        ax.set_ylim([0,ysize])
                        ax.set_zlim([0,ysize])
                        ax.view_init(elev=elevN,azim=azimA)
                        
                        #making surface an animation
                        
                        scat = ax.scatter([],[],[],s = 9.8,c=colFace,marker='o',edgecolor=colEdge)
                        
                        anim = animation.FuncAnimation(fig,animaUpdate,nop+1,interval=200)
                        
                        # save a gif of the animation using the writing package from magick
                        anim.save('shadow.mp4', dpi=80, writer = writer) #movie
                        #anim.save('shadow.gif', dpi=80, writer = animation.PillowWriter(80)) #gif
                        
                        plt.close()
                        
                        surfPlot = resortSurf(surfTemp,1) 
                        figName = 'AfterCasting.png'
                        plotSurf(surfPlot,figName)
                        
                    else:
                        
                        ctp = 0
                        surfPlot = resortSurf(surfTemp,1) 
                        plotSurf(surfPlot)
                        for xc,yc,zc,rr in zip(pinA[0],pinA[1],pinA[2],pinA[3]):
                            castShadow(surfTemp, [xc,yc,zc,rr])
                            print('Done shadowing particle:', ctp, 'size of list now -->', len(surfTemp),file=logOut)
                            
                            if plotShad:
                                #Making plots
                                surfPlot = resortSurf(surfTemp,1) 
                                figName = 'Beautplot'+str(ctp)+'.png'
                                plotSurf(surfPlot,figName)
                                #print('Done plotting shadowed particle')
                            
                            ctp+=1      
                    
                    print('Extracting Top surface',file=logOut)
                    extractTop(surfTemp)
                    print('Done extracting top, size of list now -->', len(surfTemp),file=logOut) 
                    surfPlot = resortSurf(surfTemp,1) 
                    figName = 'Beautplot_Top.png'
                    plotSurf(surfPlot,figName)
                    print('Done plotting top',file=logOut) 
                    print('Started Clustering',file=logOut)
                    topClusters = cluster(surfTemp)
                    print('Done clustering, size of list now -->', len(surfTemp),file=logOut)
                    print('Cluster Info\n','-'*len('Cluster Info'),sep='',file=logOut)
                    print(*topClusters,file=logOut)
                    badClusters = analyzeClusters(topClusters,pixelLimit = 100, arealSpreadLimit = 400)
                        
                    if badClusters[1]:
                        print('Failed spread test.',len(badClusters[1]),'clusters over set limits.')
                        testSpread = False        
                    else:
                        print('Passed spread Test. No clusters over set limits')
                        testSpread = True
                        
                    logFull.write('%s\t%d\n' %(testSpread,len(badClusters[1])))     
                    
                #Close log file
                logOut.close()
                
                print('Done')
            
                os.chdir(dirstd)
    
    #close full log file
    logFull.close()
