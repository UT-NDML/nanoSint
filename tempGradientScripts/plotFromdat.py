import struct
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotLaser(boxSize,array,colors):
    
    pos = [[],[],[],[]]
    
    for iz in range(boxSize[2]):
        for iy in range(boxSize[1]):
            for ix in range(boxSize[0]):
                l = iz*boxSize[0]*boxSize[1] + iy*boxSize[0] + ix
                if array[l]:
                    pos[0].append(ix)
                    pos[1].append(iy)
                    pos[2].append(iz)
                    pos[3].append(colors[str(array[l])])
    return pos

def procPlot(array,pos,colors):
    plotA = [[],[],[],[]]
    for ctr,ii in enumerate(array):
        if ii:
            for ij in range(3): plotA[ij].append(pos[ij][ctr])
            plotA[3].append(colors[str(int(ii))])
    return plotA

colors={'450':'y','500':'r'}

colbd = 'black'
colfc = 'white'

lim = 0.1

xsize = 94
ysize = 104
zsize = 104

fullPlots = False
sidePlots = True

psizey = 8
psizez = 8

zheight = int(zsize/psizez)
yheight = int(ysize/psizey)
psize = psizey*psizez

area = yheight*xsize
areay = zheight*xsize
fullarea = xsize*ysize

plotName = 'plotldn.png'
nameStart = "templdn"


if sidePlots:
    fig = plt.figure()
    ax = Axes3D(fig)
    
    ax.set_xlim([0,ysize])
    ax.set_ylim([0,ysize])
    ax.set_zlim([0,ysize])
    ax.view_init(elev=40,azim=50)
        
    for prank in range(psize):
        zpos = prank//psizey
        ypos = prank%psizey
        
        ffile = nameStart+str(prank)+".txt"
        
        array = []; pos = [[],[],[]]
        with open(ffile,'r') as fil:
            filA = fil.read()
            
        for ii in filA.split('\n')[1:-1]:
            ij = ii.split()[-1]
            array.append(eval(ij))
        
        if nameStart == "tempdown":
            zz = zpos*zheight-1
            for y in range(yheight):
                yy = y + ypos*yheight
                for x in range(xsize):
                    pos[0].append(zz)
                    pos[1].append(yy)
                    pos[2].append(x)
        elif nameStart == "tempup":
            zz = (zpos+1)*zheight
            for y in range(yheight):
                yy = y + ypos*yheight
                for x in range(xsize):
                    pos[0].append(zz)
                    pos[1].append(yy)
                    pos[2].append(x)
        elif nameStart == "templeft":
            yy = ypos*yheight-1
            for z in range(zheight):
                zz = z + zpos*zheight
                for x in range(xsize):
                    pos[0].append(zz)
                    pos[1].append(yy)
                    pos[2].append(x)
        elif nameStart == "tempright":
            yy = (ypos+1)*yheight
            for z in range(zheight):
                zz = z + zpos*zheight
                for x in range(xsize):
                    pos[0].append(zz)
                    pos[1].append(yy)
                    pos[2].append(x)                    
        elif nameStart == "templdn":
            yy = ypos*yheight-1
            zz = zpos*zheight-1
            for x in range(xsize):
                pos[0].append(zz)
                pos[1].append(yy)
                pos[2].append(x) 
        elif nameStart == "temprdn":
            yy = (ypos+1)*yheight
            zz = zpos*zheight-1
            for x in range(xsize):
                pos[0].append(zz)
                pos[1].append(yy)
                pos[2].append(x)  
        elif nameStart == "templup":
            yy = ypos*yheight-1
            zz = (zpos+1)*zheight
            for x in range(xsize):
                pos[0].append(zz)
                pos[1].append(yy)
                pos[2].append(x) 
        elif nameStart == "temprup":
            yy = (ypos+1)*yheight
            zz = (zpos+1)*zheight
            for x in range(xsize):
                pos[0].append(zz)
                pos[1].append(yy)
                pos[2].append(x) 
                  
        toPlot = procPlot(array,pos,colors)
        ax.scatter(toPlot[0],toPlot[1],toPlot[2],s=9.8,c= colfc,marker='o',edgecolor = toPlot[3])
    
        print(prank,"done")
        
    fig.savefig(plotName)
    plt.close()


if fullPlots:
    fileName = "tempAll0.dat"#"rho0SN504.dat"#"tempProfile.dat"#
    plotName = "tempCheck.png"#"rhoCheck.png"#"tempOriginalCheck.png"#
    
    bnf = open(fileName,"rb")
    fcont = bnf.read()
    fullsize = int(len(fcont)/8)
    Vd = struct.unpack('d'*fullsize,fcont[:fullsize*8]) #contains float information for the data file
    bnf.close()
    
    Zd = []; Yd = []; Xd = []
    colArr = []
    
    ct = 0
    for z in range(zsize):
        for y in range(ysize):
            for x in range(xsize):
                if Vd[ct]>lim:
                    Zd.append(x)
                    Yd.append(y)
                    Xd.append(z)
                    if str(int(Vd[ct])) in colors.keys():
                        colArr.append(colors[str(int(Vd[ct]))])
                    
                ct += 1
    
    fig = plt.figure()
    ax = Axes3D(fig)
    
    ax.set_xlim([0,ysize])
    ax.set_ylim([0,ysize])
    ax.set_zlim([0,ysize])
    
    if colArr: colbd = colArr[:]
    
    ax.scatter(Xd,Yd,Zd,s=9.8,c= colfc,marker='o',edgecolor = colbd)
    ax.view_init(elev=40,azim=50)
    
    fig.savefig(plotName)
    plt.close()

'''
prank = 11
zpos = prank//psizey
ypos = prank%psizey

disp = ypos*area + (zpos*zheight-1)*fullarea

arr = []
for ii in range(area):
    arr.append(disp+ii)

prank = 3
zpos = prank//psizey
ypos = prank%psizey

checkarr = []
for z in range(zheight):
    ti = xsize*yheight*z
    disp = z*fullarea + ypos*area + zpos*ysize*areay
    for ii in range(area):
        checkarr.append(disp+ii)    
'''