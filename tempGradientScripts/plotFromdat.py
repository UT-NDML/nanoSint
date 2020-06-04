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

colbd = 'black'
colfc = 'white'

lim = 0.1

xsize = 94
ysize = 104
zsize = 104

fileName = "tempProfile.dat"#"tempAll0.dat"#"rho0SN504.dat"
plotName = "tempOriginalCheck3.png"
colors={'450':'y','500':'r'}

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