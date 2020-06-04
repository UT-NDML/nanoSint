import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from array import array as ARRAY

class laserProps:
    def __init__(self,R,x,y,z,T):
        self.R = R
        self.x = x
        self.y = y
        self.z = z
        self.T = T

def fillBox(array,boxSize,laser):
    
    for iz in range(boxSize[2]):
        for iy in range(boxSize[1]):
            for ix in range(boxSize[0]):
                l = iz*boxSize[0]*boxSize[1] + iy*boxSize[0] + ix
                if (ix-laser.x)**2+(iy-laser.y)**2+(iz-laser.z)**2<=laser.R*laser.R:
                    array[l] = laser.T
    
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

def writeToFile(fileName,array):
    with open(fileName,'w') as fill:
        for ii in array:
            fill.write('%d\n' %ii)

def readFromFile(fileName):
    array = []
    with open(fileName,'r') as fill:
        for ii in fill.readlines():
            array.append(eval(ii))
    return array

boxSize = [94,104,104]
laser1 = laserProps(30,94,50,50,500)
laser2 = laserProps(40,94,50,50,450)

if True:
    array = [0 for x in range(boxSize[0]*boxSize[1]*boxSize[2])]
    
    print('Start')
    fillBox(array,boxSize,laser2)
    fillBox(array,boxSize,laser1)
    
    print('Done Filling')
    colors={'450':'y','500':'r'}
    
    fig = plt.figure()
    ax = Axes3D(fig)
    outP = plotLaser(boxSize,array,colors)
    ax.set_xlim([0,100])
    ax.set_ylim([0,100])
    print('Start Plotting')
    ax.scatter(outP[0],outP[1],outP[2],s = 9.8,c='w',marker='o',edgecolor=outP[3])
    ax.view_init(elev=50,azim=45)
    fig.savefig('TemperatureProfDat.png')
    plt.close()

    #if writing to text file    
    #writeToFile('tempProfile2.txt',array)
    
    #if writing to .dat file
    fill = open('tempProfile.dat','wb')
    s = ARRAY('d',array)
    s.tofile(fill)
    fill.close()    
    print('Done')
    
else:    
    array = readFromFile('tempProfile.txt')
    
    colors={'450':'y','500':'r'}
    
    fig = plt.figure()
    ax = Axes3D(fig)
    outP = plotLaser(boxSize,array,colors)
    ax.set_xlim([0,100])
    ax.set_ylim([0,100])
    print('Start Plotting')
    ax.scatter(outP[0],outP[1],outP[2],s = 9.8,c='w',marker='o',edgecolor=outP[3])
    ax.view_init(elev=50,azim=45)
    fig.savefig('TemperatureProfCheck.png')
    plt.close()
    
    print('Done')