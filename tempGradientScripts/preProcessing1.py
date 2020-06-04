from numpy import sqrt
import time

def readFromFile(fileName):
    array = []
    with open(fileName,'r') as fill:
        for ii in fill.readlines():
            array.append(eval(ii))
    return array

def styleFromFile(fileName):
    array = [[],[]]
    with open(fileName,'r') as fill:
        for ctr,ij in enumerate(fill.readlines()):
            for ii in ij.split():
                if ctr==0: array[ctr].append(eval(ii))
                else: array[1].append(eval(ii))
    return array


starttime = time.time()

boxSize = [100,100,100]

typeProf = 'style' #options ('mem','style','create')

if typeProf=='mem':
    print('ALL PIXELS IN ARRAYS')
    array = readFromFile('tempProfile.txt')
elif typeProf=='style':
    print('HYBRID IN ARRAYS')
    array = styleFromFile('tempStyle.txt')
else:
    print('NO ARRAYS CREATE IN LOOP')
    
midtime = time.time()
print('Time to initialize:')
print(midtime-starttime,'s')

for z in range(boxSize[2]):
    for y in range(boxSize[1]):
        for x in range(boxSize[0]):
            l = z*boxSize[0]*boxSize[1] + y*boxSize[0] + x
            if typeProf=='mem':
                temp = array[l]
            elif typeProf=='style':
                rad = sqrt((x-array[0][0])**2+(y-array[0][1])**2+(z-array[0][2])**2)
                temp = array[1][int(rad>array[0][3])+int(rad>array[0][4])]
            else:
                rad = sqrt((x-50)**2+(y-50)**2+(z-100)**2)
                if rad<=30: temp = 500
                elif rad<=40: temp = 450
                else: temp = 0

print('Time to retreive temperatures:')
print(time.time()-midtime,'s')
