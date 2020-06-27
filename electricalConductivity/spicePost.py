from sys import argv

resP = int(argv[1])

R = []

ctStart = 0
with open('spiceOutX.txt','r') as fill:
    for ii in fill.readlines():
        if ctStart>1:
            iline = ii
            break
        if ii[0:5] == '-----':
            ctStart += 1            
        
iline = iline.split()

V = eval(iline[1])   
I = eval(iline[2])
R.append(abs(V/I))   

ctStart = 0
with open('spiceOutY.txt','r') as fill:
    for ii in fill.readlines():
        if ctStart>1:
            iline = ii
            break
        if ii[0:5] == '-----':
            ctStart += 1            
        
iline = iline.split()

V = eval(iline[1])   
I = eval(iline[2])
R.append(abs(V/I))     

if not resP: print(max(R)) 
else: print(R[0],R[1])
