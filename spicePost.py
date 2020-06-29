from sys import argv

resP = int(argv[1])

R = []

ctStart = 0
with open('spiceOut.txt','r') as fill:
    for ii in fill.readlines():
        if ctStart>1:
            iline = ii.split()
            V = eval(iline[1])   
            I = eval(iline[2])
            R.append(abs(V/I))
            ctStart = 0
        if ii[0:5] == '-----':
            ctStart += 1            
        
if resP: print(max(R)) 
else: 
    for ii in R:
        print(ii,end = ' ')
