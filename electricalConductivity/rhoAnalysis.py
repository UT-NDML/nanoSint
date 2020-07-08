Alist = list(range(0,12,2)) + list(range(20,210,10)) + list(range(400,1200,200)) + list(range(1250,3250,250)) + [3150,3300]

#open file to write density information into
writeName = 'densityInfo.txt'
fill = open(writeName,'w')
fill.write('timeStep\tdensity(kg/m^3)\n')
    
for ai in Alist:
    str2write = str(ai) + '\t'
    
    fileName = 'fullT'+str(ai)+'Out.log'
    infoDic = {}
    with open(fileName,'r') as fild:
        for ii in fild.readlines():
            if ai == 0:
                if ii[0:4] == 'init':
                    lin = ii.split('=')
                    init_density = eval(lin[-1])
                    str2write += str(init_density) + '\n'
                if ii[0:6] == 'sumrho':
                    lin = ii.split('=')
                    sumrho_init = eval(lin[-1])
               
            else:
                if ii[0:6] == 'sumrho':
                    lin = ii.split('=')
                    sumrho = eval(lin[-1])
                    rel_rho_change = (sumrho - sumrho_init)/sumrho
                    density = init_density/(1-rel_rho_change)
                    str2write += str(density) + '\n'
                    
            if ii[0] == 'p':
                lin = ii.split('|')
                fl1 = lin[0]; fl2 = lin[-1]; fl3 = lin[-2]
                fl1 = fl1.replace('rad = ','')
                fl1 = fl1.split()
                keyN = fl1[0]+fl1[3]
                fl2 = fl2.replace('neck radius','').replace('neck length','').replace('=','').split()
                fl3 = fl3.replace('cenD =','')
                #props --> neck radius, neck length, radius p0, radius p1, cenD
                props = [eval(fl2[0]),eval(fl2[1]),eval(fl1[1]),eval(fl1[-1]),eval(fl3)]
                infoDic[keyN] = props
   
    fill.write('%s' %str2write)