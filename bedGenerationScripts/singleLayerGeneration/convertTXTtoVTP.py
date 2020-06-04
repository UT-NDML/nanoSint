with open('particle_input.txt','r') as part:
	allparts = part.readlines()

Allarray = []

for pt in allparts:
	pt = pt.split('\t')
	del(pt[8:],pt[4])
	pt = [float(ii) for ii in pt]
	
	Allarray.append(pt)

outf = open('input.vtp','w')

timet = 0.0
noparts = len(Allarray)

outf.write('<?xml version="1.0"?>\n')
outf.write('<!--Time = %fs -->\n' %timet)
outf.write('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
outf.write("<PolyData>\n")
outf.write('<Piece NumberOfPoints="%d" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="0">\n' %noparts)
outf.write("<Points>\n")
outf.write('<DataArray type="Float32" Name="Position" NumberOfComponents="3" format="ascii">\n')

for ii in Allarray:
    outf.write('%f\t%f\t%f\t' %(ii[0],ii[1],ii[2]));

outf.write("\n")

outf.write("</DataArray>\n</Points>\n")
outf.write('<PointData Scalars="Diameter" Vectors="Velocity">\n')
outf.write('<DataArray type="Float32" Name="Diameter" format="ascii">\n')

for ii in Allarray:
    outf.write('%f\t' %(2*ii[3]));

outf.write("\n")

outf.write("</DataArray>\n")
outf.write('<DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="ascii">\n')

for ii in Allarray:
    outf.write('%f\t%f\t%f\t' %(ii[4],ii[5],ii[6]));

outf.write("\n")

outf.write("</DataArray>\n")
outf.write('<DataArray type="Float32" Name="ID" format="ascii">\n')

for ii in range(1,noparts+1):
    outf.write('%d\t' %ii);

outf.write("\n</DataArray>\n")
outf.write("</PointData>\n<CellData></CellData>\n<Verts></Verts>\n<Lines></Lines>\n")
outf.write("<Strips></Strips>\n<Polys></Polys>\n</Piece>\n</PolyData>\n</VTKFile>\n")

outf.close()
