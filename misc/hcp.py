import sys,os,inspect
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
newpath = "/".join(newpath.split("/")[:-1])
if not newpath in sys.path:
		sys.path.append(newpath)
from xyzMath import Vec
from math import sqrt

def hcp(i,j,k,r=1.0):
	return r * Vec(
		2.0*i + (j%2)     + (k%2),
		  (j + (k%2)/3.0)  * sqrt(3.0),
		k * sqrt(6.0) / 3.0 * 2.0
	)

def printatom(v,i=1):
	return "ATOM %6i  CD  GLN D %3i     %7.3f %7.3f %7.3f  1.00 75.17           C"%(i,i,v.x,v.y,v.z)

print os.getcwd()
with open("test.pdb","w") as o:
	N = 9
	count = 1
	for i in range(N):
		for j in range(N):
			for k in range(N):
				o.write( printatom( hcp(i,j,k,3.0), count ) + "\n" )
				count += 1