from xyzMath import *
import math
with open('/Users/sheffler/tmp/euler_of_uniform_rots.txt','w') as o:
	for m in randrot(1000000):
		e = m.euler_angles()
		o.write("%6.3f %6.3f %6.3f\n"%(e.x,e.y,e.z))

