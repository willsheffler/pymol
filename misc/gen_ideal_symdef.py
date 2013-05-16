#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pymol_util as pu
import xyzMath as xyz
from math import sqrt

"""
ISOCAHEDRAL:
icos virts are 5folds:
 (0,±1,±ø),(±1,±ø,0),(±ø,0,±1)], ø is golden ratio
dodec virts are 3folds:
 (±1, ±1, ±1),(0, ±1/ø, ±ø),(±1/ø, ±ø, 0),(±ø, 0, ±1/ø)
icosidode virts are 2folds:
 (0,0,±ø)
 (0,±ø,0)
 (±ø,0,0)
 (±1/2, ±ø/2, ±(1+ø)/2)
 (±ø/2, ±(1+ø)/2, ±1/2)
 (±(1+ø)/2, ±1/2, ±ø/2)
"""

r = (1.+sqrt(5.))/2. # golden ratio

def ics2():
	edgelen = 2.0
	for x,y,z in ((2*r,0,0),(-2*r,0,0),(0,2*r,0),(0,-2*r,0),(0,0,2*r),(0,0,-2*r),
	              (+1,+r,+(1+r)) , (+r,+(1+r),+1) , (+(1+r),+1,+r) ,
	              (+1,+r,-(1+r)) , (+r,+(1+r),-1) , (+(1+r),+1,-r) ,
	              (+1,-r,+(1+r)) , (+r,-(1+r),+1) , (+(1+r),-1,+r) ,
	              (+1,-r,-(1+r)) , (+r,-(1+r),-1) , (+(1+r),-1,-r) ,
	              (-1,+r,+(1+r)) , (-r,+(1+r),+1) , (-(1+r),+1,+r) ,
	              (-1,+r,-(1+r)) , (-r,+(1+r),-1) , (-(1+r),+1,-r) ,
	              (-1,-r,+(1+r)) , (-r,-(1+r),+1) , (-(1+r),-1,+r) ,
	              (-1,-r,-(1+r)) , (-r,-(1+r),-1) , (-(1+r),-1,-r) ):	              
		yield xyz.Vec(x,y,z)/edgelen

def ics3():
	edgelen = sqrt(5.)-1.
	for x,y,z in ((+1,+1,+1) , (+1,+1,-1) , (+1,-1,+1) , (+1,-1,-1),
	              (-1,+1,+1) , (-1,+1,-1) , (-1,-1,+1) , (-1,-1,-1),
		          (0,+1/r,-r) , (+1/r,-r,0) , (-r,0,+1/r),
	              (0,+1/r,+r) , (+1/r,+r,0) , (+r,0,+1/r),
	              (0,-1/r,-r) , (-1/r,-r,0) , (-r,0,-1/r),
	              (0,-1/r,+r) , (-1/r,+r,0) , (+r,0,-1/r)):
		yield xyz.Vec(x,y,z)/edgelen

# get 2-3 pseudo-edge dis
# m = 999
# for v in ics3():
# 	for u in ics2():
# 		m = min(u.distance(v),m)
# print m

def ics5():
	edgelen = 2.0
	for x,y,z in ((0,+1,+r) , (+1,+r,0) , (+r,0,+1),
	              (0,+1,-r) , (+1,-r,0) , (-r,0,+1),
	              (0,-1,+r) , (-1,+r,0) , (+r,0,-1),
	              (0,-1,-r) , (-1,-r,0) , (-r,0,-1)):
		yield xyz.Vec(x,y,z)/edgelen


def nbrs(p1,crd,edgelen=1.0):
	for p2 in crd:
		if edgelen-0.001 <= p1.distance(p2) <= edgelen+0.001:
			yield p2

# def showvec(c,lab=None):
# def showvecfrompoint(a, c, lbl):

def showsym(virts,nbrvirts=None,edgelen=1.0,lab="",R=xyz.Imat):
	if not nbrvirts: nbrvirts = virts
	print "virtual_coordinates_start"
	for i,virtex in enumerate(virts()):
		ax1 = R*virtex.normalized()
		# print ax1
		# if pu.inpymol(): pu.showvec(4.*ax1,lab+"_ax1_"+str(i))
		for j,nbr in enumerate(nbrs(virtex,nbrvirts(),edgelen)):
			ax2 = R*xyz.projperp(virtex,nbr-virtex).normalized()
			if 0==j: print "xyz C%02i   %+16.14f,%+16.14f,%+16.14f %+16.14f,%+16.14f,%+16.14f 0,0,0"%(i,ax1.x,ax1.y,ax1.z,ax2.x,ax2.y,ax2.z)
			if 0==j: print "xyz P%02i   %+16.14f,%+16.14f,%+16.14f %+16.14f,%+16.14f,%+16.14f 0,0,0"%(i,ax1.x,ax1.y,ax1.z,ax2.x,ax2.y,ax2.z)
			print "xyz B%02i_%i %+16.14f,%+16.14f,%+16.14f %+16.14f,%+16.14f,%+16.14f 0,0,0"%(i,j,ax1.x,ax1.y,ax1.z,ax2.x,ax2.y,ax2.z)
			if pu.inpymol():
				pu.showvecfrompoint(ax1,1.5*ax1,lab+"_ax1_"+str(i)+"_"+str(j))
				pu.showvecfrompoint(ax2,1.5*ax1,lab+"_ax2_"+str(i)+"_"+str(j))
	print "virtual_coordinates_stop"

print "============================= I5 ============================="
showsym(ics5,lab="I5",R=xyz.rotation_matrix_degrees(xyz.Vec(1,1,1),-44.477512185929)) # align dodec to icos
print "============================= I3 ============================="
showsym(ics3,lab="I3")
print "============================= I2 ============================="
showsym(ics2,ics3,0.587785252292,"I2")
pu.cmd.zoom()












