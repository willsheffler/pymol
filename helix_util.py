# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import sys,os,inspect,functools
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path: sys.path.append(newpath)
import string,re,gzip,itertools
import collections
from sym_util import *
import operator as op
from xyzMath import *
from itertools import product,ifilter


def xform_covers_all_coms( axis, cen, coms, xform, nfold, show=False ):
	#TODO: add NFOLD!!!!!!!!!!!!!!
	syms = [ RAD( axis, i*360.0/nfold, cen ) for i in range(nfold) ]
	seenit = [False] * len(coms)
	for i,in_com in enumerate(coms):
		hcom = coms[0]
		for j in range(222):
			for sym in syms:
				testcom = sym*hcom
				if show: showsphere( testcom )
				if in_com.distance(testcom) < 0.1:
					# print in_com.distance(hcom)
					seenit[i] = True
					break
			hcom = xform * hcom
	return all(seenit)

def get_correction_angle( axis, cen, coms, unit_xform, error ):
	correction_ang = None
	nturns = None
	correction_ang_axis_dis = 0.0
	for i in range(len(coms)-1,1,-1):
		in_com = coms[i]
		print "try to get correction_ang", i, in_com, len(coms)
		hcom = coms[0]
		axis_dis = in_com.dot(axis)
		if axis_dis < correction_ang_axis_dis:
			continue
		for j in range(222):
			# if i==len(coms)-1: showsphere( hcom )
			if in_com.distance(hcom) < 2.0:
				# print coms[-1].distance(hcom)
	 			v1_0 = in_com - cen
	 			v2_0 = hcom - cen
	 			v1 = ( v1_0-v1_0.dot(axis)*axis ).normalized()
				v2 = ( v2_0-v2_0.dot(axis)*axis ).normalized()
				# print v1.dot(axis)
				# print v2.dot(axis)			
				correction_ang = acos( v1.dot(v2) )
				correction_ang_axis_dis = axis_dis
				# TODO: check this logic!
				if v1.cross(v2).dot(axis) < 0:
					correction_ang *= -1.0
				print "correction_ang",correction_ang,i,axis_dis
				nturns = j
				break
			hcom = unit_xform * hcom
	return correction_ang, nturns, error


def determine_helix_geometry( sele='vis', show=0 ):
	"""
delete all; load /Users/sheffler/Downloads/N4_C1_DR53_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=40)
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N2_C1_DR04_065.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N2_C3_DR04_003.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C2_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C4_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N4_C3_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N5_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N4_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)	
"""

	error = None

	axis, ang1, cen = getrelframe('((%s) and name ca and chain B)'%sele,'((%s) and name ca and chain A)'%sele
		).rotation_axis_center()
	com_all = com( "((%s) and name ca)"%sele )
	if axis.dot(com_all-cen) < 0:
		axis = Vec(0,0,0)-axis
		ang1 = 2.0*math.pi - ang1

	coms = [ com( "((%s) and name ca and chain %s)"%(sele,c) ) for c in cmd.get_chains(sele) ]	
		# make sure sorted along axis
		# coms_tosort = [ ( axis.dot(xyz-cen),xyz) for xyz in coms ]
		# coms = [ xyz for t,xyz in sorted( coms_tosort ) ]

	print "axis",axis, "sub2 angle",ang1*180.0/math.pi

	# get unitary translation along axix
	mintrans = 9e9
	for c1 in coms:
		for c2 in coms:
			if c1 is not c2:
				trans = abs( axis.dot(c2-c1) )
				if trans > 0.5: # assume is < 0.5, is cyclic sym related subunit
					mintrans = min( mintrans, trans )
	print "mintrans along axis is", mintrans


	unit_trans = 0
	for xyz in coms[1:]:
		mult = round( (xyz-coms[0]).dot(axis) / mintrans )
		unit_trans += (xyz-coms[0]).dot(axis) / mult
		# print (xyz-coms[0]).dot(axis), mult, unit_trans, (xyz-coms[0]).dot(axis) / mult
	unit_trans /= len(coms)-1 
	sub2_trans = (coms[1]-coms[0]).dot(axis)
	print mintrans, unit_trans, sub2_trans, sub2_trans/unit_trans

	# check that each translation is even multiple of unit_trans
	for xyz in coms[1:]:
		test_mult = (xyz-coms[0]).dot(axis) / unit_trans
		if abs( test_mult - round( test_mult ) ) > 0.01:
			print "(xyz-coms[0]).dot(axis) / unit_trans == ", test_mult
			assert False
	unit_ang1 = ang1 / (sub2_trans/unit_trans)
	unit_xform1 = rotation_around( axis, unit_ang1, cen )
	unit_xform1.t += unit_trans*axis

	# now recompute rotation angle based on last subunit
	# this corrects small errors from pdb coordinates or whatever...
	# or maybe unit_xform1.t += unit_trans*axis is somehow wrong and I'm dumb...
	correction_ang, nturns, error = get_correction_angle( axis, cen, coms, unit_xform1, error )
	if not correction_ang and not error:
		error = "correction_ang is None"
	if correction_ang:
		print "correcting unit_ang1 and ang1"
		unit_ang1 = unit_ang1-correction_ang/nturns
		ang1 = unit_ang1 * (sub2_trans/unit_trans)
	unit_xform1 = rotation_around( axis, unit_ang1, cen )
	unit_xform1.t += unit_trans*axis

	# angs        = [ang1]
	unit_angs   = [unit_ang1]
	unit_xforms = [unit_xform1]
	for i in range(1,4):
		ang2 = 2.0*i*math.pi + ang1
		# this is crappy duplicated code....
		unit_ang2 = ang2 / (sub2_trans/unit_trans)
		unit_xform2 = rotation_around( axis, unit_ang2, cen )
		unit_xform2.t += unit_trans*axis
		correction_ang2, nturns2, error2 = get_correction_angle( axis, cen, coms, unit_xform2, error )
		if not correction_ang2 and not error2:
			error = "correction_ang2 is None"
		if correction_ang2:
			print "correcting unit_ang and ang"
			unit_ang2 = unit_ang2-correction_ang2/nturns2
			ang2 = unit_ang2 * (sub2_trans/unit_trans)
		unit_xform2 = rotation_around( axis, unit_ang2, cen )
		unit_xform2.t += unit_trans*axis
		# angs       .append(      ang2 )
		unit_angs  .append( unit_ang2 )
		unit_xforms.append( unit_xform2)

	# now figure out Nfold based on coverage
	unit_xform = None
	unit_ang = None
	# ang = None
	if not error:
		error = "error determining helix symmetry (and check for ang errors)"

	# check up to nfold 10
	for nfold in range(1,11):
		for i in range( len(unit_xforms) ):
			unit_xform = unit_xforms[i]
			unit_ang   = unit_angs[i]
			# ang = angs[i]
			print "checking symmetry Nfold",nfold,"ang_option",i
			if xform_covers_all_coms( axis, cen, coms, unit_xform, nfold ):
				error = None
				break
		if not error:
			break
	if not error: print "nfold is",nfold

	# optional graphics
	if show > 0:
		cgo,cgo2 = [],[]
		xyz = coms[0]
		for i in range(int(math.ceil(show/nfold))):
			syms = [RAD( axis, i*360.0/nfold, cen ) for i in range(nfold)]
			for sym in syms:
				cgo.extend( cgo_sphere(sym*xyz) )
			xyz = unit_xform * xyz
		for c in coms:
			cgo2.extend( cgo_sphere( c, col=(1,0,0) ) )			
		cgo.extend( cgo_sphere( cen ) )
		cgo.extend( cgo_segment( cen, cen+axis*(xyz.dot(axis)) ) )
		cmd.load_cgo(cgo,"determine_helix_geometry")
		cmd.load_cgo(cgo2,"existing_coms")		

	print "unitary rotation angle", unit_ang*180.0/math.pi

	if error:
		print error
		return None

	return axis, cen, unit_ang, nfold, unit_trans, unit_xform


def make_helix_symdef( sele='vis', show=0 ):
	"""
delete all; load /Users/sheffler/Downloads/N4_C1_DR53_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=40)
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N2_C1_DR04_065.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N2_C3_DR04_003.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C2_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C4_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N4_C3_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N5_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)	
delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N4_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)	
	"""
	cgo = []

	axis, cen, unit_ang, nfold, unit_trans, unit_xform = determine_helix_geometry( sele, show=show )

	coms = [ com( "((%s) and name ca and chain %s)"%(sele,c) ) for c in cmd.get_chains(sele) ]	

	# now generate symfile with xyzs on axis in intervals linked to subunits of same height
	# will need to generate below and above

	maxtrans = max( c.dot(axis) for c in coms )

	hcom,hcom2 = coms[0],coms[0]
	for j in range(100):
		if j*unit_trans - 0.1 > maxtrans:
			break
		syms = [RAD( axis, i*360.0/nfold, cen ) for i in range(nfold)]
		for sym in syms:
			a  = cen+j*unit_trans*axis
			a2 = cen-j*unit_trans*axis		
			x = axis
			cgo.extend( cgo_sphere( a ) )
			if j > 0: cgo.extend( cgo_sphere( a2 ) )
			for i,in_com in enumerate(coms):
				if in_com.distance(sym*hcom) < 0.1:
					cgo.extend( cgo_sphere( sym*hcom ) )
					y = sym*hcom - a
					y2 = hcom2 - a2
					# print "xyz", 
					cgo.extend( cgo_segment( a, a+y ) )				
					if j > 0:
						cgo.extend( cgo_segment( a2, a2+y2 ) )								
						cgo.extend( cgo_sphere( a2+y2 ) )
					break
		hcom = unit_xform * hcom
		hcom2 = unit_xform.inverse() * hcom2		

	if show: cmd.load_cgo( cgo, "make_helix_symdef" )
	print "make_helix_symdef done"



def makeh(sele='vis',n=30):
	axis, cen, unit_ang, nfold, unit_trans, unit_xform = determine_helix_geometry( sele, show=n )

	n = n / nfold
	cmd.delete('helix')
	v = cmd.get_view()
	cmd.create('tmp',sele+' and chain A and name n+ca+c')
	count = 0
	for nf in range(nfold):
		x = Xform()
		print nf, x.pretty()
		for i in range(n):
			cmd.create('Htmp%i'%count,'tmp')
			xform('Htmp%i'%count,x)
			cmd.alter('Htmp%i'%count,"chain='%s'"%ROSETTA_CHAINS[count])
			count += 1
			x = x * unit_xform
	cmd.create("HELIX",'Htmp*')
	cmd.delete("Htmp*")
	cmd.delete('tmp')
	cmd.hide('ev','HELIX')
	cmd.show('lines','helix')
	util.cbc('HELIX')
	cmd.set_view(v)
cmd.extend('makeh',makeh)




def load_tests(loader, tests, ignore):
	tests.addTests(doctest.DocTestSuite())
	return tests

if __name__ == '__main__':
   import doctest
   r = doctest.testmod()
   print r


