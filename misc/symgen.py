# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import os,sys,inspect
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path: sys.path.append(newpath)
import string,re,gzip,itertools
from sym_util import *
from xyzMath import *

# delete all; run /Users/sheffler/pymol/symgen.py; symgentest(SYM="P6",NFOLD=(2,3,6),WHICHRODS=(1,0,1,),MOVASYM=1.0,SYMGEN=6,SKEW=0.6);

def closest_around(mov,tgt,axs,nfold,cen):
	mind = 9e9
	mini = -1
	for i in range(nfold):
		x = rotation_around_degrees(axs,i*360.0/nfold,cen)
		d = tgt.distance_squared(x*mov)
		print i,d
		if d < mind:
			mind = d
			mini = i
	print mini,mind
	print
	return rotation_around_degrees(axs,mini*360.0/nfold,cen)*mov

def symgentest(SYM,NFOLD,WHICHRODS=(1,1,1),MOVASYM=0.7,SYMGEN=4,SKEW=0.5,SKEW0=0.5):
	STRUTRAD = 0.02
	r0,r1,r2 = ([randvecball().length()*(abs(SKEW)+SKEW0) for i in range(3)])
	if SKEW < 0: r1 *= -1
	print r0
	print r1
	print r2
	scale = sorted((r0,r1,r2))[2]*1.0
	r0,r2,r1 = r0/scale,r1/scale,r2/scale
	a0,a1,a2 =    SYMAXIS[SYM,NFOLD[0]],                        SYMAXIS[SYM,NFOLD[1]],                        SYMAXIS[SYM,NFOLD[2]]
	c0,c1,c2 = r0*SYMAXIS[SYM,NFOLD[0]]+SYMCEN[SYM,NFOLD[0]],r1*SYMAXIS[SYM,NFOLD[1]]+SYMCEN[SYM,NFOLD[1]],r2*SYMAXIS[SYM,NFOLD[2]]+SYMCEN[SYM,NFOLD[2]]

	bo0,bo1,bo2 = (random.uniform(0.05,0.10)*projperp(v,randvec()).normalized() for v in (c0,c1,c2))

	bc0,bc1,bc2 = c0+bo0,c1+bo1,c2+bo2
	br0,br1,br2 = (random.uniform(1,1)*r.length() for r in (bo0,bo1,bo2))
	# print br0,br1,br2
	hoa0,hoa1,hoa2 = (randvecball()*(r-STRUTRAD) for r in (br0,br1,br2))
	hob0,hob1,hob2 = (randvecball()*(r-STRUTRAD) for r in (br0,br1,br2))
	hoa0,hoa1,hoa2 = bc0+hoa0,bc1+hoa1,bc2+hoa2
	hob0,hob1,hob2 = bc0+hob0,bc1+hob1,bc2+hob2
	hoa0 = closest_around(hoa0,c1,hob0,NFOLD[0],c0)
	hoa1 = closest_around(hoa1,c2,hob1,NFOLD[1],c1)
	hoa2 = closest_around(hoa2,c0,hob2,NFOLD[2],c2)
	hob0 = closest_around(hob0,c2,hoa0,NFOLD[0],c0)
	hob1 = closest_around(hob1,c0,hoa1,NFOLD[1],c1)
	hob2 = closest_around(hob2,c1,hoa2,NFOLD[2],c2)

	framecen = (c0+c1+c2).normalized()
	showsphere( bc0+framecen*MOVASYM                                           , br0 , col=(1,0,0)              , lbl="asym_jnt1" )
	showsphere( bc1+framecen*MOVASYM                                           , br1 , col=(0,1,0)              , lbl="asym_jnt2" )
	showsphere( bc2+framecen*MOVASYM                                           , br2 , col=(0,0,1)              , lbl="asym_jnt3" )
	if WHICHRODS[0]: showcyl( hoa0+framecen*MOVASYM, hob1+framecen*MOVASYM, STRUTRAD, col=(1,0,0), col2=(0,1,0), lbl="asym_rod1" )
	if WHICHRODS[1]: showcyl( hoa1+framecen*MOVASYM, hob2+framecen*MOVASYM, STRUTRAD, col=(0,1,0), col2=(0,0,1), lbl="asym_rod2" )
	if WHICHRODS[2]: showcyl( hoa2+framecen*MOVASYM, hob0+framecen*MOVASYM, STRUTRAD, col=(0,0,1), col2=(1,0,0), lbl="asym_rod3" )

	G0 = rotation_around_degrees(SYMAXIS[SYM,NFOLD[0]],360.0/NFOLD[0],SYMCEN[SYM,NFOLD[0]])
	G1 = rotation_around_degrees(SYMAXIS[SYM,NFOLD[1]],360.0/NFOLD[1],SYMCEN[SYM,NFOLD[1]])
	G2 = rotation_around_degrees(SYMAXIS[SYM,NFOLD[2]],360.0/NFOLD[2],SYMCEN[SYM,NFOLD[2]])
	# print G0.pretty(),c0
	# print G1.pretty(),c1
	# print G2.pretty(),c2
	for i,xs in enumerate(expand_xforms((G1,G2),SYMGEN)):
	# for i,xs in enumerate(SYMICS):
		# print xs.pretty()
		showsphere( xs*bc0 , 0.08 , col=(1,0,0) , lbl="sym%02i_jnt1"%i )
		showsphere( xs*bc1 , 0.08 , col=(0,1,0) , lbl="sym%02i_jnt2"%i )
		showsphere( xs*bc2 , 0.08 , col=(0,0,1) , lbl="sym%02i_jnt3"%i )

		if WHICHRODS[0]: showcyl( xs*hoa0 , xs*(c0)         , STRUTRAD, col=(1,0,0),lbl="sym%02i_cena1"%i)
		if WHICHRODS[1]: showcyl( xs*hoa1 , xs*(c1)         , STRUTRAD, col=(0,1,0),lbl="sym%02i_cena2"%i)
		if WHICHRODS[2]: showcyl( xs*hoa2 , xs*(c2)         , STRUTRAD, col=(0,0,1),lbl="sym%02i_cena3"%i)

		if WHICHRODS[2]: showcyl( xs*hob0 , xs*(c0)         , STRUTRAD, col=(1,0,0),lbl="sym%02i_cenb1"%i)
		if WHICHRODS[0]: showcyl( xs*hob1 , xs*(c1)         , STRUTRAD, col=(0,1,0),lbl="sym%02i_cenb2"%i)
		if WHICHRODS[1]: showcyl( xs*hob2 , xs*(c2)         , STRUTRAD, col=(0,0,1),lbl="sym%02i_cenb3"%i)

		if WHICHRODS[0] and WHICHRODS[2]: showcyl( xs*hoa0 , xs*hob0 , STRUTRAD, col=(1,0,0),lbl="sym%02i_cenc1"%i)
		if WHICHRODS[1] and WHICHRODS[0]: showcyl( xs*hoa1 , xs*hob1 , STRUTRAD, col=(0,1,0),lbl="sym%02i_cenc2"%i)
		if WHICHRODS[2] and WHICHRODS[1]: showcyl( xs*hoa2 , xs*hob2 , STRUTRAD, col=(0,0,1),lbl="sym%02i_cenc3"%i)

		if WHICHRODS[0]: showcyl( xs*hoa0 , xs*hob1 , STRUTRAD, col=(1,0,0),col2=(0,1,0),lbl="sym%02i_rod1"%i)
		if WHICHRODS[1]: showcyl( xs*hoa1 , xs*hob2 , STRUTRAD, col=(0,1,0),col2=(0,0,1),lbl="sym%02i_rod2"%i)
		if WHICHRODS[2]: showcyl( xs*hoa2 , xs*hob0 , STRUTRAD, col=(0,0,1),col2=(1,0,0),lbl="sym%02i_rod3"%i)


def makeimg():
	for i in range(10):
		print i
		cmd.delete("all")
		symgentest()
		cmd.ray()
		cmd.png("~/project/struts/ics3comp_%02i.png"%i)

def nulltest():
	"""
	>>> print "foo"
	foo
	"""

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite())
    return tests



if __name__ == '__main__':
   import doctest
   r = doctest.testmod()
   print r

