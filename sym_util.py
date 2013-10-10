# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import sys,os,inspect,functools
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path: sys.path.append(newpath)
import string,re,gzip,itertools
from pymol_util import *
import operator as op
from xyzMath import *

nsymmetrizecx = 0

def guesscxaxis(sele,nfold=None,chains=list(),extrasel="name CA"):
	sele = "(("+sele+") and ("+extrasel+"))"
	if not chains: chains.extend(cmd.get_chains(sele))
	if not nfold: nfold = len(chains)
	# print chains
	if len(chains) != nfold:
		print chains
		print "num chains != n-fold"
		return None
	atoms = cmd.get_model(sele).atom
	idx = {}
	for i,c in enumerate(chains): idx[c] = i
	coords = [list() for c in chains]
	print len(coords),[len(x) for x in coords]
	for a in atoms:
		if a.chain in idx:
			coords[idx[a.chain]].append(Vec(a.coord))
	return cyclic_axis(coords)

def aligncx(sele,nfold,alignsele=None,tgtaxis=Uz,chains=list(),extrasel="name CA"):
	if not alignsele: alignsele = sele
	tmp = guesscxaxis(alignsele,nfold,chains,extrasel)
	if not tmp: return None
	axis,cen,diserr,angerr = tmp
	trans(sele,-cen)
	alignaxis(sele,tgtaxis,axis,xyz.Vec(0,0,0))
	return tmp


def alignd2(sele='all',chains=list()):
	alignsele = "(("+sele+") and (name CA))"
	if not chains: chains.extend(cmd.get_chains(alignsele))
	if 4 is not len(chains): raise NotImplementedError("D2 must have chains")
	ga1 = guesscxaxis(alignsele,2,(chains[0],chains[1]))
	ga2 = guesscxaxis(alignsele,2,(chains[0],chains[2]))
	err = 90.0 - line_line_angle_degrees(ga1[0],ga2[0])
	x = alignvectors(ga1[0],ga2[0],Uz,Uy)
	xform(sele,x)
	trans(sele,-com(alignsele))
	return err

def symmetrize(sele="not symmetrized_*",alignsele=None,chains=list(),delete=True):
	global nsymmetrizecx
	if delete: cmd.delete("symmetrized_*")
	tmp = guesscxaxis(sele,None,chains)
	if not tmp: return None
	axis,cen,diserr,angerr = tmp
	# print "symmetrize TMP__C%i, distance err %f, algle error %f"%(len(chains),diserr,angerr)
	for i,c in enumerate(chains):
		newobj = "symmetrized_%i_%s"%(nsymmetrizecx,c)
		cmd.create(newobj,"(%s) and chain %s"%(sele,chains[0]))
		cmd.alter(newobj,"chain='%s'"%c)
		rot(newobj,axis,360.0*float(i)/float(len(chains)),cen)
		print "rot(",newobj,',',axis,',',360.0*float(i)/float(len(chains)),',',cen,')'
	newobj = "symmetrized_%i"%(nsymmetrizecx)
	cmd.create(newobj,"symmetrized_%i_*"%(nsymmetrizecx))
	cmd.delete("symmetrized_%i_*"%nsymmetrizecx)
	# print "rms",cmd.align(newobj,sele)
	nsymmetrizecx += 1
	return tmp

guessc2axis = functools.partial(guesscxaxis,nfold=2)
guessc3axis = functools.partial(guesscxaxis,nfold=3)
guessc4axis = functools.partial(guesscxaxis,nfold=4)
guessc5axis = functools.partial(guesscxaxis,nfold=5)
guessc6axis = functools.partial(guesscxaxis,nfold=6)
alignc2 = functools.partial(aligncx,nfold=2)
alignc3 = functools.partial(aligncx,nfold=3)
alignc4 = functools.partial(aligncx,nfold=4)
alignc5 = functools.partial(aligncx,nfold=5)
alignc6 = functools.partial(aligncx,nfold=6)

def makecx(sel = 'all',name="TMP",n = 5,axis=Uz):
	v = cmd.get_view()
	cmd.delete("TMP__C%i_*"%n)
	chains = ROSETTA_CHAINS
	for i in range(n): cmd.create("TMP__C%i_%i"%(n, i), sel+" and (not TMP__C%i_*)"%n)
	for i in range(n): rot("TMP__C%i_%i"%(n, i), axis, -360.0*float(i)/float(n))
	for i in range(n): cmd.alter("TMP__C%i_%i"%(n, i), "chain = '%s'"%chains[i])
	util.cbc("TMP__C*")
	# for i in range(n): cmd.alter("TMP__C%i_%i"%(n, i),"resi=str(int(resi)+%i)"%(1000*i));
	# util.cbc("TMP__C*")
	cmd.create(name,"TMP__*")
	cmd.delete("TMP__*")
	cmd.set_view(v)

def makedx(sel = 'all', n = 5):
	v = cmd.get_view()
	cmd.delete("D%i_*"%n)
	chains = ROSETTA_CHAINS
	for i in range(n):
		dsel  = "D%i_%i"%(n,  i)
		dsel2 = "D%i_%i"%(n,n+i)
		cmd.create(dsel , sel+" and (not D%i_*)"%n)
		rot       (dsel , Uz, 360.0*float(i)/float(n))
		cmd.create(dsel2, dsel )
		rot       (dsel2, Uy, 180.0 )
		cmd.alter (dsel , "chain = '%s'"%chains[i])
		cmd.alter (dsel2, "chain = '%s'"%chains[i+n])
	util.cbc("D*")
	cmd.set_view(v)

def maketet(sel='chain A+B and name n+ca+c',name="TET",n=12):
	v = cmd.get_view()
	cmd.delete("TMP__C%i_*"%n)
	cmd.delete(name)
	chains = ROSETTA_CHAINS
	for i in range(n): cmd.create("TMP__C%i_%i"%(n, i), sel+" and (not TMP__C%i_*)"%n)
	for i in range(n): xform("TMP__C%i_%i"%(n, i), SYMTET[i])
	for i in range(n): cmd.alter("TMP__C%i_%i"%(n, i),"resi=str(int(resi)+%i)"%(1000*i));
	util.cbc("TMP__C*")
	cmd.create(name,"TMP__*")
	cmd.delete("TMP__*")
	cmd.set_view(v)

def makeoct(sel='chain A+B and name n+ca+c',name="OCT",n=24):
	v = cmd.get_view()
	cmd.delete("TMP__C%i_*"%n)
	cmd.delete(name)
	chains = ROSETTA_CHAINS
	for i in range(n): cmd.create("TMP__C%i_%i"%(n, i), sel+" and (not TMP__C%i_*)"%n)
	for i in range(n): xform("TMP__C%i_%i"%(n, i), SYMOCT[i])
	for i in range(n): cmd.alter("TMP__C%i_%i"%(n, i),"resi=str(int(resi)+%i)"%(1000*i));
	util.cbc("TMP__C*")
	cmd.create(name,"TMP__*")
	cmd.delete("TMP__*")
	cmd.set_view(v)

def makeicos(sel='chain A+B and name n+ca+c and visible',name="ICOS",n=60):
	v = cmd.get_view()
	cmd.delete("TMP__C%i_*"%n)
	cmd.delete(name)
	chains = ROSETTA_CHAINS
	for i in range(n): cmd.create("TMP__C%i_%i"%(n, i), sel+" and (not TMP__C%i_*)"%n)
	for i in range(n): xform("TMP__C%i_%i"%(n, i), SYMICS[i])
	for i in range(n): cmd.alter("TMP__C%i_%i"%(n, i),"resi=str(int(resi)+%i)"%(1000*i));
	# util.cbc("TMP__C*")
	cmd.create(name,"TMP__*")
	cmd.delete("TMP__*")
	cmd.set_view(v)



def showcxaxis(sele,nfold=None,chains=list(),length=30,col=(1,1,1),lbl="Cx Axis"):
	g = guesscxaxis(sele,nfold,chains)
	showvecfrompoint(g[0]*2*length,g[1]-g[0]*length,col=col,lbl=lbl)

def myint(s):
   i = len(s)
   while i > 0 and not s[:i].isdigit(): i -= 1
   if not i: return None
   return int(s[:i])

def mki213(N, sel = 'all'):
	v = cmd.get_view()
	cmd.delete("i213_*")
	cmd.delete('base80345769083457')
	cmd.delete('tmp80345769083457')
	c2 = com(sel)
	c3 = xyz.Vec(0, 0, 0)
	cmd.create( 'tmp80345769083457', sel)
	a2 = c2axis('tmp80345769083457')
	cmd.delete( 'tmp80345769083457')
	a3 = xyz.Vec(0, 0, 1)
	cmd.create('base80345769083457', sel+" and chain A and visible")
	seenit = []
	R2 = [xyz.rotation_matrix_degrees(a2, 0), xyz.rotation_matrix_degrees(a2, 180), ]
	R3 = [xyz.rotation_matrix_degrees(a3, 0), xyz.rotation_matrix_degrees(a3, 120), xyz.rotation_matrix_degrees(a3, 240), ]
	C = alphabet
	print a2, c2, a3, c3
	for i21 in range(2):
		for i32 in range(3 if N > 1 else 1):
			for i22 in range(2 if N > 2 else 1):
				for i33 in range(3 if N > 3 else 1):
					for i23 in range(2 if N > 4 else 1):
						for i34 in range(3 if N > 5 else 1):
							for i24 in range(2 if N > 6 else 1):
								for i35 in range(3 if N > 7 else 1):
									for i25 in range(2 if N > 8 else 1):
										test = xyz.Vec(0, 0, 0)
										test = R2[i21]*(test-c2)+c2
										test = R3[i32]*(test-c3)+c3
										test = R2[i22]*(test-c2)+c2
										test = R3[i33]*(test-c3)+c3
										test = R2[i23]*(test-c2)+c2
										test = R3[i34]*(test-c3)+c3
										test = R2[i24]*(test-c2)+c2
										test = R3[i35]*(test-c3)+c3
										test = R2[i25]*(test-c2)+c2
										#print test
										seen = False
										for xs in seenit:
											if (xs-test).length() < 0.1:
												seen = True
												break
										if seen: continue
										else: seenit.append(test)
										n = "i213_%i%i%i%i%i%i%i%i%i"%(i25, i35, i24, i34, i23, i33, i22, i32, i21)
										cmd.create(n, 'base80345769083457')
										rot(n, a2, i21*180.0, c2)
										rot(n, a3, i32*120.0, c3)
										rot(n, a2, i22*180.0, c2)
										rot(n, a3, i33*120.0, c3)
										rot(n, a2, i23*180.0, c2)
										rot(n, a3, i34*120.0, c3)
										rot(n, a2, i24*180.0, c2)
										rot(n, a3, i35*120.0, c3)
										rot(n, a2, i25*180.0, c2)
	print len(seenit)
	cmd.delete('base80345769083457')
	cmd.set_view(v)

def viewi213(sel = "all"):
	cmd.hide('ev')
	cmd.show('rib')
	mki213(sel)
	cmd.show('car', 'not i213*')
	cmd.hide('rib', 'not i213*')
	cmd.show('lines', '(byres (%s and not i213* and chain A) within 7.0 of (%s and not i213* and chain B))'%(sel, sel))
	cmd.show('lines', '(byres (%s and not i213* and chain B) within 7.0 of (%s and not i213* and chain A))'%(sel, sel))


def mkp23(N, R=43.5, i=0, sel = 'all'):
	v = cmd.get_view()
	cmd.delete("p23_*")
	cmd.delete('base80345769083457')
	cmd.delete('tmp80345769083457')
	c2 = xyz.Vec(0, 0, 0)
	c3 = xyz.Vec(R,R,-R)
	cmd.create( 'tmp80345769083457', sel)
	cmd.delete( 'tmp80345769083457')
	a3 = [xyz.Vec(0,0,0),xyz.Vec(1,1,1),xyz.Vec(-1,-1,-1)]
	a2 = [xyz.Vec(0,0,0),xyz.Vec(1,0,0),xyz.Vec(0,1,0),xyz.Vec(0,0,1)]
	cmd.create('base80345769083457', sel+" and visible")
	seenit = []
	R2 = [xyz.rotation_matrix_degrees(a2[1],  0), # hack
		  xyz.rotation_matrix_degrees(a2[1],180),
		  xyz.rotation_matrix_degrees(a2[2],180),
		  xyz.rotation_matrix_degrees(a2[3],180) ]
	R3 = [xyz.rotation_matrix_degrees(a3[1],  0), # hack!
		  xyz.rotation_matrix_degrees(a3[1],120),
		  xyz.rotation_matrix_degrees(a3[2],120), ]
	C = alphabet
	print a2, c2, a3, c3
	for i21 in range(4):
		for i32 in range(3 if N > 1 else 1):
			for i22 in range(4 if N > 2 else 1):
				for i33 in range(3 if N > 3 else 1):
					for i23 in range(4 if N > 4 else 1):
						for i34 in range(3 if N > 5 else 1):
							for i24 in range(4 if N > 6 else 1):
								for i35 in range(3 if N > 7 else 1):
									for i25 in range(4 if N > 8 else 1):
										test = xyz.Vec(1, 1, 1)
										test = R2[i21]*(test-c2)+c2
										test = R3[i32]*(test-c3)+c3
										test = R2[i22]*(test-c2)+c2
										test = R3[i33]*(test-c3)+c3
										test = R2[i23]*(test-c2)+c2
										test = R3[i34]*(test-c3)+c3
										test = R2[i24]*(test-c2)+c2
										test = R3[i35]*(test-c3)+c3
										test = R2[i25]*(test-c2)+c2
										#print test
										seen = False
										for xs in seenit:
											if (xs-test).length() < 0.1:
												seen = True
												break
										if seen: continue
										else: seenit.append(test)
										n = "p23_%i%i%i%i%i%i%i%i%i"%(i25, i35, i24, i34, i23, i33, i22, i32, i21)
										cmd.create(n, 'base80345769083457')
										if i21 > 0: rot(n, a2[i21], 180.0, c2)
										if i32 > 0: rot(n, a3[i32], 120.0, c3)
										if i22 > 0: rot(n, a2[i22], 180.0, c2)
										if i33 > 0: rot(n, a3[i33], 120.0, c3)
										if i23 > 0: rot(n, a2[i23], 180.0, c2)
										if i34 > 0: rot(n, a3[i34], 120.0, c3)
										if i24 > 0: rot(n, a2[i24], 180.0, c2)
										if i35 > 0: rot(n, a3[i35], 120.0, c3)
										if i25 > 0: rot(n, a2[i25], 180.0, c2)
	print "seen:",len(seenit)
	cmd.delete('base80345769083457')
	cmd.set_view(v)

def selbycomp(trn=0):
	cmd.select("TRI1","TRI and chain A+B+C")
	cmd.select("TRI2","TRI and chain D+E+F")
	cmd.select("TRI3","TRI and chain G+H+I")
	cmd.select("TRI4","TRI and chain J+K+L")
	cmd.select("TRI5","TRI and chain xyz.Mat+N+O")
	cmd.select("TRI6","TRI and chain P+Q+R")
	cmd.select("TRI7","TRI and chain S+T+U")
	cmd.select("TRI8","TRI and chain xyz.Vec+W+Ux")
	cmd.select("DIM1","DIM and chain A+D")
	cmd.select("DIM2","DIM and chain B+G")
	cmd.select("DIM3","DIM and chain C+J")
	cmd.select("DIM4","DIM and chain E+U")
	cmd.select("DIM5","DIM and chain F+R")
	cmd.select("DIM6","DIM and chain H+T")
	cmd.select("DIM7","DIM and chain I+O")
	cmd.select("DIM8","DIM and chain K+Q")
	cmd.select("DIM9","DIM and chain L+N")
	cmd.select("DIM10","DIM and chain xyz.Mat+xyz.Vec")
	cmd.select("DIM11","DIM and chain P+W")
	cmd.select("DIM12","DIM and chain Ux+S")
	cmd.delete("LINE*")
	cmd.delete("serf*")

	cmd.do("""alter all, b=50
	alter all, q=1
	set gaussian_resolution,8""")
	ISO="""map_new map%s, gaussian, 2, %s, 10
	isosurface surf%s, map%s"""


	for i in range(1, 9):
		cmd.do(ISO%(("TRI%i"%i,)*4))
		cmd.color(COLORS[i-1],"surfTRI%i"%i)
		c = com("TRI%i"%i)
		# trans("TRI%i"%i,trn*c.normalized())
		obj = [
			cgo.CYLINDER,
			0.0, 0.0, 0.0,
			1.6*c.x, 1.6*c.y, 1.6*c.z,
			1.5,
			0.1,0.1,0.1,0.1,0.1,0.1,
		]
		cmd.load_cgo(obj,'LINETRI%i'%i)
	for i in range(1,13):
		cmd.do(ISO%(("DIM%i"%i,)*4))
		cmd.color(COLORS[i+7],"surfDIM%i"%i)
		c = com("DIM%i"%i)
		# trans("DIM%i"%i,trn*com("DIM%i"%i).normalized())
		obj = [
			cgo.CYLINDER,
			0.0, 0.0, 0.0,
			1.3*c.x, 1.3*c.y, 1.3*c.z,
			1.0,
			0,0,1,0,0,1
		]
		cmd.load_cgo(obj,'LINEDIM%i'%i)

def getframe(obj):
	m = cmd.get_model(obj)
	x = xyz.Vec(m.atom[       0     ].coord)
	y = xyz.Vec(m.atom[len(m.atom)/2].coord)
	z = xyz.Vec(m.atom[      -1     ].coord)
	return xyz.stub(x,y,z)

def getrelframe(newobj,refobj,Forigin=None):
	"assume the obj's are identical"
	if Forigin is None: Forigin = xyz.Xform(xyz.Imat,xyz.Vec(0,0,0))
	Fref = getframe(refobj)
	Fnew = getframe(newobj)
	Fdelta = Fnew * ~Fref
	#print (Fdelta * Forigin).pretty()
	return Fdelta * Forigin



def rechain(sel,nres):
	chains = ROSETTA_CHAINS
	ntot = len(getres(sel))
	assert ntot % nres == 0
	for i in range(ntot/nres):
		cmd.alter("resi %i-%i"%( nres*i+1,nres*(i+1)),"chain='%s'"%chains[i])

for i in range(2,21):
		globals()['makec%i'%i] = partial(makecx,n=i)

for i in range(2,21):
		globals()['maked%i'%i] = partial(makedx,n=i)

def makecxauto():
	for o in cmd.get_object_list():
		n = int(re.search("_C\d+_", o).group(0)[2:-1])
		makecx(o, n)

def makekinwire(sel,movres,fixres):
	v = cmd.get_view()
	cmd.delete("ha"); cmd.create("ha",sel); cmd.alter("ha","chain='A'")
	cmd.delete("hb"); cmd.create("hb",sel); cmd.alter("hb","chain='B'")
	cmd.delete("hc"); cmd.create("hc",sel); cmd.alter("hc","chain='C'")
	cmd.delete("hd"); cmd.create("hd",sel); cmd.alter("hd","chain='D'")
	cmd.delete("he"); cmd.create("he",sel); cmd.alter("he","chain='E'")
	cmd.align("hb and resi %i"%movres,"ha and resi %i"%fixres);
	cmd.align("hc and resi %i"%movres,"hb and resi %i"%fixres);
	cmd.align("hd and resi %i"%movres,"hc and resi %i"%fixres);
	cmd.align("he and resi %i"%movres,"hd and resi %i"%fixres);
	util.cbc('elem C')
	v = cmd.set_view(v)

def get_contigs(x,n=7):
	"""
	>>> test = list(range(1,8)) + list(range(20,33)) + list(range(40,44)) + list(range(49,50))+ list(range(0,8))
	>>> print test
	[1, 2, 3, 4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 40, 41, 42, 43, 49, 0, 1, 2, 3, 4, 5, 6, 7]

	>>> print get_contigs( test )
	[[1, 2, 3, 4, 5, 6, 7], [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], [0, 1, 2, 3, 4, 5, 6, 7]]
	"""
	if type(x) is type(""):
		x = [x[1] for x in getres(x)]
	x.append(-123456789)
	contigs = [[],]
	for i in range(len(x)-1):
		contigs[-1].append(x[i])
		if x[i]+1 is not x[i+1]:
			contigs.append(list())
	return [c for c in contigs if len(c) >= n]

# def get_contigs_termini(x,n=7):
# 	"""
# 	>>> test = list(range(1,8)) + list(range(20,33)) + list(range(40,44)) + list(range(49,50))+ list(range(0,8))
# 	>>> print test
# 	[1, 2, 3, 4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 40, 41, 42, 43, 49, 0, 1, 2, 3, 4, 5, 6, 7]

# 	>>> print get_contigs_termini( test )
# 	[[1, 2, 3, 4, 5, 6, 7], [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32], [0, 1, 2, 3, 4, 5, 6, 7]]
# 	"""
# 	pass #cend = []

def get_fixed_size_contigs(x,n=7):
	"""
	>>> test = list(range(1,8)) + list(range(20,33)) + list(range(40,44)) + list(range(49,50))+ list(range(0,8))
	>>> print test
	[1, 2, 3, 4, 5, 6, 7, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 40, 41, 42, 43, 49, 0, 1, 2, 3, 4, 5, 6, 7]

	>>> for f in get_fixed_size_contigs(test,7): print f
	[1, 2, 3, 4, 5, 6, 7]
	[20, 21, 22, 23, 24, 25, 26]
	[21, 22, 23, 24, 25, 26, 27]
	[22, 23, 24, 25, 26, 27, 28]
	[23, 24, 25, 26, 27, 28, 29]
	[24, 25, 26, 27, 28, 29, 30]
	[25, 26, 27, 28, 29, 30, 31]
	[26, 27, 28, 29, 30, 31, 32]
	[0, 1, 2, 3, 4, 5, 6]
	[1, 2, 3, 4, 5, 6, 7]

	>>> for f in get_fixed_size_contigs(test,9): print f
	[20, 21, 22, 23, 24, 25, 26, 27, 28]
	[21, 22, 23, 24, 25, 26, 27, 28, 29]
	[22, 23, 24, 25, 26, 27, 28, 29, 30]
	[23, 24, 25, 26, 27, 28, 29, 30, 31]
	[24, 25, 26, 27, 28, 29, 30, 31, 32]

	>>> print len(get_fixed_size_contigs(test,1))
	28

	>>> for f in get_fixed_size_contigs(test,4): print f
	[1, 2, 3, 4]
	[2, 3, 4, 5]
	[3, 4, 5, 6]
	[4, 5, 6, 7]
	[20, 21, 22, 23]
	[21, 22, 23, 24]
	[22, 23, 24, 25]
	[23, 24, 25, 26]
	[24, 25, 26, 27]
	[25, 26, 27, 28]
	[26, 27, 28, 29]
	[27, 28, 29, 30]
	[28, 29, 30, 31]
	[29, 30, 31, 32]
	[0, 1, 2, 3]
	[1, 2, 3, 4]
	[2, 3, 4, 5]
	[3, 4, 5, 6]
	[4, 5, 6, 7]
	"""
	f = []
	for c in get_contigs(x):
		for i in range(0,len(c)-n+1):
			f.append(range(c[i],c[i]+n))
	return f

def tmpname():
	return "TEMPORARY_"+str(random.random())

def gen_helical_alignments(sele1,sele2,pref="HALN"):
	cmd.delete(pref+"_*")
	cmd.alter('(%s) or (%s)'%(sele1,sele2),'resn="ALA"')
	cmd.alter(sele2+' and chain A','chain="Z"')
	cmd.alter(sele2+' and chain B','chain="Y"')
	cmd.alter(sele2+' and chain C','chain="X"')
	cmd.alter(sele2+' and chain D','chain="W"')
	cmd.alter(sele2+' and chain E','chain="V"')
	cmd.alter(sele2+' and chain F','chain="U"')
	cmd.alter(sele2+' and chain G','chain="T"')
	cmd.alter(sele2+' and chain H','chain="S"')
	cmd.alter(sele2+' and chain I','chain="R"')
	cmd.alter(sele2+' and chain J','chain="Q"')
	chunks1 = ["chain A and resi "+str(h)[1:-1].replace(', ','+') for h in get_fixed_size_contigs("chain A and %s and ss H"%sele1)]
	chunks2 = ["chain Z and resi "+str(h)[1:-1].replace(', ','+') for h in get_fixed_size_contigs("chain Z and %s and ss H"%sele2)]
	for i,hsel1 in enumerate(chunks1):
		name1 = pref+"_"+tmpname()
		algn1 = name1+" and "+hsel1
		cmd.create(name1,sele1)
		for j,hsel2 in enumerate(chunks2):
			name2 = pref+"_"+tmpname()
			algn2 = name2+" and "+hsel2
			cmd.create(name2,sele2)
			print algn2,algn1
			print name1+" and chain A and "+hsel1
			print name2+" and chain Z and "+hsel2
			# now align them
			cmd.align( name2+" and chain Z and "+hsel2, name1+" and chain A and "+hsel1 )
			name3 = pref+"_%03i_%03i"%(i,j)
			cmd.create(name3,name1+" or "+name2)
			util.cbc(name3+" and elem C")
			cmd.delete(name2)
		cmd.delete(name1)

def colorI53(sel="visible"):
	a1 = '('+sel+') and (elem C) and (chain A+C+E+G+I)'
	a2 = '('+sel+') and (elem C) and (chain 2+M+Q+U+Y)'
	b1 = '('+sel+') and (elem C) and (chain B+L+N)'
	b2 = '('+sel+') and (elem C) and (chain 3+D+b)'
	cmd.color('cyan'    ,a1)
	cmd.color('green'   ,b1)
	cmd.color('magenta' ,a2)
	cmd.color('yellow'  ,b2)


def alignsym(sel="all",arch="I32",ax1=Vec(0,0,1),ax2=Vec(0.356825,0.000002,0.934171)):
	if arch == "T32":
		tgt1 = Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()
		tgt2 = Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()
	if arch == "T33":
		tgt1 = Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()
		tgt2 = Vec( 1.00000000000000, 1.00000000000000,-1.00000000000000 ).normalized()
	if arch == "O32":
		tgt1 = Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()
		tgt2 = Vec( 1.00000000000000, 1.00000000000000, 0.00000000000000 ).normalized()
	if arch == "O42":
		tgt1 = Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()
		tgt2 = Vec( 1.00000000000000, 1.00000000000000, 0.00000000000000 ).normalized()
	if arch == "O43":
		tgt1 = Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()
		tgt2 = Vec( 1.00000000000000, 1.00000000000000, 1.00000000000000 ).normalized()
	if arch == "I32":
		tgt1 = Vec( 0.93417235896272, 0.00000000000000, 0.35682208977309 ).normalized()
		tgt2 = Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()
	if arch == "I52":
		tgt1 = Vec( 0.85065080835204, 0.52573111211914, 0.00000000000000 ).normalized()
		tgt2 = Vec( 1.00000000000000, 0.00000000000000, 0.00000000000000 ).normalized()
	if arch == "I53":
		tgt1 = Vec( 0.85065080835204, 0.52573111211914, 0.00000000000000 ).normalized()
		tgt2 = Vec( 0.93417235896272, 0.00000000000000, 0.35682208977309 ).normalized()
	if abs(ax1.angle(ax2) - tgt1.angle(tgt2) ) > 0.001:
		print "your axes aren't spaced correctly for",arch,"angle should be",tgt1.angle(tgt2)
		return
	x = alignvectors( ax1, ax2, tgt1, tgt2 )
	xform(sel,x)

def xtal_frames(tgt=None,skip=tuple(),r=100):
	axes = list()
	objs = cmd.get_object_list()
	if not tgt: tgt = objs[0]
	assert tgt in objs
	c = com(tgt)
	covered = list()
	for o in objs:
		if o == tgt: continue
		x = getrelframe(tgt,o)

		# seenit = False
		# for xc in covered:
		# 	if x.t.distance(xc.t) < 0.001:
		# 		seenit = True
		# if seenit: continue

		axis,ang = x.rotation_axis()
		if ang < 1.0: continue # hack, nfold <= 6
		mov = proj( axis, x.t).length()
		if abs(mov) > 0.001: continue
		nf = 2*math.pi/ang
		if nf % 1.0 > 0.001: continue
		nf = int(nf)
		if nf in skip: continue
		if nf > 6 or nf == 5 or nf == 1: continue
		ctot = Vec(0,0,0)
		xtot = Xform()
		for i in range(nf):
			ctot += c
			covered.append(xtot)
			c = x * c
			xtot *= x
		ctot /= nf
		# beg = ctot - r*axis
		# end = ctot + r*axis
		beg = ray_sphere_intersection( axis,ctot,c,r)
		end = ray_sphere_intersection(-axis,ctot,c,r)
		if not beg or not end: continue
		if nf is not 2: showcyl(beg,end,0.3,col=(1.0,1.0,1.0))
		else:           showcyl(beg,end,0.2,col=(1.0,0.5,0.2))
		print round(nf),ctot,o

def makeh(sele='vis',n=30):
	cmd.delete('helix')
	v = cmd.get_view()
	x0 = getrelframe(sele+' and chain B',sele+' and chain A')
	x = Xform()
	cmd.create('tmp',sele+' and chain A and name n+ca+c')
	for i in range(n):
		cmd.create('Htmp%i'%i,'tmp')
		xform('Htmp%i'%i,x)
		cmd.alter('Htmp%i'%i,"chain='%s'"%ROSETTA_CHAINS[i])
		print ROSETTA_CHAINS[i]
		x = x * x0
	cmd.create("HELIX",'Htmp*')
	cmd.delete("Htmp*")
	cmd.delete('tmp')
	cmd.hide('ev','HELIX')
	cmd.show('lines','helix')
	util.cbc('HELIX')
	cmd.set_view(v)

cmd.extend('makeh',makeh)

def make_ab_components(dir):
	if not os.path.exists(dir+"_AB"):
		os.mkdir(dir+"_AB")
	for fn in os.listdir(dir):
		if not fn.endswith(".pdb") and not fn.endswith(".pdb.gz"): continue
		cmd.delete("all")
		cmd.load(dir+"/"+fn,"a")
		makec6("a",name="c6")
		cmd.save(dir+"_AB/"+fn,"c6 and chain A+B")

def nulltest():
	"""
	>>> print "foo"
	foo
	"""
	return None

def load_tests(loader, tests, ignore):
	tests.addTests(doctest.DocTestSuite())
	return tests

if __name__ == '__main__':
   import doctest
   r = doctest.testmod()
   print r


