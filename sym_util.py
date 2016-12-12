# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import sys,os,inspect,functools
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path: sys.path.append(newpath)
import string,re,gzip,itertools
import collections
from pymol_util import *
import operator as op
from xyzMath import *
from itertools import product,ifilter
try:
	from cluster import HierarchicalClustering
except:
	print "couldn't import HierarchicalClustering from cluster"


nsymmetrizecx = 0




def get_xforms_by_chain(sele="all",verbose=False,userms=False):
	v = cmd.get_view()
	cen = com("("+sele+") and (name CA and not HET)")
	chains = cmd.get_chains(sele)
	xforms = dict()
	maxrms = 0.0
	for c1,c2 in filter( lambda t: True, product(chains,chains) ):
		refsele = "((%s) and chain %s and name CA and not HET)"%(sele,c1)
		movsele = "((%s) and chain %s and name CA and not HET)"%(sele,c2)
		if userms: x,rms = getrelframe_rmsalign( movsele, refsele, Xform(-cen) )
		else:      x,rms = getrelframe( movsele, refsele, Xform(-cen)), 0.0
		maxrms = max(maxrms,rms)
		xforms[c1,c2] = x
		#if c1 in "AB" and c2 in "AB":
		#	print movsele
		#	print refsele
		#	print x.pretty()
		#	print
	# raise Exception
	cmd.set_view(v)
	return xforms, maxrms

def find_symelems(sele_or_xforms="all",verbose=False):
	xforms = sele_or_xforms
	if isinstance(sele_or_xforms,basestring): xforms, maxrms = get_xforms_by_chain(sele_or_xforms,verbose=True)
	elif not isinstance(sele_or_xforms,dict): raise ValueError
	symelems = list()
	maxangerr = 0.0
	for c,x in xforms.items():
		assert len(c)==2
		assert isinstance(x,Xform)
		if c[0]==c[1]: continue
		dis = x.t.length()
		if dis > 5.0: continue
		axis,ang = x.rotation_axis()
		nfold = round(math.pi*2.0/ang)
		angerr = abs(ang-math.pi*2.0/nfold)*180.0/math.pi
		if verbose: print "candidate symelem:",nfold, c, angerr, axis
		if angerr > 360.0/nfold/8.0: continue # require unambiguous symelems
		maxangerr = max(maxangerr,angerr*nfold)
		symelems.append( (nfold,axis,c,angerr) )
	symelemdis = lambda x,y: line_line_angle_degrees(x[1],y[1]) if x[0]==y[0] else 9e9
	if verbose:
		for se1,se2 in filter( lambda t: t[0]<t[1], product(symelems,symelems) ):
			if se1[0]==se2[0]:
				print se1
				print se2
				print symelemdis(se1,se2), "degrees"
				print
	hier = HierarchicalClustering(symelems, symelemdis )
	thresh = 6.0
	clusters = hier.getlevel(thresh);
	print "number of symmetry element clusters at threshold",thresh,"degrees is" , len(clusters)
	centers0 = list()
	maxaxiserr = 0.0
	for clust in clusters:
		print "symelem cluster:",clust
		center = list(clust[0])
		center[2] = list((center[2],))
		for i in range(1,len(clust)):
			ax = clust[i][1]
			center[1] = center[1] + ( ax if ax.dot(center[1]) > 0 else -ax )
			center[2].append(clust[i][2])
			center[3] = max(center[3],clust[i][3])
		center[1].normalize()
		centers0.append(center)
		axiserr = 0.0
		for c in clust:	axiserr = max( axiserr, 1.0-abs(center[1].dot(c[1])) )
		maxaxiserr = max(maxaxiserr,axiserr)
	# sort on nfold, then on number of chain pairs in cluster
	centers0 = sorted( centers0, cmp = lambda x,y: cmp(y[0],x[0]) if x[0]!=y[0] else cmp(len(y[2]),len(x[2])) )
	centers = list()
	for center in centers0:
		if verbose: print "DEBUG prune center:",center
		seenit = False
		for censeen in centers:
			remainder = abs( ( censeen[0] / center[0] ) % 1.0)
			if verbose: print "   ",remainder,censeen
			if remainder > 0.01: continue # not a symmetry multiple
			if 1.0-abs(center[1].dot(censeen[1])) < 0.01:
				seenit = True # axis are same
		if not seenit:
			centers.append(center)
	print "centers:"
	cen_of_geom = com("("+sele_or_xforms+") and (name CA and not HET)")
	for center in centers:
		print center
		# if center[0]>2.1: continue
		#showvecfrompoint(50*center[1],cen_of_geom)
	return centers, maxrms, maxangerr, maxaxiserr

def guessdxaxes(sele="all",verbose=False):
	nfold = len(cmd.get_chains(sele))
	assert nfold % 2 is 0
	nfold /= 2
	symelems, maxrms, angerr, axiserr = find_symelems(sele,verbose=verbose)
	assert len(symelems) > 1
	assert symelems[0][0] == float(nfold)
	assert symelems[1][0] == float(2)
	axis_high = symelems[0][1]
	axis_low  = symelems[1][1]
	return axis_high, axis_low, maxrms, angerr, axiserr

def aligndx(sele='all',verbose=False):
	trans(sele,-com(sele+" and name CA and not HET"))
	haxis, laxis, maxrms, angerr, axiserr = guessdxaxes(sele,verbose=verbose)
	xalign = alignvectors( haxis, laxis, Uz, Ux )
	xform(sele,xalign)
	return maxrms, angerr, axiserr

def guesscxaxis(sele,nfold=None,chains0=list(),extrasel="name CA"):
	sele = "(("+sele+") and ("+extrasel+") and (not het))"
	check = False
	if not chains0:
		chains0.extend(cmd.get_chains(sele))
		check = True
	if not nfold:
		nfold = len(chains0)
		check = True
	# print chains0
	if check and len(chains0) != nfold:
		print chains0
		print "num chains != n-fold"
		return None
	print "chains0:", chains0
	chains = list()
	for i,c in enumerate(chains0):
		if isinstance(c,basestring):
			chains.append( (c,) )
		elif isinstance(c,collections.Iterable):
			chains.append( c )
		else:
			raise ValueError("chain must be string or list of strings")
	atoms = cmd.get_model(sele).atom
	chain_index = {}
	for i,clist in enumerate(chains):
		for c in clist:
			chain_index[c] = i
	coords = [list() for c in chains]
	print len(coords),[len(x) for x in coords]
	for a in atoms:
		if a.chain in chain_index:
			coords[chain_index[a.chain]].append(Vec(a.coord))
	for c in coords:
		print len(c)
	return cyclic_axis(coords)

def aligncx(sele,nfold,alignsele=None,tgtaxis=Uz,chains=list(),extrasel="name CA"):
	if not alignsele: alignsele = sele
	tmp = guesscxaxis(alignsele,nfold,chains,extrasel)
	if not tmp: return None
	axis,cen,diserr,angerr = tmp
	trans(sele,-cen)
	alignaxis(sele,tgtaxis,axis,xyz.Vec(0,0,0))
	return tmp

def align_helix( sele, nrepeat, tgt_axis=Vec(1,0,0) ):
	resi0 = int( cmd.get_model(sele).atom[0].resi )
	sel1 = "( %s ) and resi %i-%i" % ( sele, resi0+0*nrepeat, resi0+1*nrepeat-1 )
	sel2 = "( %s ) and resi %i-%i" % ( sele, resi0+1*nrepeat, resi0+2*nrepeat-1 )	
	xform = getrelframe_rmsalign( sel1, sel2 )[0]
	axis, ang, cen = xform.rotation_axis_center()
	print axis, ang, cen
	trans( sele, -cen )
	alignaxis( sele, tgt_axis, axis )
	trans( sele, Vec(0,0,-com(sele).z) )

# def alignd2(sele='all',chains=list()):
# 	alignsele = "(("+sele+") and (name CA))"
# 	if not chains: chains.extend(cmd.get_chains(alignsele))
# 	if 4 is not len(chains): raise NotImplementedError("D2 must have chains")
# 	ga1 = guesscxaxis( alignsele, 2,[ (chains[0],chains[1]), (chains[2],chains[3]) ] )
# 	ga2 = guesscxaxis( alignsele, 2,[ (chains[0],chains[2]), (chains[1],chains[3]) ] )
# 	assert ga1 is not None and ga2 is not None
# 	err = 90.0 - line_line_angle_degrees(ga1[0],ga2[0])
# 	x = alignvectors(ga1[0],ga2[0],Uz,Uy)
# 	xform(sele,x)
# 	trans(sele,-com(alignsele))
# 	return err

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



def showcxaxis(sele,nfold=None,chains=list(),length=30,col=(1,1,1),lbl="Cx Axis"):
	g = guesscxaxis(sele,nfold,chains)
	showvecfrompoint(g[0]*2*length,g[1]-g[0]*length,col=col,lbl=lbl)

def myint(s):
   i = len(s)
   while i > 0 and not s[:i].isdigit(): i -= 1
   if not i: return None
   return int(s[:i])


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

def rechain(sel,nres):
	chains = ROSETTA_CHAINS
	ntot = len(getres(sel))
	assert ntot % nres == 0
	for i in range(ntot/nres):
		cmd.alter("resi %i-%i"%( nres*i+1,nres*(i+1)),"chain='%s'"%chains[i])


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
		print o
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

def make_helix_old(sele='vis',n=30,nfold=1):
	n = n / nfold
	cmd.delete('helix')
	v = cmd.get_view()
	x0 = getrelframe(sele+' and chain B',sele+' and chain A')
	axis, ang, cen = x0.rotation_axis_center()
	# print axis, ang, cen
	cmd.create('tmp',sele+' and chain A and name n+ca+c')
	count = 0
	for nf in range(nfold):
		x = RAD( axis, nf * 360.0/nfold, cen )
		print nf, x.pretty()
		for i in range(n):
			cmd.create('Htmp%i'%count,'tmp')
			xform('Htmp%i'%count,x)
			cmd.alter('Htmp%i'%count,"chain='%s'"%ROSETTA_CHAINS[count])
			count += 1
			x = x * x0
	cmd.create("HELIX",'Htmp*')
	cmd.delete("Htmp*")
	cmd.delete('tmp')
	cmd.hide('ev','HELIX')
	cmd.show('lines','helix')
	util.cbc('HELIX')
	cmd.set_view(v)

cmd.extend('make_helix_old',make_helix_old)

def color_by_2component(col1="green",col2="cyan"):
	chains = r"""ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz!@#$&.<>?]{}|-_\~=%"""
	chains1 = [chains[i] for i in range(0,len(chains),2)]
	chains2 = [chains[i] for i in range(1,len(chains),2)]
	for c in chains1:
		print c+c+c+c+c
		cmd.color(col1,"chain "+c)
	for c in chains2:
		print c+c+c+c+c
		cmd.color(col2,"chain "+c)

def make_ab_components(dir):
	if not os.path.exists(dir+"_AB"):
		os.mkdir(dir+"_AB")
	for fn in os.listdir(dir):
		if not fn.endswith(".pdb") and not fn.endswith(".pdb.gz"): continue
		cmd.delete("all")
		cmd.load(dir+"/"+fn,"a")
		makec6("a",name="c6")
		cmd.save(dir+"_AB/"+fn,"c6 and chain A+B")

def trim_sym(sel='all',na=1,nb=1):
	a = [x[1] for x in getres(sel + " and chain A")]
	b = [x[1] for x in getres(sel + " and chain B")]
	for ia in range( len(a)/na, len(a) ):
		cmd.remove( sel + " and chain A and resi " + str(a[ia]) )
	for ib in range( len(b)/nb, len(b) ):
		cmd.remove( sel + " and chain B and resi " + str(b[ib]) )

cmd.extend('trim_sym',trim_sym)


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


