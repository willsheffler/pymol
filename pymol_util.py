# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import sys,os,inspect
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path:
		sys.path.append(newpath)
import string
import math
import re
import gzip
import itertools
import glob
import sets
from pymol import cmd
from pymol.cgo import *
from random import randrange
from math import sqrt
import xyzMath as xyz
from xyzMath import Ux,Uy,Uz,Imat
from functools import partial

numcom = 0
numvec = 0

COLORS = ('cyan', 'lightmagenta', 'yellow', 'salmon', 'hydrogen', 'slate', 'orange',
 'lime',
 'deepteal',
 'hotpink',
 'yelloworange',
 'violetpurple',
 'grey70',
 'marine',
 'olive',
 'smudge',
 'teal',
 'dirtyviolet',
 'wheat',
 'deepsalmon',
 'lightpink',
 'aquamarine',
 'paleyellow',
 'limegreen',
 'skyblue',
 'warmpink',
 'limon',
 'violet',
 'bluewhite',
 'greencyan',
 'sand',
 'forest',
 'lightteal',
 'darksalmon',
 'splitpea',
 'raspberry',
 'grey50',
 'deepblue',
 'brown')
COLORS += COLORS + COLORS + COLORS

aa_1_3 = {'A': 'ALA',
		  'C': 'CYS',
		  'D': 'ASP',
		  'E': 'GLU',
		  'F': 'PHE',
		  'G': 'GLY',
		  'H': 'HIS',
		  'I': 'ILE',
		  'K': 'LYS',
		  'L': 'LEU',
		  'xyz.Mat': 'MET',
		  'N': 'ASN',
		  'P': 'PRO',
		  'R': 'ARG',
		  'Q': 'GLN',
		  'S': 'SER',
		  'T': 'THR',
		  'xyz.Vec': 'VAL',
		  'W': 'TRP',
		  'Uy': 'TYR',
}
aa_3_1 = {'ALA': 'A',
		  'ASP': 'D',
		  'GLU': 'E',
		  'CYS': 'C',
		  'PHE': 'F',
		  'GLY': 'G',
		  'HIS': 'H',
		  'ILE': 'I',
		  'LYS': 'K',
		  'LEU': 'L',
		  'MET': 'xyz.Mat',
		  'ASN': 'N',
		  'PRO': 'P',
		  'GLN': 'Q',
		  'ARG': 'R',
		  'SER': 'S',
		  'THR': 'T',
		  'VAL': 'xyz.Vec',
		  'TRP': 'W',
		  'TYR': 'Uy',
}
aa_types = {
  'A': 'hydrophobic',
  'C': 'cysteine',
  'D': 'negative',
  'E': 'negative',
  'F': 'aromatic',
  'G': 'hydrophobic',
  'H': 'polar',
  'I': 'hydrophobic',
  'K': 'positive',
  'L': 'hydrophobic',
  'xyz.Mat': 'hydrophobic',
  'N': 'polar',
  'P': 'proline',
  'Q': 'polar',
  'R': 'positive',
  'S': 'polar',
  'T': 'polar',
  'xyz.Vec': 'phydrophobic',
  'W': 'aromatic',
  'Uy': 'aromatic',
}

def showaxes():
    v = cmd.get_view()
    obj = [
        BEGIN, LINES,
        COLOR, 1.0, 0.0, 0.0, 
        VERTEX,   0.0, 0.0, 0.0,
        VERTEX,  20.0, 0.0, 0.0,
        COLOR, 0.0, 1.0, 0.0, 
        VERTEX, 0.0,   0.0, 0.0,
        VERTEX, 0.0,  20.0, 0.0, 
        COLOR, 0.0, 0.0, 1.0, 
        VERTEX, 0.0, 0.0,   0.0,
        VERTEX, 00, 0.0,   20.0,
        END
        ]                                                                                            
    cmd.load_cgo(obj,'axes')
    cmd.set_view(v)

def showaxes2():
    v = cmd.get_view()
    obj = [
        BEGIN, LINES,
        COLOR, 1.0, 0.0, 0.0, 
        VERTEX, -20.0, 0.0, 0.0,
        VERTEX,  20.0, 0.0, 0.0,
        COLOR, 0.0, 1.0, 0.0, 
        VERTEX, 0.0, -20.0, 0.0,
        VERTEX, 0.0,  20.0, 0.0, 
        COLOR, 0.0, 0.0, 1.0, 
        VERTEX, 0.0, 0.0, -20.0,
        VERTEX, 00, 0.0,   20.0,
        END
        ]                                                                                            
    cmd.load_cgo(obj,'axes')
    cmd.set_view(v)

def getchain(sele):
	try:
		c = list(sets.Set([x.chain for x in cmd.get_model(sele).atom]))
		c.sort()
		return c
	except:
		return []


def getres(sele):
	try:
		r = list(sets.Set([(x.chain, int(x.resi)) for x in cmd.get_model(sele).atom]))
		r.sort()
		return r
	except:
		return []


def com(sel="all", state=1):
	## assumes equal weights (best called with "and name ca" suffix)
	model = cmd.get_model(sel, state)
	c = xyz.Vec(0)
	for a in model.atom:
		c += xyz.Vec(a.coord)
	c = c / len(model.atom)
	return c


def showcom(sel="all"):
	global numcom
	c = com(sel)
	print "Center of mass: ", c
	cgo = [pymol.cgo.COLOR, 1.0, 1.0, 1.0, SPHERE, c.x, c.y, c.z, 1.0] ## white sphere with 3A radius
	cmd.load_cgo(cgo, "com%i" % numcom)
	numcom += 1


def showvec(c,lab=None):
	if lab is None:
		global numvec
		lab = "vec%i"%numvec
		numvec += 1
	cgo = [COLOR, 1.0, 1.0, 1.0, SPHERE, c.x, c.y, c.z, 1.0] ## white sphere with 3A radius
	cmd.load_cgo(cgo,lab)


def showvecfrompoint(a, c, lbl):
	cmd.load_cgo([COLOR, 1.0, 1.0, 1.0,
				  SPHERE,   c.x,       c.y,       c.z,    0.2,
				  CYLINDER, c.x,       c.y,       c.z,
							c.x + a.x, c.y + a.y, c.z + a.z, 0.1,
					1, 1, 1, 1, 1, 1, ], lbl)


class ResBB(object):
	"""docstring for ResBB"""
	def __init__(self, n, ca = None, c = None, ss = None):
		super(ResBB, self).__init__()
		if type(n) is ResBB:
			self.n  = n.n
			self.ca = n.ca
			self.c  = n.c
			self.ss = n.ss
			return
		elif type(n) is type(""):
			assert len(getres(n)) is 1
			m = cmd.get_model(n)
			self.n  = xyz.Vec( m.atom[0].coord )
			self.ca = xyz.Vec( m.atom[1].coord )
			self.c  = xyz.Vec( m.atom[2].coord )
			self.ss = m.atom[1].ss
			return
		assert type(n) is xyz.Vec
		assert type(ca) is xyz.Vec
		assert type(c) is xyz.Vec
		self.n  = n
		self.ca = ca
		self.c  = c
		self.ss = ss
	def rms(s, o):
		assert type(o) is ResBB
		return math.sqrt( (s.n-o.n)*(s.n-o.n) + (s.ca-o.ca)*(s.ca-o.ca) + (s.c-o.c)*(s.c-o.c) )
	def __rmul__(r, m):
		return ResBB( m*r.n, m*r.ca, m*r.c, r.ss )
	def __add__(r, v):
		return ResBB( r.n+v, r.ca+v, r.c+v, r.ss )
	def __sub__(r, v):
		return r + (-v)
	def __str__(r):
		return "n "+str(r.n)+", ca "+str(r.ca)+", c "+str(r.c)+", ss "+r.ss
	def stub(r):
		return xyz.Xform.from_four_points(xyz.Xform(None, None), r.ca, r.n, r.ca, r.c)


class DisulfLib(object):
	"""docstring for DisulfLib"""
	def __init__(self, fn):
		super(DisulfLib, self).__init__()
		self.fn = fn
		self.lib = {"EE":[], "EH":[], "EL":[], "HE":[], "HH":[], "HL":[], "LE":[], "LH":[], "LL":[]}
		for l in open(fn).readlines():
			ss1, ss2, sep, rt, xx, xy, xz, yx, yy, yz, zx, zy, zz, x, y, z = l.split()
			self.lib[ss1+ss2].append(Jump(xyz.Mat(xx, xy, xz, yx, yy, yz, zx, zy, zz), xyz.Vec(x, y, z)))

	def disulf_rms(self, r1, r2):
		assert type(r1) is ResBB
		assert type(r2) is ResBB
		minrms = 9e9
		if (r1.ca-r2.ca).length() > 10: return minrms
		ss1, ss2 = r1.ss, r2.ss
		if r1.ss not in "HEL" or len(r1.ss) != 1: ss1 = "L"
		if r2.ss not in "HEL" or len(r2.ss) != 1: ss2 = "L"
		for j in self.lib[ss1+ss2][:1]:
			s = r1.stub()
			r1_in_s = s.to_frame(r1)
			r3_in_s = r1_in_s + xyz.Vec(1, 1, 1) #( j.t * r1_in_s ) + j.trans
			r3      = s.from_frame(r3_in_s)
			rms = r2.rms(r3)
			if rms < minrms: minrms = rms
		return minrms


def chirality(fe1, fe2, fe3, fe4):
	a = fe2-fe1
	b = fe3-fe1
	c = a.cross(b)
	return c.dot(fe4-fe1)


def dimeraxis1(m1, m2):
	if len(m1.atom) != len(m2.atom):
		print "selections must contain the same number of atoms!"
		return
	i = j = randrange(len(m1.atom))
	while j == i: j = randrange(len(m1.atom))
	a1, a2 = xyz.Vec(m1.atom[i].coord), xyz.Vec(m2.atom[i].coord)
	u = (a1+a2)/3
	a1, a2 = xyz.Vec(m1.atom[j].coord), xyz.Vec(m2.atom[j].coord)
	v = (a1+a2)/3
	if v.x-u.x < 0:
		u, v = v, u
	return (v-u).normalized()


def dimeraxis(sel1, sel2, nsamp=100):
	m1 = cmd.get_model(sel1)
	m2 = cmd.get_model(sel2)
	a, wtot = xyz.Vec(0), 0
	for i in range(int(nsamp)):
		u = dimeraxis1(m1, m2)
		# print u
		a  += u
		wtot += u.length()
	a /= wtot
	return a


def axisofrot1(m1, m2, m3):
	if len(m1.atom) != len(m2.atom) or len(m2.atom) != len(m3.atom) or len(m3.atom) != len(m1.atom):
		print "selections must contain the same number of atoms!"
		return
	i = randrange(len(m1.atom))
	a1, a2, a3 = xyz.Vec(m1.atom[i].coord), xyz.Vec(m2.atom[i].coord), xyz.Vec(m3.atom[i].coord)
	u = (a1+a2+a3)/3
	i = randrange(len(m1.atom))
	a1, a2, a3 = xyz.Vec(m1.atom[i].coord), xyz.Vec(m2.atom[i].coord), xyz.Vec(m3.atom[i].coord)
	v = (a1+a2+a3)/3
	if u.length() > v.length():
		u, v = v, u
	return v-u


def axisofrot(sel1, sel2, sel3, nsamp=100):
	m1 = cmd.get_model(sel1).atom[0]
	m2 = cmd.get_model(sel2).atom[0]
	m3 = cmd.get_model(sel3).atom[0]
	m1 = xyz.Vec(m1.coord)
	m2 = xyz.Vec(m2.coord)
	m3 = xyz.Vec(m3.coord)
	return ((m2-m1).cross(m2-m3)).normalized()
	# a, wtot = xyz.Vec(0), 0
	# for i in range(int(nsamp)):
	#    u = axisofrot1(m1, m2, m3)
	#    a  += u
	#    wtot += u.length()
	# a /= wtot
	return a

#   0.373875069479 0.733798169028 0.567236881341 24.7756576007

def bond_zns(sel):
	cmd.unbond(sel + " and resn ZNS", sel + " and not resn ZNS")
	# for c, i in getres("resn ZNS"):
	#    cmd.bond(sel+" and chain %s and resi %i and name ZN1"%(c, i),
	#             sel+" and chain %s and resi %i and name S*"%(c, i))
	for c, i in getres("name ZN1"):
		cmd.bond(sel + " and chain %s and resi %i and name ZN1" % (c, i),
					sel + " and chain %s and resi %i and name S*, N2, N4" % (c, i))
	allsg = cmd.get_model(sel + " and (resn CYS and (name CB))").atom
	allsf = cmd.get_model(sel + " and (resn ZHC and (name S*, C1, C4))").atom
	while allsg:
		closest = [9e9, None, None]
		for sga in allsg:
			for sfa in allsf:
				sg = xyz.Vec(sga.coord)
				sf = xyz.Vec(sfa.coord)
				if (sg - sf).length() < closest[0]:
					closest = [(sg-sf).length(), sga, sfa]
		sga, sfa = closest[1:]
		allsg.remove(sga)
		allsf.remove(sfa)
		if closest[0] > 10.0: break
		cmd.bond(sel+" and resi %s and chain %s and name %s"%(sga.resi, sga.chain, sga.name),

					sel+" and resi %s and chain %s and name %s"%(sfa.resi, sfa.chain, sfa.name) )


def trans(sel,v):
	if   xyz.isvec(v): cmd.translate([v.x,v.y,v.z], sel, 0, 0)
	elif xyz.isnum(v): cmd.translate([ v , v , v ], sel, 0, 0)
	else: raise NotImplementedError
def transx(sel,x):
	if xyz.isnum(x): cmd.translate([x,0,0], sel, 0, 0)
	else: raise NotImplementedError
def transy(sel,y):
	if xyz.isnum(y): cmd.translate([0,y,0], sel, 0, 0)
	else: raise NotImplementedError
def transz(sel,z):
	if xyz.isnum(z): cmd.translate([0,0,z], sel, 0, 0)
	else: raise NotImplementedError

def rot(sel, axis, ang, cen=xyz.V0):
	if not xyz.isvec(axis): raise NotImplementedError
	if not xyz.isnum(ang): raise NotImplementedError
	if not xyz.isvec(cen): raise NotImplementedError	
	# if abs(axis.x) < 0.00001: axis.x = 0.0
	# if abs(axis.y) < 0.00001: axis.y = 0.0
	# if abs(axis.z) < 0.00001: axis.z = 0.0
	cmd.rotate([round(axis.x,6), round(axis.y,6), round(axis.z,6)], ang, sel, 0, 0, None, [cen.x, cen.y, cen.z])
def rotx(sel, ang, cen=xyz.V0):
	if not xyz.isnum(ang): raise NotImplementedError
	if not xyz.isvec(cen): raise NotImplementedError	
	cmd.rotate([1,0,0], ang, sel, 0, 0, None, [cen.x, cen.y, cen.z])
def roty(sel, ang, cen=xyz.V0):
	if not xyz.isnum(ang): raise NotImplementedError
	if not xyz.isvec(cen): raise NotImplementedError	
	cmd.rotate([0,1,0], ang, sel, 0, 0, None, [cen.x, cen.y, cen.z])
def rotz(sel, ang, cen=xyz.V0):
	if not xyz.isnum(ang): raise NotImplementedError
	if not xyz.isvec(cen): raise NotImplementedError	
	cmd.rotate([0,0,1], ang, sel, 0, 0, None, [cen.x, cen.y, cen.z])

def xform(sel,x):
	if not xyz.isxform(x): raise NotImplementedError
	axis,ang = x.rotation_axis()
	rot(sel,axis,xyz.degrees(ang))
	trans(sel,x.t)
	# a,r,c = x.rotation_center()
	# cmd.rotate([a.x,a.y,a.z], r, sel, 0, 0, None, [c.x, c.y, c.z])

def rot_by_matrix(sel, R, cen = xyz.Vec(0, 0, 0)):
  m = cmd.get_model(sel)
  for i in range(len(m.atom)):
		tmp = ( R * (xyz.Vec(m.atom[i].coord)-cen) ) + cen
		m.atom[i].coord = [tmp.x, tmp.y, tmp.z]
  cmd.load_model(m, sel, 1, discrete = 1)

def rotrad(sel, axis, ang, cen = None):
	return rot(sel, axis, ang*180.0/math.pi, cen)

def pointaxis(sel):
	u = xyz.Vec(1, 1, 1)
	m = cmd.get_model(sel)
	c = com(sel)
	axis = xyz.Vec(0)
	for a in m.atom:
		v = xyz.Vec(a.coord)-c
		if v.dot(u) < 0: v = -v
		axis += v
	return axis.normalized()


def alignaxis(sel, newaxis, oldaxis = None, cen = xyz.Vec(0, 0, 0)):
	if oldaxis is None: oldaxis = pointaxis(sel)
	if cen is None: cen = com(sel)
	newaxis.normalize()
	oldaxis.normalize()
	if abs(oldaxis.dot(newaxis)) > 0.99999: return
	axis = newaxis.cross(oldaxis).normalized()
	ang  = -math.acos(max(-1.0, min(1.0, newaxis.dot(oldaxis))))*180/math.pi
	# print "rot around", axis, "by", ang
	if str(ang) == "nan": return
	rot(sel, axis, ang, cen)
	return axis, ang

# def showvec(x, y, z, l = 100, c = None):
#    if c is None: c = com("all")
#    v = xyz.Vec(x, y, z).normalized()
#    a1 = c + v*l/2
#    a2 = c - v*l/2
#    obj = [
#          BEGIN, LINES,
#          COLOR, 1.0, 1.0, 1.0,
#          VERTEX, a1.x, a1.y, a1.z,
#          VERTEX, a2.x, a2.y, a2.z,
#          END
#    ]
#    cmd.load_cgo(obj, 'vec')
#


def mysetview(look, up, pos = None, cen = None, ncp = None, fcp = None):
	Xaxis = -look.cross(up).normalized()
	Yaxis = xyz.projperp(look, up).normalized()
	Zaxis = -look.normalized()
	oldv = cmd.get_view()
	v = list(oldv)
	r = xyz.Mat(Xaxis, Yaxis, Zaxis)
	v[0] = r.xx; v[1] = r.xy; v[2] = r.xz
	v[3] = r.yx; v[4] = r.yy; v[5] = r.yz
	v[6] = r.zx; v[7] = r.zy; v[8] = r.zz
	if pos is not None: v[ 9:12] = pos.x, pos.y, pos.z
	if cen is not None: v[12:15] = cen.x, cen.y, cen.z
	if ncp is not None: v[15]    = ncp
	if fcp is not None: v[16]    = fcp
	cmd.set_view(v)
	return Yaxis

def meancoords(sel1, sel2, n = "mix", w = 0.5):
	cmd.delete(n)
	a = [xyz.Vec(x.coord) for x in cmd.get_model(sel1).atom]
	b = [xyz.Vec(x.coord) for x in cmd.get_model(sel2).atom]
	d = [(a[i]-b[i]).length() for i in range(len(a))]
	#print max(d), argmax(d), cmd.get_model("179L").atom[argmax(d)].resi
	#print min(d)
	cmd.create(n, sel1)
	m = cmd.get_model(sel1)
	for j in range(len(m.atom)):
		m.atom[j].coord[0] = w*b[j].x + (1.0-w)*a[j].x
		m.atom[j].coord[1] = w*b[j].y + (1.0-w)*a[j].y
		m.atom[j].coord[2] = w*b[j].z + (1.0-w)*a[j].z
	cmd.load_model(m, n, 1)


def mygetview():
	v = cmd.get_view()
	return xyz.Vec(v[2], v[5], v[8])

def swell():
	a = [xyz.Vec(x.coord) for x in cmd.get_model("177L").atom]
	b = [xyz.Vec(x.coord) for x in cmd.get_model("179L").atom]
	d = [(a[i]-b[i]).length() for i in range(len(a))]
	#print max(d), argmax(d), cmd.get_model("179L").atom[argmax(d)].resi
	#print min(d)
	for i in range(0, 101):
		cmd.create("mix%03i"%i, "177L")
		m = cmd.get_model("177L")
		r = float(i)/100.0
		print "mixing", r
		for j in range(len(m.atom)):
			# for k in range(len(m.atom)):
			#   d = (a[j]-a[k]).length() - (b[j]-b[k]).length()
			# print j, m.atom[j].coord
			m.atom[j].coord[0] = r*b[j].x + (1.0-r)*a[j].x
			m.atom[j].coord[1] = r*b[j].y + (1.0-r)*a[j].y
			m.atom[j].coord[2] = r*b[j].z + (1.0-r)*a[j].z
			# print j, m.atom[j].coord
		cmd.load_model(m, "mix%03i"%i, 1)
		cmd.save("177L_179L_raw_mix%03i.pdb"%i, "mix%03i"%i)


def mkhelix(sel, t, r, n):
	v = cmd.get_view()
	for i in range(1, n):
		cmd.delete("h%i"%i)
		cmd.create("h%i"%i, sel)
		rot  ("h%i"%i, xyz.Vec(0, 0, 1), i*r, xyz.Vec(0, 0, 0))
		trans("h%i"%i, xyz.Vec(0, 0, i*t))
	cmd.set_view(v)

def mkhelix4(sel, t, r, n):
	for i in range(n):
		cmd.delete("h%i"%i)
		cmd.create("h%i"%i, sel)
		rt = 90.0 * (i%4)
		print i, rt
		if i%4==0 and i > 0:
			cmd.create("hem%i"%i, "basehem")
			trans("hem%i"%i, xyz.Vec(0, 0, t*(int(i/4))))
			rot  ("hem%i"%i, xyz.Vec(0, 0, 1), r*(int(i/4)), xyz.Vec(0, 0, 0))
			if i in (4, 12, 20, 28, 36, 44, 52):
				rot  ("hem%i"%i, xyz.Vec(0, 0, 1), 90, xyz.Vec(0, 0, 0))
		rot  ("h%i"%i, xyz.Vec(0, 0, 1), rt, xyz.Vec(0, 0, 0))
		trans("h%i"%i, xyz.Vec(0, 0, t*(int(i/4))))
		rot  ("h%i"%i, xyz.Vec(0, 0, 1), r*(int(i/4)), xyz.Vec(0, 0, 0))
	cmd.hide("ev", "not (name N, CA, C, O, CB) and (h1, h3, h4, h6, h9, h11, h12, h14, h17, h19, h20, h22, h25, h27)")
	cmd.hide("ev", "name SG and not (h0, h5, h8, h13, h16, h21, h24, h29, h32, h37, h40, h45, h48, base)")
#PyMOL>run /Users/sheffler/pymol_scripts/util.py; mkhelix4("4hel* and vis", 10.2, -23.9, 12)


def mirror(sel, nname = "mirror", crd = 0):
	cmd.delete(nname)
	a = [xyz.Vec(x.coord) for x in cmd.get_model(sel).atom]
	#print max(d), argmax(d), cmd.get_model("179L").atom[argmax(d)].resi
	#print min(d)
	cmd.create(nname, sel)
	m = cmd.get_model(sel)
	for j in range(len(m.atom)):
		m.atom[j].coord[crd] *= -1
	cmd.load_model(m, nname, 1)

def inversion(sel, nname = "inv"):
	cmd.delete(nname)
	a = [xyz.Vec(x.coord) for x in cmd.get_model(sel).atom]
	#print max(d), argmax(d), cmd.get_model("179L").atom[argmax(d)].resi
	#print min(d)
	cmd.create(nname, sel)
	m = cmd.get_model(sel)
	for j in range(len(m.atom)):
		m.atom[j].coord[0] *= -1
		m.atom[j].coord[1] *= -1
		m.atom[j].coord[2] *= -1
	cmd.load_model(m, nname, 1)


def mkc4(sel, a = xyz.Vec(1, 0, 0), c = xyz.Vec(0, 0, 0)):
	cmd.delete("c1")
	cmd.delete("c2")
	cmd.delete("c3")
	cmd.delete("c4")
	cmd.create("c1", sel)
	cmd.create("c2", sel)
	cmd.create("c3", sel)
	cmd.create("c4", sel)
	rot("c2", a, 90, c)
	rot("c3", a, 180, c)
	rot("c4", a, 270, c)

def mkc3(sel, a = xyz.Vec(0, 0, 1), c = xyz.Vec(0, 0, 0)):
	cmd.delete("c1")
	cmd.delete("c2")
	cmd.delete("c3")
	cmd.create("c1", sel)
	cmd.create("c2", sel)
	cmd.create("c3", sel)
	rot("c2", a, 120, c)
	rot("c3", a, 240, c)
	cmd.alter('c1', 'chain = "A"')
	cmd.alter('c2', 'chain = "B"')
	cmd.alter('c3', 'chain = "C"')

def mkc2(sel, a = xyz.Vec(0, 0, 1), c = xyz.Vec(0, 0, 0)):
	cmd.delete("c1")
	cmd.delete("c2")
	cmd.create("c1", sel)
	cmd.create("c2", sel)
	rot("c2", a, 180, c)
	cmd.alter('c1', 'chain = "A"')
	cmd.alter('c2', 'chain = "B"')


def mkd2(sel = "all"):
	cmd.create("w", sel)
	cmd.create("x", "w")
	cmd.create("y", "w")
	cmd.create("z", "w")
	rot('x', Ux, 180, xyz.Vec(0, 0, 0))
	rot('y', Uy, 180, xyz.Vec(0, 0, 0))
	rot('z', Uz, 180, xyz.Vec(0, 0, 0))
	cmd.alter('w', 'chain = "A"')
	cmd.alter('x', 'chain = "B"')
	cmd.alter('y', 'chain = "C"')
	cmd.alter('z', 'chain = "D"')



def alignall(sel = "all", obj = None):
	l = cmd.get_object_list()
	if not obj: obj = l[0]
	if obj not in l:
		print "ERROR object", obj, "not found!!!"
		return
	for o in l:
		if o == obj: continue
		cmd.do("align "+o+" and ("+sel+"), "+obj+" and ("+sel+")")
	return

def centerall(sel = "all"):
	l = cmd.get_object_list()
	for o in l:
		s = "%s and (%s)"%(o, sel)
		trans(s, -com(s))
	return


def showcst(fname):
	for l in open(fname).readlines():
		kind, a1, r1, a2, r2 = l.split()[:5]
		a1 = "(resi %s and name %s)"%(r1, a1)
		a2 = "(resi %s and name %s)"%(r2, a2)
		cmd.bond(a1, a2)
		cmd.show("lines", a1)
		cmd.show("lines", a2)


def bondzn():
	for o in cmd.get_object_list():
		print "add zn bonds for", o
		for r, c in getres(o+" and elem ZN"):
			cmd.bond("(%s and resi %s and chain %s)"%(o, r, c), "(%s and resn HIS and elem N) within 2.5 of (%s and resi %s and chain %s)"%(o, o, r, c) )
		break

def loadmov(d):
	files = glob.glob(d+"/*.pdb")
	for f in files: cmd.load(f, "mov")

def drawlines(p, d, lab = "lines", COL = (1, 1, 1), SIZE = 20.0):
	if type(p) is type(xyz.Vec(1)): p = [p, ]
	if type(d) is type(xyz.Vec(1)): d = [d, ]
	assert len(p) == len(d)
	obj = [ BEGIN, LINES, COLOR, COL[0], COL[1], COL[2], ]
	for i in range(len(p)):
		p1 = p[i] + -SIZE*d[i]
		p2 = p[i] +  SIZE*d[i]
		obj.extend([
				VERTEX, p1.x, p1.y, p1.z,
				VERTEX, p2.x, p2.y, p2.z,
			])
	obj.append(END)
	cmd.load_cgo(obj, lab)


################### disulf cone tests ##################################

def drawtestconeva(v, a):
	cmd.delete("cone hyperb")
	# v = com(sele+" and name SG")
	# a = (v - com(sele+" and name CB")).normalized()
	# print v, a
	ang = 5.0
	R = xyz.rotation_matrix_degrees(a, ang)
	perp = a.cross(randvec()).normalized()
	p = v + perp
	d = xyz.rotation_matrix_degrees(perp, 45)*a
	P, D = [], []
	for i in range(int(360/ang)):
		P.append(p)
		D.append(d)
		p = R*(p-v) + v
		d = R*d
	drawlines(P, D, "hyperb", COL = (0, 0, 1))
	p = v
	d = xyz.rotation_matrix_degrees( a.cross(randvec()), 45 ) *a
	P, D = [], []
	for i in range(360/ang):
		P.append(p)
		D.append(d)
		d = R*d
	drawlines(P, D, "cone", COL = (1, 0, 0))


def conelineinter(p, d, v, a, t):
	t = t / 180.0 * math.pi
	xyz.Mat = a.outer(a) - math.cos(t)*math.cos(t)*xyz.Mat(1, 0, 0, 0, 1, 0, 0, 0, 1)
	D = p-v
	c2 =   d.dot(xyz.Mat*d)
	c1 = 2*d.dot(xyz.Mat*D)
	c0 =   D.dot(xyz.Mat*D)
	disc = c1*c1 - 4*c0*c2
	if disc < 0: return ()
	if disc == 0: return ( p + (-c1)/(2.0*c2)*d, )
	disc = sqrt(disc)
	return ( p+(-c1+disc)/(2.0*c2)*d , p+(-c1-disc)/(2.0*c2)*d )


def createpoint(sele, p, lab):
	cmd.create(lab, sele+" and name CA")
	trans(lab, p-com(lab))

def test_conelineinter(sele):
	v = com(sele+" and name SG")
	a = (v - com(sele+" and name CB")).normalized()
	print v, a
	p = v+5*randvec()
	d = randvec().normalized()
	# v = xyz.Vec(0, 0, 0)
	# a = xyz.Vec(1, 0, 0)
	t = 45
	Ux = conelineinter(p, d, v, a, t)
	print "p", p
	print "d", d
	print "v", v
	print "a", a
	print "t", t
	print "Ux", Ux
	cmd.delete("lines")
	cmd.delete("Ux*")
	cmd.delete("L*")
	cmd.delete("A*")
	cmd.delete("B*")
	drawlines(p, d)
	for i, x in enumerate(Ux):
		createpoint(sele, x, "Ux"+str(i))
		o = xyz.projperp(a, x-v)
		L = o.length()
		o = o.normalized()
		ang = 90.0 - math.atan( 1.0 / L )*180.0/math.pi
		print ang
		o1 = xyz.rotation_matrix_degrees(a,  ang )*o
		o2 = xyz.rotation_matrix_degrees(a, -ang )*o
		createpoint(sele, v+o1, "A"+str(i))
		createpoint(sele, v+o2, "B"+str(i))
		drawlines(x, (x-(v+o1)).normalized(), "xyzMath"+str(i))
		drawlines(x, (x-(v+o2)).normalized(), "LB"+str(i))


######################### rosettaholes stuff ###################

def useRosettaRadii():
	cmd.alter("element C", "vdw = 2.00")
	cmd.alter("element N", "vdw = 1.75")
	cmd.alter("element O", "vdw = 1.55")
	cmd.alter("element H", "vdw = 1.00")
	cmd.alter("element P", "vdw = 2.15")
	cmd.alter("element S", "vdw = 1.90")
	cmd.alter("element RE", "vdw = 1.40")
	cmd.alter("element CU", "vdw = 1.40")
	cmd.set("sphere_scale", 1.0)

def expandRadii(delta = 1.0, sel = 'all'):
	for a in cmd.get_model(sel).atom:
		r = float(a.vdw) + float(delta)
		cmd.alter("index "+`a.index`, 'vdw = '+`r`)
	cmd.rebuild(sel, "spheres")

def contractRadii(delta = 1.0, sel = 'all'):
	for a in cmd.get_model(sel).atom:
		r = float(a.vdw) - float(delta)
		cmd.alter("index "+`a.index`, 'vdw = '+`r`)
	cmd.rebuild(sel, "spheres")

def useOccColors(sel = "all"):
	d = {}
	for a in cmd.get_model().atom:
		d[a.q] = True
	colors = rainbow
	# random.shuffle(colors)
	# colors *= len(d)/len(colors)+1
	for ii in range(len(d.keys())):
		cmd.color( colors[ii] , "%s and q = %i"%(sel, d.keys()[ii]))

def useTempColors(sel = "all"):
	for a in cmd.get_model(sel).atom:
		q = a.b
		c = intcolors[ int(floor(q))%len(intcolors) ]
		cmd.color( c , "%s and resi %s and name %s"%(sel, a.resi, a.name))


def useOccRadii(sel = "all"):
	for a in cmd.get_model(sel).atom:
		q = a.q
		if q >= 3:
			print "shrik radius"
			q <- 0.1
		cmd.alter("%s and resi %s and name %s"%(sel, a.resi, a.name), "vdw = %f"%(q))
	cmd.rebuild()

def useTempRadii(sel = "all"):
	for ii in range(30):
		radius = "%0.1f"%(float(ii+1)/10)
		cmd.alter(sel+" and b = "+radius, "vdw = "+radius)
	cmd.rebuild()



######## general utils #############

def natom(sel = "all"):
	return cmd.select(sel)

def nca(sel = "all"):
	return natom("( "+sel+" and (name ca))")

def chaincount(sel = "all"):
	cc = []
	for c in getchain(sel):
		cc.append((cmd.select("( "+sel+") and (chain %s) and (not (HET and not resn MSE+CSW))"%c), c))
	cc.sort()
	return cc

name1 = {"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F", "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", "LEU":"L",
		 "MET":"xyz.Mat", "ASN":"N", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S", "THR":"T", "VAL":"xyz.Vec", "TRP":"W", "TYR":"Uy",
		 "MSE":"xyz.Mat", "CSW":"C"}

def lcs(S, T):
	 L = {}
	 z = 0
	 ret = 0, -1, 0, -1
	 for i in range(len(S)):
		 for j in range(len(T)):
			 L[(i, j)] = 0
			 if S[i] == T[j]:
				 if i == 0 or j == 0:
					 L[(i, j)] = 1
				 else:
					 L[(i, j)] = L[(i-1, j-1)] + 1
				 if L[(i, j)] > z:
					 z = L[i, j]
					 ret = i-z+1, i, j-z+1, j
	 return ret

################## symmetry utils ###################


#############################################################################################################

# unscramble DX symmetry, seems to have a bug
# computes COMs, tries to find the double-ring pattern
# NO GOOD: think of rotational offset
def sort_dcoms(subs):
  N = len(subs)/2
  for i in range(1, N+1): assert "tmp%iA"%i in cmd.get_object_list()
  for i in range(1, N+1): assert "tmp%iB"%i in cmd.get_object_list()
  coms = [com(s) for s in subs]
  tmp = [(coms[0].distance(x), i) for i, x in enumerate(coms)]
  tmp.sort()
  td = tmp[2][0]-tmp[1][0] < tmp[3][0]-tmp[2][0]
  ring = [ tmp[1 if td else 2][1] , 0 , tmp[2 if td else 3][1] ]
  while len(ring) < N:
	 tmp = [(coms[ring[-1]].distance(x), i) for i, x in enumerate(coms)]
	 tmp.sort()
	 assert (tmp[2][0]-tmp[1][0] < tmp[3][0]-tmp[2][0]) == td
	 v1 = tmp[1 if td else 2][1]
	 v2 = tmp[2 if td else 3][1]
	 assert ring[-2] in (v1, v2)
	 ring.append( v1 if v2 == ring[-2] else v2 )
  #print ring
  #print [subs[i] for i in ring]
  namemap = {}
  for i, r in enumerate(ring):
	 assert not subs[r] in namemap
	 namemap[subs[r]] = "sub%iA"%(i+1)
	 tmp = [(coms[r].distance(x), j) for j, x in enumerate(coms)]
	 tmp.sort()
	 #print r, [subs[x[1]] for x in tmp]
	 sub_partner = subs[tmp[3 if td else 1][1]]
	 assert not sub_partner in namemap
	 namemap[sub_partner] = "sub%iB"%(i+1)
  assert len(set(namemap.keys())) == 2*N
  for i in range(1, N+1): assert "tmp%iA"%i in cmd.get_object_list()
  for i in range(1, N+1): assert "tmp%iB"%i in cmd.get_object_list()
  chains = {}
  for i in range(N):
	 chains["sub%iA"%(i+1)] = string.uppercase[2*N+0]
	 chains["sub%iB"%(i+1)] = string.uppercase[2*N+1]
  for k in namemap:
	 cmd.set_name(k, namemap[k])
	 cmd.alter(namemap[k], "chain = '%s'"%chains[namemap[k]])
  for i in chains: assert i in cmd.get_object_list()

def procD5dat(lfile = None, biod = "/data/biounit", outd = None):
  N = 5
  if outd is None: outd = "/Users/sheffler/project/sym_comp/D%i"%N
  print outd
  Nnobio = 0; Nok = 0; Ncontact = 0; Nbig = 0; Nnsym = 0; Nnomxatm = 0; Nhomogen = 0
  files = open(lfile).readlines() if lfile else map(os.path.basename, glob.glob(biod+"/*"))
  for fn in files:
	 fn = fn.strip()
	 pdb = os.path.basename(fn)[:4].lower()
	 bnum = int(fn[-1:])
	 print fn, pdb, bnum
	 if os.path.exists(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb") or os.path.exists(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb.gz"):
		Nok += 1
		continue
	 fname = biod+"/"+fn
	 if not os.path.exists(fname): fname += ".gz"
	 if not os.path.exists(fname):
		Nnobio += 1
		continue
	 cmd.delete("all")
	 cmd.load(fname, 'm')
	 cmd.remove("resn HOH")
	 cmd.remove('not alt a+""')
	 hf = cmd.select("(HET and not resn MSE+CSW)", state = 1) / cmd.select("ALL", state = 1)
	 if hf > 0.1: continue
	 cmd.remove("(HET and not resn MSE+CSW) and not name n+ca+c+o+cb")
	 assert 0 == cmd.select("all", state = 2*N)
	 assert 0 == cmd.select("all", state=  N)
	 if cmd.select('all', state = 2) != 0:
		cc = chaincount('m')
		if len(cc) < N:
		  Nnsym += 1
		  continue
		for i in range(1, N+1):
		  cmd.create("sub%iA"%i, "m and chain %s"%(cc[i-1][1]), 1, 1)
		  cmd.create("sub%iB"%i, "m and chain %s"%(cc[i-1][1]), 2, 1)
		continue
	 else:
		cc = chaincount("m")
		if len(cc) < 2*N:
		  Nnsym += 1
		  continue
		for i, l in enumerate([str(i)+c for i in range(1, N+1) for c in "AB"]):
		  cmd.create("tmp%s"%l, "m and chain %s"%(cc[i][1]), 1, 1)
		sort_dcoms(["tmp"+str(i)+c for i in range(1, N+1) for c in "AB"])
	 for slct in ["sub"+str(i)+c for i in range(1, N+1) for c in "AB"]:
		 if iscontig(slct):
			 cmd.create("mxatm", slct)
			 break
	 if cmd.select("name CA and mxatm") > 250:
		 Nbig += 1
		 continue
	 if cmd.select("mxatm") < 50:
		 Nnomxatm += 1
		 continue
	 chains = ["sub"+str(i)+c for i in range(1, N+1) for c in "AB"]
	 done = False
	 count = 0
	 # while not done and count < 50:
	 #    done = True
	 #    random.shuffle(chains)
	 #    for i in range(len(chains)):
	 #       for j in range(i+1, len(chains)):
	 #          done = done and homogenizechains(chains[i], chains[j])
	 #    count += 1
	 # if count >= 50:
	 #    Nhomogen += 1
	 #    continue
	 # if cmd.select("mxatm") < 50:
	 #    Nnomxatm += 1
	 #    continue
	 # cm = com("sub*")
	 # for s in chains: trans(s, -cm)
	 # sel = " or ".join(chains)
	 # a1 = c5axis("sub*", None, ["A", "C", "E", "G", "I"])
	 # alignaxis(sel, xyz.Vec(0, 0, 1), a1, xyz.Vec(0, 0, 0))
	 # a2 = c2axis("sub*", None, ["A", "B"])
	 # a2.z = 0.0
	 # a2.normalize()
	 # alignaxis(sel, xyz.Vec(0, 1, 0), a2, xyz.Vec(0, 0, 0))
	 # cmd.align("mxatm", "sub1A")
	 # if not os.path.exists(outd): os.mkdir(outd)
	 # cmd.save(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb", "mxatm")
	 # Nok += 1
  print Nok, Nbig, Nnsym, Ncontact, Nnobio, Nnomxatm, Nhomogen


def untangle_sidechains(sele):
	for c, i in getres(sele):
		if c != "": cmd.unbond("%s and resi %i and chain %s and not name N+C"%(sele, i, c), "not (%s and resi %i and chain %s)"%(sele, i, c))
		else:       cmd.unbond("%s and resi %i and              not name N+C"%(sele, i  ), "not (%s and resi %i             )"%(sele, i  ))

def orb_cyl(lab = ""):
	cmd.delete(lab+"o1")
	cmd.delete(lab+"o2")
	p1 = com('pk1')
	p2 = com('pk2')
	p3 = com('pk3')
	dr = (p1-p2).normalized()
	ax = (p3-p2).cross(p1-p2).normalized()
	o1 = xyz.rotation_matrix_degrees(ax, 60.0)*dr
	o2 = xyz.rotation_matrix_degrees(ax,-60.0)*dr
	cmd.load_cgo(o1.cgofrompoint(p1), lab+"o1")
	cmd.load_cgo(o2.cgofrompoint(p1), lab+"o2")


def drawring( p1 = None, p2 = None, p3 = None, col = [1, 1, 1], lab = "ring"):
	if p1 is None: p1 = com('pk1')
	if p2 is None: p2 = com('pk2')
	if p3 is None: p3 = com('pk3')
	cmd.delete(lab)
	axs = (p2-p1).normalized()

	obj = [ BEGIN, LINE_LOOP,   COLOR, col[0], col[1], col[2],  ]
	for i in range(0, 360, 5):
		st = xyz.rotation_matrix_degrees(axs, i  )*(p3-p2)+p2
		obj.extend( [ VERTEX, st.x, st.y, st.z, ] )

	obj.append( END )

	cmd.load_cgo(obj, lab)


def drawringcar( c, a, r, col = [1, 1, 1], lab = "ring"):
	cmd.delete(lab)
	p1 = c
	p2 = c+a
	p3 = c + r*xyz.projperp(a, xyz.Vec(1, 2, 3)).normalized()
	drawring(p1, p2, p3, col, lab)

def drawsph(col = [1, 1, 1], lab = "sph"):
	cmd.delete(lab)
	p1 = com('pk1')
	p2 = com('pk2')
	p3 = com('pk3')
	p4 = com('pk4')
	axs1 = (p2-p1).normalized()
	axs2 = (p3-p2).normalized()

	obj = [ BEGIN, ]
	for i in range(0, 360, 10):
		obj.extend( [ LINE_LOOP,   COLOR, col[0], col[1], col[2],  ] )
		axs = xyz.rotation_matrix_degrees(axs1, i)*axs2
		for j in range(0, 360, 3):
			st = xyz.rotation_matrix_degrees(axs, j)*(p4-p2)+p2
			obj.extend( [ VERTEX, st.x, st.y, st.z, ] )
	obj.append( END )

	cmd.load_cgo(obj, lab)

def dsf(CA1, CB1, CA2, CB2, lab = ''):
	c = CA1+(CB1-CA1).normalized()*2.27887
	a = (CB1-CA1).normalized()
	r = 1.6454076
	drawringcar(c, a, r, [1, 1, 0], 'cr'+lab)
	d = a.dot(c-CB2)
	r2 = sqrt( 3.0*3.0 - d*d )
	c2 = CB2 + d*a
	a2 = a
	drawringcar(c2, a2, r2, [1, 1, 1], 'cd'+lab)
	d = (c-c2).length()
	d2 = (r*r+d*d-r2*r2)/2/d
	x = d2 * (c2-c).normalized() + c
	h = sqrt(r*r - d2*d2)
	a3 = a.cross(c2-c).normalized()
	(x+h*a3).show('p1'+lab)
	(x-h*a3).show('p2'+lab)


def alignallrms(sele):
	r = {}
	for i in cmd.get_object_list():
		r[i] = cmd.align('fr52re', i)[0]
	for k, v in r.items():
		if v < 2: print k, v

def mkpntx(s1, s2):
	x = com(s1)-com(s2)
	if abs(x.x) > 0.1 or abs(x.y) > 0.1:
		print "DIE!", x
		return
	z = x.z
	c = xyz.Vec(0, z/2/math.tan(36.0*math.pi/180), z/2)
	cmd.delete('p1');   cmd.delete('p2');   cmd.delete('p3');   cmd.delete('p4');   cmd.delete('p5'   )
	cmd.create('p1', s1);cmd.create('p2', s1);cmd.create('p3', s1);cmd.create('p4', s1);cmd.create('p5', s1)
	rot('p1', xyz.Vec(1, 0, 0), 0*72, c)
	rot('p2', xyz.Vec(1, 0, 0), 1*72, c)
	rot('p3', xyz.Vec(1, 0, 0), 2*72, c)
	rot('p4', xyz.Vec(1, 0, 0), 3*72, c)
	rot('p5', xyz.Vec(1, 0, 0), 4*72, c)
	cmd.color('green'  , 'p1 and elem C')
	cmd.color('cyan'   , 'p2 and elem C')
	cmd.color('yellow' , 'p3 and elem C')
	cmd.color('magenta', 'p4 and elem C')
	cmd.color('orange' , 'p5 and elem C')


def ifsphab():
	cmd.show("spheres", "((vis and chain A) within 6 of (vis and chain B)) or ((vis and chain B) within 6 of (vis and chain A))")

def ifsphab():
	cmd.show("sticks", "byres (((vis and chain C) within 8 of (vis and chain A)) or ((vis and chain B) within 8 of (vis and chain A)))")
	cmd.show("sticks", "chain A")


def charge(sel):
	p = cmd.select(sel+" and name CA and resn LYS+ARG")
	n = cmd.select(sel+" and name CA and resn GLU+ASP")
	print p-n

def redoA(sel = "not sub", N = None):
	cmd.delete("sub*")
	cmd.hide("ev", sel+" and not chain A")
	v = cmd.get_view()
	chains = list(set(a.chain for a in cmd.get_model("name ca").atom))
	chains.sort()
	if N: chains = chains[:N]
	for c in chains:
		if c == "A": continue
		cmd.create("sub"+c, sel+" and chain A")
		cmd.align( "sub"+c, sel+" and chain "+c)
		cmd.alter( "sub"+c, "chain = '%s'"%c)
	color_by_chain()
	cmd.set_view(v)
	print charge(sel+" and chain A")

def redotoi(sel = "not sub*"):
	v = cmd.get_view()
	cmd.delete("sub*")
	cmd.create("subB", sel)
	cmd.create("subC", sel)
	axs = com('chain A and name ZN')
	rot("subB", axs, 120)
	rot("subC", axs, 240)
	cmd.alter("subB and chain A", "chain = 'E'")
	cmd.alter("subB and chain B", "chain = 'F'")
	cmd.alter("subB and chain C", "chain = 'G'")
	cmd.alter("subB and chain D", "chain = 'H'")
	cmd.alter("subC and chain A", "chain = 'I'")
	cmd.alter("subC and chain B", "chain = 'J'")
	cmd.alter("subC and chain C", "chain = 'K'")
	cmd.alter("subC and chain D", "chain = 'L'")
	color_by_chain()
	cmd.set_view(v)
	print charge("chain A")

def redopent(sel):
	v = cmd.get_view()
	cmd.remove("chain B+C+D+E")
	cmd.delete("pnt*")
	cmd.create("pntB", sel)
	cmd.create("pntC", sel)
	cmd.create("pntD", sel)
	cmd.create("pntE", sel)
	rot("pntB", xyz.Vec(1, 0, 0), 72)
	rot("pntC", xyz.Vec(1, 0, 0), 144)
	rot("pntD", xyz.Vec(1, 0, 0), 216)
	rot("pntE", xyz.Vec(1, 0, 0), 288)
	cmd.color("green"   , sel+" and elem C")
	cmd.color("cyan"   , "pntB and elem C")
	cmd.color("magenta", "pntC and elem C")
	cmd.color("yellow" , "pntD and elem C")
	cmd.color("pink"   , "pntE and elem C")
	cmd.alter("pntB", "chain = 'B'")
	cmd.alter("pntC", "chain = 'C'")
	cmd.alter("pntD", "chain = 'D'")
	cmd.alter("pntE", "chain = 'E'")
	cmd.set_view(v)
	print charge("chain A")

def getaa(c, r, o = 'all'):
	m = cmd.get_model("%s and chain %s and resi %d"%(o, c, int(r)))
	return aa_3_1[m.atom[0].resn]

def mkifaceresfile(fn = None):
	sel = "all"
	cmd.delete("all")
	cmd.load(fn)
	fn += ".resfile"
	o = sys.stdout
	if fn: o = open(fn, 'w')
	print >>o, """AUTO\nNATRO\n\nstart"""
	for c, r in getres("(chain A and (%s)) within 5 of (chain B and (%s))"%(sel, sel)):
		assert c == "A"
		print >>o, "%4d A PIKAA %s"%(r, getaa(c, r))
	if fn: o.close()


def color_by_chain():
	chains = set(a.chain for a in cmd.get_model("name ca").atom)
	print chains
	for i, c in enumerate(chains):
		cmd.color(COLORS[i], "chain %s and elem C"%c)



def getnative():
	v = cmd.get_view()
	nats = []
	for obj in cmd.get_object_list():
		if len(obj) == 4:
			nats.append(obj)
			continue
		pid = obj[:4]
		if pid in nats: continue
		print "fetching native", pid
		cmd.fetch(pid)
		cmd.remove(pid+" and not chain A")
		cmd.remove(pid+" and resn HOH")
		cmd.align(pid, obj)
	cmd.set_view(v)



def makecx(sel = 'all', n = 5):
	cmd.delete("C%i_*"%n)
	chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
	for i in range(n): cmd.create("C%i_%i"%(n, i), sel+" and (not C%i_*)"%n)
	for i in range(n): rot("C%i_%i"%(n, i), Uz, 360.0*float(i)/float(n))
	for i in range(n): cmd.alter("C%i_%i"%(n, i), "chain = '%s'"%chains[i])

for i in range(2,21):
		globals()['makec%i'%i] = partial(makecx,n=i)

def makecxauto():
	for o in cmd.get_object_list():
		n = int(re.search("_C\d+_", o).group(0)[2:-1])
		makecx(o, n)


def floats2vecs(i):
	while 1: yield xyz.Vec(i.next(), i.next(), i.next())

def nnb(v, s, r):
	nnb = 0
	for x in s:
		if v.distance(x) < r: nnb += 1
	return nnb

def testhsphere(rratio = 2.0):
	def by3(i):
		while 1: yield i.next(), i.next(), i.next()
	cmd.delete("hsphere*")
	lvl = list()
	crd = list()
	nbr = list()
	cld = list()
	pnt = list()
	if not hasattr(rratio, '__iter__'):
		vmult = [( 1.0 + (rratio-1.0)*(1.0-(float(l)-1.0)/6.0) ) for l in range(8)]
	elif len(rratio) is 7:
		vmult = [0]+rratio
	with gzip.GzipFile("/Users/sheffler/Dropbox/project/sphere_hierarchy/hsphere.dat.gz") as f:
		for l, m, n in by3(f):
			l, p, x, y, z = l.split()
			lvl.append(int(l))
			pnt.append(int(p))
			crd.append( xyz.Vec(float(x), float(y), float(z)) * vmult[int(l)] )
			cld.append([int(x) for x in m.split() if x != "65535"])
			nbr.append([int(x) for x in n.split() if x != "65535"])
			# print lvl
			# print pnt
			# print crd
			# print cld
			# print nbr
			# return

	obj = list()
	obj.extend( [ BEGIN, POINTS , COLOR, 1.0, 0.0, 0.0 ] )
	for v in crd:
		if v.z >= 0:
			obj.extend([VERTEX, 1.001*v.x, 1.001*v.y, 1.001*v.z ])
	obj.append(END)
	obj.extend( [ BEGIN, LINES, COLOR, 0.2, 0.2, 0.2, ]  )
	for v, c, l in itertools.izip(crd, cld, lvl):
			if l > 5: continue
			#if v.z >= 0:
			for vc in (crd[x-1] for x in c):
					#if vc.z > 0:
					obj.extend([VERTEX,  v.x,  v.y , v.z ])
					obj.extend([VERTEX, vc.x, vc.y, vc.z ])
	obj.append(END)
	cmd.load_cgo(obj, "hsphere")

def testsphere():
	cmd.delete("s*")
	cmd.delete('icos')
	spheres = []
	icos = []
#   for i, n in enumerate((72, 132, 252, 492, 972, 1932, 3812, 7002, 13232, 25232, 50432)):
	#for i, n in enumerate((72, 92, 122, 132, 162, 192, 212, 252, 272, 282, 312, 362, 372, 392, 432, 482, 492, 522, 572, 612, 632, 642, 672, 732, 752, 762, 792, 812, 842, 912, 932, 972, )):
	for n in (32, 122, 482, 1922, 7682, 30722):
	#for n in (32, 92, 272, 812, 2432):#, 482, 1922, 7682, 30722):
		f = gzip.GzipFile("/Users/sheffler/Dropbox/project/sphere_hierarchy/att/sphere_%i.dat.gz"%n)
		s = list(floats2vecs(float(x) for x in f))
		f.close()
		obj = [ BEGIN, POINTS , COLOR, 1.0, 1.0, 1.0 ]
		doneicos = True if icos else False
		for v in s:
			#if (v.x < 0 or v.y < 0 or v.z < 0): continue
			if not doneicos and nnb(v, s, 0.5) == 6: icos.append(v)
			#if v.x >= 0 and v.y >= 0 and v.z >= 0: obj.extend([VERTEX,  v.x,  v.y,  v.z ])
			if v.z >= 0: obj.extend([VERTEX,  v.x,  v.y,  v.z ])
			#obj.extend([ VERTEX, v.x*0.86, v.y*0.86, v.z*0.86 ])
		obj.append(END)
		cmd.load_cgo(obj, "s%i"%n)
		spheres.append(s)
		#break

	obj = [ BEGIN, LINES, COLOR, 1.0, 1.0, 1.0, ]
	for v in icos:
		for u in icos:
			if v.distance(u) < 1.1:
				obj.extend([VERTEX, v.x, v.y, v.z])
				obj.extend([VERTEX, u.x, u.y, u.z])
				# if v.x >= 0 and v.y >= 0 and v.z >= 0: obj.extend([VERTEX, v.x, v.y, v.z])
				# if u.x >= 0 and u.y >= 0 and u.z >= 0: obj.extend([VERTEX, u.x, u.y, u.z])
	obj.append(END)
	cmd.load_cgo(obj, "icos")

	# samp = spheres[0]
	# d = min(m.distance(n) for m in samp for n in samp if m is not n) * 1.5
	# snew = list()
	# for v in samp:
	#   nb = [n for n in samp if 0.001 < v.distance(n) < d]
	#   for c in (((v+m+n)/3.0).normalized() for m in nb for n in nb if 0.001 < m.distance(n) < d):
	#       if not snew or min(x.distance(c) for x in snew) > 0.001:
	#           snew.append(c)
	# snew.extend(samp)
	# obj = [ BEGIN, POINTS, COLOR, 0.0, 0.0, 1.0, ]
	# print len(snew)
	# o = open("/Users/sheffler/Dropbox/project/sphere_hierarchy/sphere_%i.dat"%len(snew), 'w')
	# for v in snew:
	#   print >>o, v.x
	#   print >>o, v.y
	#   print >>o, v.z
	#   if v.z > 0: obj.extend([VERTEX, v.x, v.y, v.z])
	# o.close()
	# obj.append(END)
	# cmd.load_cgo(obj, "snew")

	# for s in sph:
	#   print len(s), angle_degrees(s[0],s[1])
	# test = list()
	# for i in range(len(sph[-2])):
	#   mn = 9e9; mni = -1
	#   for j in range(len(sph[-1])):
	#       d = sph[-2][i].distance(sph[-1][j])
	#       if d < mn:
	#           mn = d
	#           mni = j
	#   test.append( sph[-1][mni] )
	#   #sph[-1]    = sph[-1][:mni]+sph[-1][mni+1:]
	# obj = [ BEGIN, POINTS, COLOR, 1.0, 0.0, 0.0, ]
	# for v in test:
	#   obj.extend([ VERTEX, v.x, v.y, v.z ])
	# obj.append(END)
	# cmd.load_cgo(obj, "sphere1")
	# obj = [ BEGIN, POINTS, COLOR, 0.0, 1.0, 0.0, ]
	# for v in sph[-1]:
	#   obj.extend([ VERTEX, v.x, v.y, v.z ])
	# obj.append(END)
	# cmd.load_cgo(obj, "sphere2")

def stubalign(s = 'all', s1 = 'pk1', s2 = 'pk2', s3 = 'pk3'):
	trans(s, -com(s1))
	if -0.99999 < (com(s1)-com(s2)).dot(xyz.Vec(1, 0, 0)) < 0.99999: alignaxis(s, xyz.Vec(1, 0, 0), com(s1)-com(s2))
	else: print "alignaxis FUBAR!"
	y = com(s3).y
	z = com(s3).z
	if   abs(z) < 0.00001: ax =   0.0 if y > 0 else 180.0
	elif abs(y) < 0.00001: ax = 270.0 if z > 0 else  90.0
	else:   ax = -math.atan(z/y)*180.0/math.pi
	rot(s, Ux, ax)



def tmpdoit():
	axs = {
		"TET":(xyz.Vec(0, 0, 1)+xyz.Vec( 0.942809043336065, 0                , -0.333333328372267)).normalized(),
		"OCT":(xyz.Vec(0, 0, 1)+xyz.Vec( 0.666666666666667, 0.666666666666667, 0.333333333333333)).normalized(),
		"ICS":(xyz.Vec(0, 0, 1)+xyz.Vec(-0.333333320519719, 0.57735021589134,  0.745356039515023)).normalized(),
	}
	for fn in glob.glob("/Users/sheffler/project/111108_BPY_TOI/picksdone/*.pse"):
		print fn
		cmd.delete("all")
		cmd.load(fn)
		sym = ""
		for o in cmd.get_object_list():
			if not (o.startswith("TET") or o.startswith("OCT") or o.startswith("ICS")):
				cmd.delete(o)
			else:
				sym = o[:3]
		if sym == "TET":
			cmd.remove("not vis")
			v = xyz.Vec(0.816496579408716, 0, 1.57735027133783)/2.0
			rot("all", v, 180.0)
			cmd.save(fn.replace(".pse", "_A.pdb").replace("picksdone", "for_neil"))

def symmetrizec2(sel1, sel2):
	v = cmd.get_view()
	cmd.delete("sA")
	cmd.delete("sB")
	print cmd.select(sel1)
	print cmd.select(sel2)
	cmd.create("sA", sel1)
	cmd.create("sB", sel2)
	cmd.alter("sA", "chain = 'A'")
	cmd.alter("sB", "chain = 'B'")
	a = c2axis("(sA or sB)")
	alignaxis("(sA or sB)", xyz.Vec(0, 0, 1), a, xyz.Vec(0, 0, 0))
	trans("sA", xyz.Vec(0, 0, -com("sA").z))
	trans("sB", xyz.Vec(0, 0, -com("sB").z))
	cmd.delete("sB")
	cmd.create("sB", "sA")
	cmd.alter("sB", "chain = 'B'")
	rot("sB", xyz.Vec(0, 0, 1), 180.0)
	cmd.set_view(v)




def corresponding_atom_names(sel1,sel2,file):
	with open(file) as f:
		s = f.read()
	m1 = cmd.get_model(sel1).atom
	m2 = cmd.get_model(sel2).atom
	fr=[]
	to=[]
	for a1 in m1:
		mn = 9e9
		mi = None
		for a2 in m2:
			d = xyz.Vec(a1.coord).distance(xyz.Vec(a2.coord))
			if d < mn:
				mi = a2
				mn = d
		fr.append(mi.name)
		to.append(a1.name.replace(" ",""))
	for i in range(len(fr)):
		s = re.sub("\\b"+fr[i]+"\\b",to[i],s)
	with open(file+".rename","w") as o:
		o.write(s)
	print "wrote",file+".rename"




def bestalign(s1,s2):
	ol1 = [x for x in cmd.get_object_list() if x.startswith(s1)]
	ol2 = [x for x in cmd.get_object_list() if x.startswith(s2)]
	keepers = []
	for o1 in ol1:
		mn = 9e9
		mo = None
		for o2 in ol2:
			r = cmd.align(o1,o2)[0]
			if r < mn:
				mn = r
				mo = o2
		cmd.align(o1,mo)
		keepers.append(mo)
		print "%7.3f    %10s %10s"%(mn,o1,mo)
	for o2 in ol2:
		if o2 not in keepers:
			cmd.delete(o2)



cmd.extend("com", com)
cmd.extend("showcom", showcom)
cmd.extend("axisofrot", axisofrot)
cmd.extend("bond_zns", bond_zns)
cmd.extend("meancoords", meancoords)
cmd.extend("swell", swell)
cmd.extend("axes",showaxes)
cmd.extend('useRosettaRadii', useRosettaRadii)
cmd.extend('expandRadii',expandRadii)
cmd.extend('contractRadii',contractRadii)
