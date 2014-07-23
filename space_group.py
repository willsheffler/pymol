import os, inspect, re, pickle, itertools
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
from xyzMath import Vec,Mat,Xform,Ux,Uy,Uz
from pymol_util import pymol,cgo_sphere,cube

class SpaceGroup(object):
	"""docstring for SpaceGroup"""
	def __init__(self, names, celltype, ops):
		super(SpaceGroup, self).__init__()
		self.name = names[0]
		self.names = names
		self.celltype = celltype
		self.ops = ops
		self.frames = list()
		for op in self.ops:
			self.frames.append( self.op2xform(op) )
	def op2vec(self,op):
		op = op.replace(" ","")
		v,t = Vec(), 0
		for o in op.replace("-","+-").strip("+").split("+"):
			if   o ==  'X': v += Ux
			elif o == '-X': v -= Ux
			elif o ==  'Y': v += Uy
			elif o == '-Y': v -= Uy
			elif o ==  'Z': v += Uz
			elif o == '-Z': v -= Uz
			elif o.count("/"):
				n,d = o.split("/")
				t += float(n) / float(d)
			else:
				raise ValueError
		return v,t
	def op2xform(self,op):
		Xop,Yop,Zop = op.split(",")
		# print op
		X,Tx = self.op2vec(Xop)
		Y,Ty = self.op2vec(Yop)
		Z,Tz = self.op2vec(Zop)
		# print X
		# print Y
		# print Z
		# print Tx,Ty,Tz
		# print
		return Xform( Mat( X.x, X.y, X.z,
		                   Y.x, Y.y, Y.z,
		                   Z.x, Z.y, Z.z
		             ), Vec(Tx,Ty,Tz) )
		# return Xform( Mat( X.x, Y.x, Z.x,
		# 	               X.y, Y.y, Z.y,
		# 	               X.z, Y.z, Z.z,
		#                    ), Vec(Tx,Ty,Tz) )
		# raise NotImplementedError
	def __str__(self):
		return "SpaceGroup %s, unit cell size: %i" % ( self.name, len(self.frames) )
	def show(self,points=[Vec(7,12,13),Vec(33,34,35)]):
		v = pymol.cmd.get_view()
		CGO = list()
		for x in self.frames:
			x.t = 100.0*x.t
			assert self.celltype == "CUBIC"
			for p in points:
				for i,j,k in itertools.product(range(-1,2),range(-1,2),range(-1,2)):
					CGO.extend( cgo_sphere( x*p + Vec(i,j,k)*100.0 ) )
		pymol.cmd.load_cgo( CGO, self.name )
		cube(Vec(0,0,0),Vec(100,100,100))
		pymol.cmd.set_view(v)

path = os.path.dirname(inspect.getfile(inspect.currentframe()))
picklefile = "data/spacegroups.pkl"
if False:#os.path.exists(picklefile):
	with open(picklefile) as infile:
		spacegroups = pickle.load(infile)
else:
	spacegroups = dict()
	with open(path+'/data/symop.lib') as f:
		while f:
			header = f.readline().strip()
			if not header: break
			header = re.sub("!+","!",header)
			# print header
			comment = ''
			if header.find("!") >= 0:
				splt = header.split("!")
				header,comment = splt[0],"".join(splt[1:])
				header = header.strip()
				comment = comment.strip()
			header = header.split()
			idx, nop, dunno, name1, name2, celltype = header[:6]
			idx, nop, dunno = int(idx), int(nop), int(dunno)
			names = " ".join(header[6:]).strip("'").split("' '")
			ops = list()
			for iop in range(nop):
				ops.append( f.readline().strip() )
			# print idx,names,len(ops),comment
			sg = SpaceGroup( names, celltype, ops )
			spacegroups[sg.name] = sg
	with open(picklefile,"w") as out:
		pickle.dump(spacegroups,out)

print "loaded",len(spacegroups),"spacegroups"
