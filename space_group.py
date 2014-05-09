import os, inspect

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
	def op2xform(self,op):
		print op
		raise NotImplementedError

path = os.path.dirname(inspect.getfile(inspect.currentframe()))
spacegroups = list()
with open(path+'/data/symop.lib') as f:
	while f:
		header = f.readline().strip()
		comment = ''
		if header.find("!") >= 0:
			header,comment = header.split(" ! ")
		header = header.split()
		idx, nop, dunno, name1, name2, celltype = header[:6]
		idx, nop, dunno = int(idx), int(nop), int(dunno)
		names = " ".join(header[6:]).strip("'").split("' '")
		ops = list()
		for iop in range(nop):
			ops.append( f.readline().strip() )
		print idx,names,len(ops),comment
		spacegroups.append( SpaceGroup( names, celltype, ops ) )
		