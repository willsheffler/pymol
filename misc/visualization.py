# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import sys,os,inspect,time
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path: sys.path.append(newpath)
import string,re,gzip,itertools
from sym_util import *

# xform_movie_old("fere.pdb_*",3,"/Users/sheffler/project/symgen_movie/fere.pdb_xforms.dat")
def xform_movie_old_generator(sele,xforms,nfold):
	chains = ROSETTA_CHAINS
	undo = [Xform()]*nfold
	for i,x in enumerate(xforms):
		for n in range(nfold):
			sel = sele+" and chain "+chains[n]
			R = xyz.rotation_matrix_degrees(Vec(0,0,1),360.0*n/nfold)
			xform( sel, undo[n] )
			xrot = R * x * R.transposed()
			xform( sel, xrot )
			undo[n] = ~xrot
		a,r = x.rotation_axis()
		yield i,a,r

def xform_movie_old(sele,nfold,fn):
	cmd.remove("not name n+ca+c+o+cb")
	cmd.alter("all","vdw=1.75")
	cmd.rebuild()
	makecx(sele,nfold)
	cmd.delete(sele)
	cmd.alter("C3_0","chain='A'")
	cmd.alter("C3_1","chain='B'")
	cmd.alter("C3_2","chain='C'")
	cmd.create("trimer","C3_*")
	cmd.delete("C3_*")
	cmd.color('green'  ,'chain A')
	cmd.color('cyan'   ,'chain B')
	cmd.color('magenta','chain C')
	cmd.hide('ev'); cmd.show("spheres")
	xforms = xyz.read_xforms(fn)
	miter = xform_movie_old_generator("trimer",xforms,nfold)
	for i in range(len(xforms)):
	 	miter.next()
	 	cmd.create("mov","trimer",1,i+1)
	cmd.delete("trimer")

def nulltest():
	"""
	>>> print "foo"
	foo
	"""

def twocomp_movie_generator(fname="/Users/sheffler/Dropbox/test/tcdock/test.docks",symxforms=xyz.SYMTET,N=500):
	"""
	>>> twocompmovie()
	"""
	arch,f1,f2 = os.popen("head -n1 %s"%fname).read().split()
	print "twocompmovie",fname,arch,f1,f2
	x = xyz.read_xforms(fname,start=1)
	cmd.load(f1,"TEMPORARY_1")
	cmd.load(f2,"TEMPORARY_2")
	cmd.remove("not name n+ca+c+o+cb")
	cmd.alter("all","vdw=2.0")
	cmd.rebuild()
	for i in range(0,min(N*2,len(x)),2):
		# cmd.delete("not (_TEMPORARY_1 or _TEMPORARY_2)")
		cmd.create("TOXFORM_1","TEMPORARY_1")
		cmd.create("TOXFORM_2","TEMPORARY_2")
		xform("TOXFORM_1",x[i+0])
		xform("TOXFORM_2",x[i+1])
		for j,sx in enumerate(symxforms[:6]):
			cmd.create("sub_A_%i"%j,"TOXFORM_1")
			cmd.create("sub_B_%i"%j,"TOXFORM_2")
			cmd.alter("sub_A_%i"%j,"chain='%s'"%ROSETTA_CHAINS[2*j+0])
			cmd.alter("sub_B_%i"%j,"chain='%s'"%ROSETTA_CHAINS[2*j+1])
			xform("sub_A_%i"%j,sx)
			xform("sub_B_%i"%j,sx)
		cmd.delete("TOXFORM_1")
		cmd.delete("TOXFORM_2")
		cmd.create("arch%i"%(i/2),"sub_*")
		cmd.delete("sub_*")
		yield "arch%i"%(i/2)

def make_movie_from_generator(mgen):
	for i,obj in enumerate(mgen):
	 	cmd.create("mov",obj,1,i+1)
	 	cmd.delete(obj)
	cmd.delete("TEMPORARY_1")
	cmd.delete("TEMPORARY_2")



def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite())
    return tests

if __name__ == '__main__':
   import doctest
   r = doctest.testmod()
   print r
