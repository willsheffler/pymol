# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import os,sys
newpath = os.path.abspath(os.path.dirname(__file__))
if not newpath in sys.path: sys.path.append(newpath)
import string,re,gzip,itertools
from pymol import cmd
from LA import *
from sym_util import *
from pymol_util import *

def iscontig(sel):
	m = cmd.get_model(sel+" and name N+CA+C").atom
	for i in range(1,len(m)):
		if ( Vec(m[i-1].coord) - Vec(m[i].coord) ).length() > 1.8:  return False
		return True

def procCdat(N=3,lfile=None,biod="/data/biounit",outd=None):
	if lfile is None: lfile=os.path.expanduser("~/Dropbox/project/sym_comp/meta/C%i.dat"%N)
	if outd  is None:  outd=os.path.expanduser(".//C%i"%N)
	print outd
	Nnobio=0; Nok=0;  Nbig=0; Nnsym=0; Nnomxatm=0; Nhomogen=0
	for fn in open(lfile).readlines():
		try:
				fn = fn.split()[0].strip()
				pdb = fn[3:7]
				bnum = int(fn[-1:])
				if os.path.exists(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb") or os.path.exists(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb.gz"): 
					Nok += 1
					continue
				fname = biod+"/"+fn
				if not os.path.exists(fname): fname += ".gz"
				if not os.path.exists(fname):
					Nnobio += 1
					print "no file",fname
					continue
				cmd.delete("all")
				cmd.load(fname,'m')
				cmd.remove("resn HOH")
				cmd.remove('not alt a+""')
				#hf = cmd.select("(HET and not resn MSE+CSW)",state=1) / cmd.select("ALL",state=1)
				#if hf > 0.1: continue
				#cmd.remove("(HET and not resn MSE+CSW)")
				if cmd.select('all',state=N) != 0:
					for i in range(1,N+1):
						cmd.create("sub%i"%i,"m and not (HET and not resn MSE+CSW)",i,1)
				else:
					cc = chaincount("m")
					if len(cc) < N: 
						Nnsym += 1
						print "ERROR:",fname," symmetry error"
						continue
					for i in range(1,N+1):
						cmd.create("sub%i"%i,"m and chain %s and not (HET and not resn MSE+CSW)"%(cc[-i][1]),1,1)
				for i in range(1,N+1):
					if iscontig("sub%i"%i):
						cmd.create("mxatm","sub%i"%i)
						break
				if cmd.select("mxatm") < 50: 
					Nnomxatm += 1
					print "ERROR: mxatm < 50"
				 	continue
				if cmd.select("name CA and mxatm") > 500: 
					Nbig += 1
					print "ERROR:",fname," more than 500 res"
					continue
				chains = ["sub%i"%i for i in range(1,N+1)]
				done = False
				count = 0
				while not done and count < 50:
					done = True
					random.shuffle(chains)
					for i in range(len(chains)):
						for j in range(i+1,len(chains)):
							done = done and homogenizechains(chains[i],chains[j])
					count += 1
				if count >= 50: 
					Nhomogen += 1
					print "ERROR: error homogenizing"
					continue
				if cmd.select("sub1") < 50: 
					Nnomxatm += 1
					print "ERROR: less than 50 atoms"
					continue
				cm = com("sub*")
				for i in range(1,N+1):
					trans("sub%i"%i,-cm)
				a = [cmd.get_model("sub%i and name CA"%i).atom for i in range(1,N+1)]
				axis = Vec(0,0,0)
				for i in range(len(a[0])):
					axis1 = Vec(0,0,0)
					for j in range(N): axis1 += Vec(a[j][i].coord)
					if axis1.length() > 0.0001 and axis.dot(axis1) < 0: axis1 *= -1
					axis += axis1		
				axis.normalize()
				for i in range(1,N+1):
					alignaxis("sub%i"%i,Vec(0,0,1),axis,Vec(0,0,0))
				#cmd.create("final1","mxatm")
				#cmd.create("final2","mxatm")
				#cmd.create("final3","mxatm")
				#cmd.align("final1","sub1")
				#cmd.align("final2","sub2")
				#cmd.align("final3","sub3")
				#return
				if not os.path.exists(outd): os.mkdir(outd)
				cmd.align("mxatm","sub1")
				cmd.save(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb","mxatm")
				Nok += 1
				print "SUCCESS on",N,fname
		except (pymol.CmdException,Exception) as e:
			print "EXCEPTION!!!!!!",e
			print fname
			pass
	print Nok, Nbig, Nnsym, Nnobio, Nnomxatm, Nhomogen

def procD2dat(lfile=None,biod="/data/biounit",outd=None):
	N = 4
	if lfile is None: lfile='/Users/sheffler/project/sym_comp/meta/D2.list'
	if outd  is None:  outd="/Users/sheffler/project/sym_comp/D2"
	print outd
	Nnobio=0; Nok=0; Ncontact=0; Nbig=0; Nnsym=0; Nnomxatm=0; Nhomogen=0
	for fn in open(lfile).readlines():
		fn = fn.strip()
		pdb = fn[3:7]
		bnum = int(fn[-1:])
		if os.path.exists(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb") or os.path.exists(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb.gz"): 
			Nok += 1
			continue
		fname = biod+"/"+fn
		if not os.path.exists(fname):
		  fname += ".gz"
		if not os.path.exists(fname):
		  Nnobio += 1
		  continue
		#print pdb,bnum,fname
		cmd.delete("all")
		cmd.load(fname,'m')      
		cmd.remove("resn HOH")
		cmd.remove('not alt a+""')
		#hf = cmd.select("(HET and not resn MSE+CSW)",state=1) / cmd.select("ALL",state=1)
		#if hf > 0.1: continue
		#cmd.remove("(HET and not resn MSE+CSW)")
		if   cmd.select('all',state=4) != 0:
			for i in range(1,N+1):
				cmd.create("sub%i"%i,"m and not (HET and not resn MSE+CSW)",i,1)
		elif cmd.select('all',state=2) != 0:
			cc = chaincount('m')
			if len(cc) < 2:
				Nnsym += 1
				continue
			cmd.create("sub1","m and chain %s and not (HET and not resn MSE+CSW)"%(cc[0][1]),1,1)
			cmd.create("sub2","m and chain %s and not (HET and not resn MSE+CSW)"%(cc[1][1]),1,1)
			cmd.create("sub3","m and chain %s and not (HET and not resn MSE+CSW)"%(cc[0][1]),2,1)
			cmd.create("sub4","m and chain %s and not (HET and not resn MSE+CSW)"%(cc[1][1]),2,1)
		else:
			cc = chaincount("m")
			if len(cc) < N:
				sym = cmd.get_symmetry("m")
				if   sym[6] == "I 2 2 2":
					trans('m',Vec(0,-sym[1],0))
					print pid
				#elif sym[6] == "P 21 21 2" and len(cc)==2:
				#   trans('m',Vec(-sym[0]/2.0,-sym[1]/2.0,0))
				#   cmd.create("sub1","m and chain %s and not (HET and not resn MSE+CSW)"%(cc[0][1]),1,1)
				#   cmd.create("sub2","m and chain %s and not (HET and not resn MSE+CSW)"%(cc[1][1]),1,1)
				#   cmd.create("sub3","m and chain %s and not (HET and not resn MSE+CSW)"%(cc[0][1]),1,1)
				#   cmd.create("sub4","m and chain %s and not (HET and not resn MSE+CSW)"%(cc[1][1]),1,1)
				#   rot("sub3",Vec(0,0,1),180,Vec(0,0,0))
				#   rot("sub4",Vec(0,0,1),180,Vec(0,0,0))
				elif sym[6] in ('C 1 2 1','P 21 21 21','P 62 2 2','P 64 2 2','P 65 2 2','P 63 2 2','P 61 2 2','C 2 2 21'):
					Nnsym += 1
					#if pid != "1y2k_2": return
					continue           
				else:
					Nnsym += 1
					continue#return
				cmd.save(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb","m")
				Nok += 1
				continue
			else:
				for i in range(1,N+1):
					cmd.create("sub%i"%i,"m and chain %s and not (HET and not resn MSE+CSW)"%(cc[-i][1]),1,1)
		for i in range(1,N+1):
			if iscontig("sub%i"%i):
				cmd.create("mxatm","sub%i"%i)
				break
		if cmd.select("name CA and mxatm") > 250: 
			Nbig += 1
			continue
		if cmd.select("mxatm") < 50: 
			Nnomxatm += 1
			continue
		chains = ["sub%i"%i for i in range(1,N+1)]
		done = False
		count = 0
		while not done and count < 50:
			done = True
			random.shuffle(chains)
			for i in range(len(chains)):
				for j in range(i+1,len(chains)):
					done = done and homogenizechains(chains[i],chains[j])
			count += 1
		if count >= 50: 
			Nhomogen += 1
			continue
		if cmd.select("sub1") < 50: 
			Nnomxatm += 1
			continue      
		cm = com("sub*")
		for i in range(1,N+1):
			trans("sub%i"%i,-cm)
		a = [cmd.get_model("sub%i and name CA"%i).atom for i in range(1,N+1)]
		a1 = Vec(0,0,0)
		for i in range(len(a[0])):
			axis1 = Vec(a[0][i].coord) + Vec(a[1][i].coord)
			if axis1.length() > 0.0001 and a1.dot(axis1) < 0: axis1 *= -1
			a1 += axis1
		a1.normalize()
		for i in range(1,N+1):
			alignaxis("sub%i"%i,Vec(1,0,0),a1,Vec(0,0,0))
		a = [cmd.get_model("sub%i and name CA"%i).atom for i in range(1,N+1)]
		a1 = Vec(0,0,0)
		for i in range(len(a[0])):
			axis1 = Vec(a[0][i].coord) + Vec(a[2][i].coord)
			if axis1.length() > 0.0001 and a1.dot(axis1) < 0: axis1 *= -1
			a1 += axis1
		a1.normalize()
		for i in range(1,N+1):
			alignaxis("sub%i"%i,Vec(0,1,0),a1,Vec(0,0,0))
		cmd.align("mxatm","sub1")
		cmd.create("final2","mxatm")
		cmd.create("final3","mxatm")
		cmd.create("final4","mxatm")
		rot('final2',Vec(1,0,0),180,Vec(0,0,0))
		rot('final3',Vec(0,1,0),180,Vec(0,0,0))
		rot('final4',Vec(0,0,1),180,Vec(0,0,0))
		n1 = cmd.select('mxatm within 4 of final2')
		n2 = cmd.select('mxatm within 4 of final3')
		n3 = cmd.select('mxatm within 4 of final4')
		if n1 < 10 and n2 < 10 and n3 < 10:
			Ncontact += 1
			continue
		if not os.path.exists(outd): os.mkdir(outd)
		cmd.save(outd+"/"+pdb+"_"+str(bnum)+"_sub1.pdb","mxatm")
		Nok += 1
		#return
	print Nok, Nbig, Nnsym, Ncontact, Nnobio, Nnomxatm, Nhomogen



if __name__ == '__main__': 
	# pymol environment
	moddir='/home/sheffler/pymol/modules'
	sys.path.insert(0, moddir)
	os.environ['PYMOL_PATH'] = "/home/sheffler/pymol"
 	import pymol
	pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
	pymol.finish_launching()
	cmd = pymol.cmd
	for i in range(2,9): procCdat(i)
