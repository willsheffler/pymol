# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-

import sys,os,inspect
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path:
		sys.path.append(newpath)
import string, math, re, gzip, itertools, glob, sets
from random import randrange
from math import sqrt
import xyzMath as xyz
from xyzMath import Ux,Uy,Uz,Imat
from functools import partial
import cProfile



class Inspector(object):
	"""docstring for Inspector"""
	def __init__(self, l="/Users/sheffler/tmp/ariel10_p6_32_reverted/list",tgt="/Users/sheffler/tmp/ariel10_p6_32_reverted/pymol_picks/"):
		super(Inspector, self).__init__()
		self.pdblist = list()
		with open(l) as infile:
			for pdb in infile.readlines():
				self.pdblist.append(pdb.strip())
		self.index = -1
		self.tgt = tgt+"/"
		self.current = None
		if not os.path.exists(self.tgt):
			os.mkdir(self.tgt)
		if not os.path.exists(self.tgt+"seenit.list"):
			self.seenit = []
		else:
			with open(self.tgt+"seenit.list") as infile:
				self.seenit = [x.strip() for x in infile.readlines()]

	def write_seenit(self):
		with open(self.tgt+"seenit.list","w") as outfile:
			for x in self.seenit: outfile.write(x+"\n")

	def next(self):
		if self.current:
			cmd.delete(self.current)
			cmd.delete(self.nat1)
			cmd.delete(self.nat2)
			cmd.delete("all")
			self.seenit.append(self.pdblist[self.index])
			self.write_seenit()

		while self.pdblist[self.index] in self.seenit:
			self.index += 1
			if self.index >= len(self.pdblist):
				print "out of structures!"
				return

		cmd.load(self.pdblist[self.index])
		self.current = cmd.get_object_list()[-1]

		cmd.select("comp1","chain A+C+E")
		cmd.select("comp2","chain B+H+L")
		cmd.select("ocompA1","chain I+G")
		cmd.select("ocompA2","chain J+F")
		cmd.select("ocompB1","chain K+M")
		cmd.select("ocompB2","chain D+N")
		cmd.select("iface","byres ((comp1 within 8 of comp2) or (comp2 within 8 of comp1) or (comp1 within 8 of ocompA1) or (ocompA1 within 8 of comp1) or (comp2 within 8 of ocompA2) or (ocompA2 within 8 of comp2))")
		print "inspecting",self.current
		util.cbc()
		util.cnc()
		cmd.show('car')
		cmd.orient()
		cmd.zoom('iface')
		cmd.show('sti','iface')
		cmd.hide('ev','hydro')
		cmd.select('None')
		v = cmd.get_view()

		self.nat1 = self.current[ 0: 9]
		self.nat2 = self.current[10:19]
		print self.nat1, self.nat2
		cmd.load("/work/sheffler/data/jacob_"+self.nat1[:2]+"/"+self.nat1+".pdb.gz")
		cmd.load("/work/sheffler/data/jacob_"+self.nat2[:2]+"/"+self.nat2+".pdb.gz")
		cmd.super(self.nat1,self.current+" and chain C")
		cmd.super(self.nat2,self.current+" and chain A")
		cmd.color("white","(%s or %s) and elem C"%(self.nat1,self.nat2))
		cmd.hide("ev","(%s or %s)"%(self.nat1,self.nat2))
		cmd.show("lines","(%s or %s)"%(self.nat1,self.nat2))
		cmd.set_view(v)

	def keep(self):
		if self.current:
			os.system("cp "+self.pdblist[self.index] + " " + self.tgt)

def yes():
	global inspect
	inspect.keep()
	inspect.next()
def no():
	global inspect
	inspect.next()

cmd.extend("yes",yes)
cmd.extend("no",no)