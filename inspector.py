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
	def __init__(self, l="/Users/sheffler/runs/I53/130930_matdes_fixbb/filt_sc_area_nobup_5k.list",tgt="/Users/sheffler/runs/I53/130930_matdes_fixbb/pymol_picks/"):
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
		if self.current: cmd.delete(self.current)

		while self.pdblist[self.index] in self.seenit:
			self.index += 1
			if self.index >= len(self.pdblist):
				print "out of structures!"
				return
		self.seenit.append(self.pdblist[self.index])
		self.write_seenit()

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