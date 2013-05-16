import sys,os,gzip,cPickle,random
sys.path.append("/Users/sheffler/pymol")
from xyzMath import *
from pymol_util import *

BIOUNIT_DIR="/data/pdb/biounit/"

FORBID_DEFAUT = set([
"1o061",
"1o062",
"1o063",
'3bbz3',
'4dac5',
'2xkm1',
])

MAXCLEN=50
MAXTLEN=50

def get_pdbids():
	thelix  = set()
	psmall  = set()
	pnosym  = set()
	phetero = set()
	pbadman = set()
	with open("/Users/sheffler/dropbox/project/pdbcatalog/biounit_term_helix_ids.txt")    as F: [thelix.add(l[:4]) for l in F]
	with open("/Users/sheffler/dropbox/project/pdbcatalog/biounit_term_helix_nres.txt")   as F: [psmall.add(l[:4]) for l in F if (int(l.split()[1]) <= MAXCLEN and int(l.split()[2]) <= MAXTLEN)]
	with open("/Users/sheffler/dropbox/project/pdbcatalog/3dcomplex_NPS.txt")             as F: [pnosym.add(l[:4]) for l in F]
	with open("/Users/sheffler/dropbox/project/pdbcatalog/pisa_heteromers_ids.txt")       as F: [phetero.add(l[:4]) for l in F]
	with open("/Users/sheffler/dropbox/project/pdbcatalog/biounit_term_helix_badman.txt") as F: [pbadman.add(l[:4]) for l in F]
	pids = set(thelix)
	pids.intersection_update(psmall)
	# pids.intersection_update(phetero)
	# pids.intersection_update(pnosym)
	pids = list(pids.difference(pbadman))
	pids.sort()
	# print "Nthelix:",len(thelix), ";  Nsize:",len(psmall), ";  Nnps3dc:",len(pnosym), ";  Npbadman:",len(pbadman), ";  Nhits:",len(pids)
	return pids,phetero,pnosym



class TH(object):
	"""
	T_H_HIT_HEADER        ID   pdb_b pchn  pres  chn Nchn clen nres  stres cenres endres                   CEN_X                   CEN_Y                   CEN_Z                   AXS_X                   AXS_Y                   AXS_Z                   ORI_X                   ORI_Y                   ORI_Z
	TERM_HELIX_HIT 101m1_HC1      101m1 A   140   -1    1  154  154    134    141    149       27.43285494578981       5.165086809753125      -1.288777309947627    -0.04107042507923700     -0.9961428002295765     0.07754187084787641     -0.9000069395055572    0.003178218272984028     -0.4358639785190436
	"""
	def __init__(self, line):
		super(TH, self).__init__()
		dat = line.split()
		# for i,d in enumerate(dat):
		# 	print i,d
		# sys.exit()
		self.id    = dat[ 1]
		self.pdb   = dat[ 2]
		self.chain = dat[ 3]
		# assert len(self.chain)==1
		self.pres  = int(dat[ 4])
		self.chain_num  = int(dat[ 5])
		self.num_chains = int(dat[ 6])
		self.chain_len  = int(dat[ 7])
		self.nres   = int(dat[ 8])
		self.stres  = int(dat[ 9])
		self.cenres = int(dat[10])
		self.endres = int(dat[11])
		CEN_X = float(dat[12])
		CEN_Y = float(dat[13])
		CEN_Z = float(dat[14])
		AXS_X = float(dat[15])
		AXS_Y = float(dat[16])
		AXS_Z = float(dat[17])
		ORI_X = float(dat[18])
		ORI_Y = float(dat[19])
		ORI_Z = float(dat[20])
		self.cen = Vec(CEN_X,CEN_Y,CEN_Z)
		self.axs = Vec(AXS_X,AXS_Y,AXS_Z)
		self.ori = Vec(ORI_X,ORI_Y,ORI_Z)
	def pdb_bounds(self):
		return self.pres-self.cenres+self.stres,  self.pres-self.cenres+self.endres,
	def show(self,x=Xform()):
		tag = self.id[:5]
		if not tag in cmd.get_object_list():
			load_biounit("%s/%s/%s.pdb%s.gz"%(BIOUNIT_DIR,self.pdb[1:3],self.pdb[:4],self.pdb[4]),tag)
			cmd.remove("het and "+tag)
			cmd.hide('ev',tag)
			xform(tag,x)
		# pofst = self.pres - self.cenres
		# print self.stres-pofst, self.endres-pofst
		# cmd.show('car',tag+" and resi %i-%i"%(self.stres-pofst,self.endres-pofst))
		# cmd.show('sti',tag+" and resi %i-%i"%(self.stres-pofst,self.endres-pofst))
		showvecfrompoint(x.R*self.axs*20,x*self.cen-x.R*self.axs*10,lbl=self.id+"_HA")
		showvecfrompoint(x.R*self.ori*05,x*self.cen                ,lbl=self.id+"_HO")
	def __getitem__(self,k): return getattr(self,k)
	def __str__(self):
		return "TH( %(id)s %(chain)s%(pres)s %(stres)s-%(cenres)s-%(endres)s %(chain_num)s %(axs)s %(ori)s %(cen)s )"%self

class HH(object):
	"""
	#  HELIX_JUNCTION   101m1_HC1 101m1_HN2   1   2   141  16 CN    20  12   0 -10   0  10   170.768906     0.746851         0.298094    -0.620917     0.724984         0.000000     0.000000     0.000000        68.704990    30.614476    -3.059844
	#  HELIX_JUNCTION
	#  hits[ihit].pdb "_H" (hits[ihit].chain<0?"C":"N") ihit
	#  hits[jhit].pdb "_H" (hits[jhit].chain<0?"C":"N") jhit
	#  ihit
	#  jhit
	#    a.cenres+(int)Ashift
	#    a.stopres-a.startres+1
	#    (a.chain<0?"C":"N") (b.chain<0?"C":"N")   b.cenres+(int)Bshift   b.stopres-b.startres+1   Aoffset   Boffset   Ashift   Bshift   ang   talonga   raxis.x()   rotcen.x()   x.t.x()
	"""
	def __init__(self,line,hconnect):
		super(HH, self).__init__()
		dat = line.split()
		# for i,d in enumerate(dat):
		# 	print i,d
		# sys.exit()
		self.h1 = hconnect[dat[1]]
		self.h2 = hconnect[dat[2]]
		if abs(self.h1.cenres-self.h2.cenres) < 20: raise RuntimeError
		self.res1 = int(dat[3])
		self.res2 = int(dat[4])
		self.id = self.h1.id+dat[3]+self.h2.id+dat[4]
		rax = float(dat[ 8]); ray = float(dat[ 9]); raz = float(dat[10])
		xtx = float(dat[11]); xty = float(dat[12]); xtz = float(dat[13])
		self.x = Xform(rotation_matrix_degrees(Vec(rax,ray,raz),float(dat[7])),Vec(xtx,xty,xtz))
	def show(self,x=Xform()):
		self.h1.show(x)
		self.h2.show(x)
		axs,cen,nf = rotcen(x*self.x)
		showvecfrompoint(axs*50,cen-axs*25)
		return axs,cen,nf
	def __str__(self):
		h1 = "{0}_{1}{2}".format(self.h1.pdb,self.h1.chain,self.h1.pres)
		h2 = "{0}_{1}{2}".format(self.h2.pdb,self.h2.chain,self.h2.pres)
		return "HH( {0} {1} {2} {3} {4} )".format(h1,h2,self.res1,self.res2,self.x.pretty())

def cachefiles():
	return "/Users/sheffler/tmp/hterm_%03i-%03i.dat.gz"%(MAXCLEN,MAXTLEN),"/Users/sheffler/tmp/hjunct_%03i-%03i.dat.gz"%(MAXCLEN,MAXTLEN)
def datafiles():
	c1,c2 = cachefiles()
	if os.path.exists(c1) and os.path.exists(c2): return c1,c2,True
	return "/Users/sheffler/dropbox/project/pdbcatalog/hterm.dat.gz","/Users/sheffler/dropbox/project/pdbcatalog/hjunct.dat.gz",False

def load_data(Nh=999999999999,Nj=999999999999,reload=True):
	pids,het,sym = get_pdbids()
	# for i in range(min(15,len(pids)/30)): print " ".join(pids[30*i:30*i+30])
	print "npdbs:",len(pids)
	have_data_cache = datafiles()[2]
	if have_data_cache: print "USE CACHED DATA"
	hconnect = {}
	with   gzip.open(datafiles ()[0]    ) as F:
	  with gzip.open(cachefiles()[0],'a') as O:
		for line in F:
			if "nan" in line: continue
			th = TH(line)
			# print line
			# print th
			# return
			if th.nres > MAXTLEN or th.chain_len > MAXCLEN:
				continue
			if th.pdb[:4] in pids:
				hconnect[th.id] = th
				if not have_data_cache: O.write(line)
				if len(hconnect) >= Nh: break
			if th.chain_num > th.num_chains:
			 	print line
			 	print th
				sys.exit("BAD CHAIN NUM!!!")
	print "hhconn:",len(hconnect)
	joints = []
	with   gzip.open(datafiles ()[1]    ) as F:
	  with gzip.open(cachefiles()[1],'a') as O:
		for line in F:
			if "nan" in line:
				continue
			try:
				hh = HH(line,hconnect)
				if hh.h1.num_chains > 2: continue
				if hh.h1.pdb in FORBID_DEFAUT: continue
				if abs(hh.h1.chain_num)==abs(hh.h2.chain_num) and hh.h1.num_chains > 1: continue
				joints.append(hh)
				if len(joints) >= Nj: break
				# if random.random() < 10.01:
					# print listne
					# print hh
					# print
					# print hh.h1
					# print hh.h2
					# print
					# return
				if not have_data_cache: O.write(line)
			except (KeyError,RuntimeError):
				pass
			# print line			# print hh			# print " ",hh.h1			# print " ",hh.h2			# break
	print "joints:",len(joints)
	# print "dumping"
	# with gzip.open("/Users/sheffler/tmp/biounit_term_helix_%03i-%03i.pickle.gz"%(MAXCLEN,MAXTLEN),'w') as out:
	# 	cPickle.dump((hconnect,joints),out,protocol=cPickle.HIGHEST_PROTOCOL)

	joints = [j for j in joints if j.h1.chain_num*j.h2.chain_num<0]

	return hconnect,joints


# def rotcen(x):
# 	axs,ang	= x.rotation_axis()
# 	nf    = int(round(2*pi/ang))
# 	cen,sub = Vec(0,0,0),Vec(0,0,0)
# 	for i in range(1,nf):
# 		sub = x*sub
# 		cen += sub
# 	cen /= nf
# 	return axs,cen,nf

def load_biounit(fn,obj):
	cmd.load(fn,obj)
	cmd.remove(obj+" and het")
	cmd.remove(obj+" and not name n+ca+c+o+cb")
	nstates = cmd.count_states(obj)
	chains = cmd.get_chains(obj)
	cindex = [ROSETTA_CHAINS.find(c.upper()) for c in chains]
	nchain = max(cindex)+1
	if nstates > 1:
		for s in range(1,nstates+1):
			for c in cindex:
				# print "load_biounit",s,c,ROSETTA_CHAINS[c]
				cmd.create(obj+"_STATE_%i_%i"%(s,c),obj,s,1)
				cmd.alter (obj+"_STATE_%i_%i"%(s,c),"chain='%s'"%ROSETTA_CHAINS[(s-1)*nchain+c])
		cmd.delete(obj)
		cmd.create(obj,obj+"_STATE_*")
		cmd.delete(obj+"_STATE_*")


def makesym(axs,cen,nf,objs=None,tag="makesym"):
	cmd.delete(tag+'_*')
	if objs is None:
		objs = [o for o in cmd.get_object_list()]
	for o in objs: cmd.hide("ev",o)
	for i in range(nf):
		state = 1
		for o in objs:
			for s in range(1,cmd.count_states(o)+1):
				# print "create",o,i,s
				cmd.create(tag+"_%i_%s_%s"%(i,o,s),o,s,state)
				# state += 1
		tmpobjs = [o for o in cmd.get_object_list() if o.startswith(tag+"_%i_"%i)]
		for j,o in enumerate(tmpobjs): cmd.alter(o,"chain='%s'"%ROSETTA_CHAINS[i*len(tmpobjs)+j])
		cmd.create(tag+"_%i"%i,tag+"_%i_*"%i,0,0)
		for j,o in enumerate(tmpobjs): cmd.delete(o)
		rot(tag+"_%i"%i,axs,i*360/nf,cen)
		# if i == 1:
		# 	clashsele = "({0}_0 within 2 of {0}_1) or ({0}_1 within 2 of {0}_0)".format(tag)
		# 	print clashsele
		# 	symclash = cmd.select(clashsele)
		# 	if symclash > 0:
		# 		print "symmetric clashes",symclash
		# 		cmd.delete(tag+"_*")
		# 		return False
		# renumber(tag+"_%i"%i)
		# chains = cmd.get_chains(tmpsel)
		# for j,c in enumerate(chains):
		# 	cmd.alter("tmp%i and chain %s"%(i,c),"chain='%s'"%ROSETTA_CHAINS[len(chains)*i+j])
	# for i in range(nf): cmd.create("tmp%i"%i,sel+" and (not tmp*)",1,1)
	for i in range(nf):	cmd.show("lines",tag+"_%i and name n+ca+c"%i)
	util.cbc()
	cmd.zoom()
	print "save /Users/sheffler/project/struts/picks/C%i_"%nf+"_".join(objs)+".pdb,"+tag+"_*"
	return True


# def symcheck(a,b,x,tol=1.0):
# 	axs,ang	= x.rotation_axis()
# 	skew = abs(x.t.dot(axs))
# 	if skew > 5*tol: return False
# 	nf    = round(2*pi/ang)
# 	# # if nf <= 2 or 30 < nf: return False
# 	# nferr = abs(2*pi/ang-nf)
# 	# angerr = abs(ang-2*pi/nf)
# 	# # if x.t.length()*sin(angerr) > tol: return False

# 	# print
# 	# print "SYMCHECK",skew
# 	# print "SYMCHECK",nferr
# 	# print "SYMCHECK",angerr
# 	# print "SYMCHECK",x.t.length()*sin(angerr)
# 	# print

# 	axs0,cen0,nf0 = symmetrize_xform(a.h1.cen,x,nf)
# 	xsym = rotation_around_degrees(axs0,360.0/float(nf0),cen0)
# 	err = (x * a.h2.cen  -   xsym * a.h2.cen).length_squared() + (x * b.h2.cen  -   xsym * b.h2.cen).length_squared()
# 	if err < tol: return True

# 	return False



def helix_joint_has_clashes(a,b):
	cmd.remove('not name n+ca+c+o+cb')
	ha,hb = a.h2,b.h1
	if   ha.chain_num < 0 and hb.chain_num > 0:
		s1 = "%s and not (chain %s and resi %i-999999)"%(ha.pdb,ha.chain,ha.pdb_bounds()[0]-4) # FIXME chain A
		s2 = "%s and not (chain %s and resi 000000-%i)"%(hb.pdb,hb.chain,hb.pdb_bounds()[1]+4) # FIXME chain A
	elif ha.chain_num > 0 and hb.chain_num < 0:
		s1 = "%s and not (chain %s and resi 000000-%i)"%(ha.pdb,ha.chain,ha.pdb_bounds()[1]+4) # FIXME chain A
		s2 = "%s and not (chain %s and resi %i-999999)"%(hb.pdb,hb.chain,hb.pdb_bounds()[0]-4) # FIXME chain A
	else:
		print "must align helices N->C or C->N!"
		raise NotImplementedError("must align helices N->C or C->N!")
	# print "helix_joint_has_clashes1",s1
	# print "helix_joint_has_clashes2",s2
	cmd.select("notH_a",s1)
	cmd.select("notH_b",s2)
	nclash = cmd.select("(notH_a within 4 of notH_b) or (notH_b within 4 of notH_a)")
	# print "helix_joint_has_clashes",nclash
	return nclash != 0

def helix_joint_has_symmetry_clashes(a,b,xaln,xsym):
	cmd.remove('not name n+ca+c+o+cb')
	ha,hb = a.h1,b.h2
	if   ha.chain_num > 0 and hb.chain_num < 0:
		s1 = "%s and not (chain %s and resi 000000-%i)"%(ha.pdb,ha.chain,ha.pdb_bounds()[1]+4) # FIXME chain A
		s2 = "%s and not (chain %s and resi %i-999999)"%(hb.pdb,hb.chain,hb.pdb_bounds()[0]-4) # FIXME chain A
	elif ha.chain_num < 0 and hb.chain_num > 0:
		s1 = "%s and not (chain %s and resi %i-999999)"%(ha.pdb,ha.chain,ha.pdb_bounds()[0]-4) # FIXME chain A
		s2 = "%s and not (chain %s and resi 000000-%i)"%(hb.pdb,hb.chain,hb.pdb_bounds()[1]+4) # FIXME chain A
	else:
		print "must align helices N->C or C->N!"
		raise NotImplementedError("must align helices N->C or C->N!")
	# cmd.delete("symclashtmp")
	cmd.create("symclashtmp","("+s1+") or ("+s2+")")
	xform("symclashtmp",xsym)
	symclashsel = "((symclashtmp within 2 of (("+s1+") or ("+s2+"))) or ((("+s1+") or ("+s2+")) within 2 of symclashtmp))"
	nclash = cmd.select(symclashsel)
	symclashsel = "((symclashtmp within 10 of (("+s1+") or ("+s2+"))) or ((("+s1+") or ("+s2+")) within 10 of symclashtmp))"
	ncontact = cmd.select(symclashsel)
	cmd.delete("symclashtmp")
	return nclash,ncontact



def gen3way(joints,forbid=set(FORBID_DEFAUT),nflb=2,nfub=99999,starti=0):
	forbid = set(forbid)
	cleanup = False
	for i in range(starti,len(joints)):
		a = joints[i]
		if a.h1.num_chains > 2: continue
		if abs(a.h1.chain_num)==abs(a.h2.chain_num) and a.h1.num_chains > 1: continue
		print "PERCENT_DONE", float(i)/float(len(joints))*100.0
		print "NEW i ==",i
		if a.h1.pdb in forbid: continue
		if a.h2.pdb in forbid: continue
		for j in range(len(joints)):
			if cleanup:
				cmd.delete('all')
				cleanup = False
			b = joints[j]
			if b.h1.num_chains > 2: continue
			if a.h1.chain_num*b.h1.chain_num < 0: continue
			if a.h1.pdb == b.h1.pdb: continue
			if b.h1.pdb in forbid: continue
			if b.h2.pdb in forbid: continue
			if abs(a.h1.chain_num)==abs(a.h2.chain_num) and abs(b.h1.chain_num)==abs(b.h2.chain_num): continue
			if abs(b.h1.chain_num)==abs(b.h2.chain_num) and b.h1.num_chains > 1: continue

			af = Xform().from_four_points(a.h1.cen,a.h1.axs,a.h1.ori,V0)
			bf = Xform().from_four_points(b.h1.cen,b.h1.axs,b.h1.ori,V0)
			xaln = a.x * af * ~bf
			x = xaln*b.x*~xaln*a.x

			axs,cen,nf = symmetrize_xform(a.h1.cen,x)
			if nflb > nf or nfub < nf: continue
			xsym = rotation_around_degrees(axs,360.0/float(nf),cen)
			ofsterr = math.sqrt(( (x * a.h2.cen  -   xsym * a.h2.cen).length_squared() + (x * b.h2.cen  -   xsym * b.h2.cen).length_squared() ) / 2.0)
			hsymdis = (xaln*b.h2.cen).distance(xsym*a.h1.cen)
			ha1 = xsym.R*a.h1.axs
			hc1 = xsym  *a.h1.cen
			ha2 = xaln.R*b.h2.axs
			hc2 = xaln  *b.h2.cen
			skew = line_line_distance(ha1,hc1,ha2,hc2)
			hang = line_line_angle_degrees(ha1,ha2)
			if skew + skew/hsymdis*10 > 1.0: continue
			if hang > 10.0: continue
			if ofsterr + ofsterr/hsymdis*10 > 1.5: continue

			# if not symcheck(a,b,x): continue

			cleanup = True
			# if a.h2.cen.distance(xaln*b.h1.cen) > 0.1: continue
			try:
				a.h1.show()
				a.h2.show()
				b.h1.show(xaln)
				b.h2.show(xaln)
			except IOError:
				continue
			# print "ah1axs,ah1ori =", repr(a.h1.axs), ',', repr(a.h1.ori)
			# print "ah2axs,ah2ori =", repr(a.h2.axs), ',', repr(a.h2.ori)
			# print "bh1axs,bh1ori =", repr(b.h1.axs), ',', repr(b.h1.ori)
			# print "bh2axs,bh2ori =", repr(b.h2.axs), ',', repr(b.h2.ori)
			# print "ax            =", repr(a.x)
			# print "bx            =", repr(b.x)
			# print "cx            =", repr(xaln)
			assert (xaln.R*b.h2.axs-x.R*a.h1.axs).length() < 0.0001

			if helix_joint_has_clashes(a,b):
				# print "CLASH",a.id+"_"+b.id
				continue
			else:
				print "NOCLASH",a.id+"_"+b.id


			# xform('all',xaln)

			nclash,ncontact = helix_joint_has_symmetry_clashes(a,b,xaln,xsym)
			if nclash > 0: continue

			# x,cen = symmetrize_xform(a.h1.cen,x,nf)
			# axs,ang = x.rotation_axis()
			showvecfrompoint(50*axs,cen)
			# print a
			# print b

			objs = [a.h1.pdb,b.h1.pdb]
			objs.sort()
			print "OBJS",objs

			makesym(axs,cen,nf,objs)

			fn = "/data/project/helix_assembly/cyclic/C%i_"%nf+"_".join(objs)+"_%i_%i.pdb"%(i,j)
			for inf in range(nf):
				trans    ("makesym_%i"%inf,-cen)
				alignaxis("makesym_%i"%inf,Uz,axs)
			cmd.save(fn,"makesym_*")
			os.system("gzip "+fn)

			nres_approx = cmd.select("makesym_0 and name ca")
			print "PERCENT_DONE", float(i)/float(len(joints))*100.0, i,j

			print "HIT",
			print fn,
			print ncontact,
			print nres_approx,
			# print skew,
			print hang,
			print ofsterr,
			print hsymdis

			message = ( yield (a,b,axs,cen,nf) )
			if message is not None:
				forbid.add(message)
				if a.h1.pdb in forbid: break
				if a.h2.pdb in forbid: break
				if b.h1.pdb in forbid: break
				if b.h2.pdb in forbid: break

	print "YOU FORBADE THESE:"
	for f in forbid:
		print f


if __name__ == '__main__':
	bh2axs,bh2ori = Vec(  0.157027,  0.088023,  0.983664 ) , Vec( -0.254474,  0.965993, -0.045819 )
	bh1axs,bh1ori = Vec( -0.058735, -0.183069,  0.981344 ) , Vec( -0.998029,  0.032544, -0.053662 )
	ah2axs,ah2ori = Vec(  0.036009,  0.336586, -0.940964 ) , Vec( -0.842618, -0.496016, -0.209672 )
	ah1axs,ah1ori = Vec(  0.287930,  0.310220, -0.906013 ) , Vec(  0.771555, -0.635564,  0.027581 )
	ax            = Xform( Mat( Vec( -0.902341,  0.379216, -0.204879 ), Vec( -0.184325, -0.769177, -0.611875 ), Vec( -0.389621, -0.514356,  0.763959 ) ), Vec(  90.645917, 24.612097, 40.190192 ) )
	bx            = Xform( Mat( Vec(  0.953091,  0.197334, -0.229515 ), Vec( -0.250173,  0.940406, -0.230325 ), Vec(  0.170387,  0.276939,  0.945660 ) ), Vec(   1.351979, -1.708957, 21.053139 ) )
	xaln          = Xform( Mat( Vec(  0.568684, -0.797590, -0.201118 ), Vec( -0.818463, -0.573041, -0.041746 ), Vec( -0.081953,  0.188348, -0.978677 ) ), Vec( -10.364693, 97.473637, 11.653443 ) )
	print ax.R*ah1axs - ah2axs
	print ah2axs - xaln.R*bh1axs
	print xaln.R*bh2axs  - (xaln.R*bx.R*~xaln.R*ax.R*ah1axs)



	import doctest
	r = doctest.testmod()
	print r
	# hconn,joints = load_data()
	# random.shuffle(joints)
	# for j in joints[:5]:
	# 	print j
	# 	print " ",j.h1
	# 	print " ",j.h2
	# 	print


def runsearch():
	global MAXCLEN,MAXTLEN
	MAXCLEN,MAXTLEN = 100,150
	print "RUNNING HELIX SEARCH"
	# tmp = sys.stdout
	# sys.stdout = open("/data/project/helix_assembly/cyclic/pymol_gen.log",'w')
	hconn,joints = load_data();
	gen = gen3way(joints,starti=304);
	for a,b,axs,cen,nf in gen:
		print "find another hit"
		cmd.delete('all')
		try:
			gen.next()
		except Exception as e:
			print "UNHANDLED EXCEPTION!",e
	# sys.stdout.close()
	# sys.stdout = tmp

runsearch()

"""
run /Users/sheffler/pymol/pick_helix_junctions.py; runsearch()

run /Users/sheffler/pymol/pick_helix_junctions.py; hconn,joints = load_data();

run /Users/sheffler/pymol/pick_helix_junctions.py; gen = gen3way(joints);

run /Users/sheffler/pymol/pick_helix_junctions.py; delete all; a,b,axs,cen,nf=gen.next(); hide rib; show line; util.cbc(quiet=True);zo; mysetview(axs,Uz)

run /Users/sheffler/pymol/pick_helix_junctions.py; makesym(axs,cen,nf)

"""
