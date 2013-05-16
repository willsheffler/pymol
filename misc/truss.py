# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import sys,os,inspect
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path: sys.path.append(newpath)
import string,re,gzip,itertools,operator
from pymol_util import *
from sym_util import *
from xyzMath import *
import random

"""
the only thing you need is probably the math, such as it is, here:
	rotation_around_dof_to_set_vec_vec_angle(dofaxis,tgt0,v1,v2):
	gets the rotation angles, if any, of v2 around dofaxis such that angle(v1,v2) == tgt0
	dofaxis = helix axis
	tgt0 = sym magic angle (54.7=asin(sr2/sr3),35.3=asin(sr1/sr3),22.9=asin(G/2/sr3))
	v1 = symaxis1 (comp1 and comp2 in arbitrary frame, but aligned together along helix axis)
	v2 = symaxis2 (comp1 and comp2 in arbitrary frame, but aligned together along helix axis)

you then "slide" along the helix axis vector to intersect the sym axis:
def slide_to_make_lines_intersect(dof,l,l0,m,m0):
	dof = helix axis
	l   = symaxis1
	l0  = any point along line of symmetry for comp1
	m   = symaxis2 (rotated according to the above around helix axis)
	m0  = any point along line of symmetry for comp2
"""

HELIX_D_PER_RES = 1.53663757571
HELIX_CA_RADIUS = 2.27024021944

def helix_ray(sele,nc):
	r = getres(sele)
	atom = cmd.get_model(sele+" and name CA").atom
	# print "".join(a.ss for a in atom)
	if any(a.ss!='H' for a in atom):
		# print "helix_ray: not all helix!"
		return None
	if len(atom) < 9:
		print "helix less than 10 res",len(atom)
		return None
	ca = [Vec(a.coord) for a in atom]
	begin = reduce(operator.add,ca[ :7])/7.0
	end   = reduce(operator.add,ca[-7:])/7.0
	axis  = (end-begin).normalized()
	if nc:
		ca0,cen = ca[+6],begin
	else:
		ca0,cen = ca[-6],end
	cen = cen+proj(axis,ca0-cen)
	vca = projperp(axis,ca0-cen).normalized()
	# showvecfrompoint(vca,cen)
	return axis,cen,vca
def c_helix_ray(sele): return helix_ray(sele,False)
def n_helix_ray(sele): return helix_ray(sele,True)
def n_terminal_helix(ss,maxend=5):
	"""
	>>> ss = "LLLHHHHHHHHHHHHHHHHLLLLLLLLLLLLLHHHHHHHHHHHHHHHHLLL"
	>>> hh = n_terminal_helix(ss)
	>>> print ss[:hh[0]],ss[hh[0]:hh[1]+1],ss[hh[1]+1:]
	LLLH HHHHHHHHHHHHHH HLLLLLLLLLLLLLHHHHHHHHHHHHHHHHLLL

	>>> hh = c_terminal_helix(ss)
	>>> print ss[:hh[0]],ss[hh[0]:hh[1]+1],ss[hh[1]+1:]
	LLLHHHHHHHHHHHHHHHHLLLLLLLLLLLLLH HHHHHHHHHHHHHH HLLL

	>>> print c_terminal_helix("")
	None

	>>> print c_terminal_helix("L")
	None

	>>> print c_terminal_helix("H")
	None
	"""
	i = 0
	while i < len(ss) and ss[i]=='L' and ss[i]!='E': i += 1
	nbeg = i
	while i < len(ss) and ss[i]=='H': i += 1
	if not i: return None
	nend = i-1 if ss[i-1]=='H' and i else i-2
	# s = "".join(a[1] for a in ss1)
	nbeg,nend = nbeg+1,nend-1
	if nbeg >= nend or nbeg > maxend: return None
	nend = min(nend,nbeg+13)
	return nbeg,nend
def c_terminal_helix(ss,maxend=5):
	if not ss: return None
	ss = list(reversed(ss))
	h = n_terminal_helix(ss,maxend)
	if not h: return None
	return len(ss)-h[1]-1,len(ss)-h[0]-1
def helix_termini_rays(sele,maxend=5):
	atom = cmd.get_model(sele+" and name ca").atom
	ss = "".join(a.ss for a in atom)
	res = [    a.resi for a in atom]
	nh = n_terminal_helix(ss,maxend)
	ch = c_terminal_helix(ss,maxend)
	nhsel = sele+" and resi "+"+".join(res[nh[0]:nh[1]+1]) if nh else "none"
	chsel = sele+" and resi "+"+".join(res[ch[0]:ch[1]+1]) if ch else "none"
	# print "helix ray '",sele,"' nh", nhsel, nh
	nr = n_helix_ray(nhsel) if nh else None
	# print "helix ray '",sele,"' ch", chsel, ch
	cr = c_helix_ray(chsel) if ch else None
	ni = [res[i] for i in nh] if nh else None
	ci = [res[i] for i in ch] if ch else None
	return ni,ci,nr,cr
def show_helix_termini_rays(sele,maxend=5):
	for r in helix_termini_rays(sele,maxend)[2:]:
		if r:
			showvecfrompoint(r[0]*10,r[1])
			showvecfrompoint(r[2]*05,r[1])

def find_helix_align_cage_xforms(nr1,cr1,nr2,cr2,nfst1,nfst2,num_sym_expand=1,showme=None):
	"""
	assumes sym axes are both Z
	"""
	symangs = (ATET ,AOCT ,AICS )
	symlabs = ("TET","OCT","ICS")
	# symangs = (45.0,)
	# symlabs = ("OCT",)
	if showme and showme.upper().startswith("NO"):
		showme = None
	# nr1,cr1 = helix_termini_rays(sele1)
	# nr2,cr2 = helix_termini_rays(sele2)
	xforms = []
	samples = 0
	for swap in (True,False):
	# for sel1,sel2,nf1,nf2,cr,nr in ((sele1,sele2,nf1st,nf2st,cr1,nr2),(sele2,sele1,nf2st,nf1st,cr2,nr1)):
		# sel1,sel2 = (sele2,sele1) if swap else (sele1,sele2)
		nf1,nf2   = (nfst2,nfst1) if swap else (nfst1,nfst2)
		cr,nr     = (cr2,nr1)     if swap else (cr1,nr2)
		if not cr or not nr: continue
		# print "JOIN",cr
		# print "  TO",nr
		a1,c1,r1 = cr
		a2,c2,r2 = nr
		symaxs1 = Uz
		symcen1 = V0

		# showvecfrompoint(a1*10,c1,lbl="a1")
		# showvecfrompoint(r1*05,c1,lbl="r1")
		# showvecfrompoint(a2*10,c2,lbl="a2")
		# showvecfrompoint(r2*05,c2,lbl="r2")
		# showvecfrompoint(symaxs1*20,V0,lbl="symaxs1")
		x = alignvectors(a2,r2,a1,r1)
		x.t = c1-x*c2
		s2 = x.R*Uz
		assert a1.angle_degrees(x.R*a2) < 0.1
		for symlab,symang in izip(symlabs,symangs):
			# print symlab
			for a in rotation_around_dof_to_set_vec_vec_angle(a1,symang,symaxs1,s2):
				samples += 1
				print "found solution for",symlab,symang
			# align helix and symmetry
				x2 = rotation_around_degrees(a1,a,c1)
				symaxs2 = x2.R*s2
				symcen2 = x2*x*V0
				assert line_line_angle_degrees(symaxs1,symaxs2)-symang < 0.1

				d = slide_to_make_lines_intersect(a1,symaxs2,symcen2,symaxs1,V0)
				if not 0 < d < 40: continue
				# print symlab,d,swap
				x2.t += a1*d
				symcen2 += a1*d
				sc = skew_lines_center(symaxs1,V0,symaxs2,symcen2)
				ax1 = Xform(Imat,-sc)
				ax2 = ax1*x2*x
				# if min(ax1.t.length(),ax2.t.length())/sindeg(symang) < 100: continue
				# if /sindeg(symang) < 100: continue
				# print ax1.t.length(),ax1.t.length()/sindeg(symang)
				# print ax2.t.length(),ax2.t.length()/sindeg(symang)
				tgtsymaxs1 = SYMAXIS[symlab,nf1]
				tgtsymaxs2 = SYMAXIS[symlab,nf2]

				assert line_line_angle_degrees(tgtsymaxs1,tgtsymaxs2)-symang < 0.1
				assert line_line_angle_degrees( ax1.R*Uz , ax2.R*Uz )-symang < 0.1
				# print line_line_angle_degrees(tgtsymaxs1,tgtsymaxs2)
				# print line_line_angle_degrees(ax1.R*Uz,ax2.R*Uz)
				# print ax1.R*Uz,"tgt",tgtsymaxs1
				# print ax2.R*Uz,"tgt",tgtsymaxs2

				orig1,orig2 = ax1.R*Uz, ax2.R*Uz
				# if orig1.dot(tgtsymaxs1) < 0.0: orig1 = -orig1
				# if orig2.dot(tgtsymaxs2) < 0.0: orig2 = -orig2
				if orig1.dot(orig2) < 0: orig2 = -orig2
				xsymaln = alignvectors_minangle(orig1,orig2,tgtsymaxs1,tgtsymaxs2)
				# print xsymaln.R
				ax1 = xsymaln*ax1
				ax2 = xsymaln*ax2
				assert line_line_angle_degrees(ax1.R*Uz,ax2.R*Uz)-symang < 0.1

				# print ax1.R*Uz,"tgt",tgtsymaxs1
				# print ax2.R*Uz,"tgt",tgtsymaxs2
				assert line_line_angle_degrees(ax1.R*Uz,tgtsymaxs1) < 0.1
				assert line_line_angle_degrees(ax2.R*Uz,tgtsymaxs2) < 0.1

				assert (ax1.R*a1).angle_degrees(ax2.R*a2) < 0.1

				# helix alignment
				out_of_plane = abs((ax1.R*Uz).cross(ax2.R*Uz).normalized().dot((ax1.R*a1).normalized()))
				# if out_of_plane > 0.5: continue
				# print "plane",out_of_plane

				alnreg = dihedral_degrees(ax1.R*r1,V0,ax1.R*a1,ax2.R*r2)
				tgtreg = (d*66.66+180)%360-180
				# print alnreg,tgtreg
				regerr = min(abs(tgtreg-alnreg),abs(tgtreg+360-alnreg),abs(tgtreg-alnreg-360))
				diserr = 0.75-abs( d % 1.5 - 0.75 )
				# print d,d%1.5-.75
				distol = ( 0.2+0.01*d)/out_of_plane
				regtol = (20.0+1.00*d)/out_of_plane
				# print "D",d,"REG",regerr,regtol,"DIS",diserr,distol
				# print (diserr)**2 + (regerr/100)**2 , (regtol/100)**2+distol**2
				# if (diserr)**2 + (regerr/100)**2 > (regtol/100)**2+distol**2: continue
				if diserr > distol: continue
				if regerr > regtol: continue

				# clashes
				symcentmp1,hlxcentmp1 = ax1*V0,ax1*c1
				symcentmp2,hlxcentmp2 = ax2*V0,ax2*c2
				sinsymang = sindeg(symang)
				# print
				# print "rough clash symcen1",symcentmp1.length()
				# print "rough clash symcen2",symcentmp2.length()
				# print "rough clash hlxcen1",hlxcentmp1.length()
				# print "rough clash hlxcen2",hlxcentmp2.length()
				# print "rough clash symcend",(symcentmp1-symcentmp2).length()
				# print "rough clash symdot ",symcentmp1.dot(symcentmp2)
				fail = []
				if symcentmp1.length()              < 10/1/sinsymang: continue #fail.append("cen1")
				if symcentmp2.length()              < 10/1/sinsymang: continue #fail.append("cen2")
				if (symcentmp1-symcentmp2).length() < 10/2/sinsymang: continue #fail.append("cend")
				if hlxcentmp1.length()              < 10/1/sinsymang: continue #fail.append("hlx1")
				if hlxcentmp2.length()              < 10/1/sinsymang: continue #fail.append("hlx2")
				if symcentmp1.dot(symcentmp2) < 0.0:                  continue
				if fail:
					continue
					# print "fail:",fail

				if showme=="ALL" or showme==symlab:
					tag = str(random.random())[2:]
					o1,o2 = 'ha_'+tag+'_1','ha_'+tag+'_2'
					cmd.create(o1,'comp2' if swap else 'comp1')
					cmd.create(o2,'comp1' if swap else 'comp2')
					cmd.alter(o1,"chain='A'")
					cmd.alter(o2,"chain='B'")
					# showvecfrompoint(ax1.R*r1*10,ax1*c1,lbl="H1REG")
					# showvecfrompoint(ax2.R*r2*10,ax2*c2,lbl="H2REG")
					# xform(o1,ax1)
					# xform(o2,ax2)
					showvecfrompoint(ax1.R*Uz*10,ax1*V0,lbl="SYM1")
					showvecfrompoint(ax2.R*Uz*10,ax2*V0,lbl="SYM2")

			# gen sym frames
				G1 = rotation_around_degrees(ax1.R*Uz,360.0/nf1,V0)
				G2 = rotation_around_degrees(ax2.R*Uz,360.0/nf2,V0)
				gx = []
				if num_sym_expand < 2:
					if swap: gx.append((ax2,ax1))
					else:    gx.append((ax1,ax2))
					if showme=="ALL" or showme==symlab: print "cant 'showme' with num_sym_expand < 2"
				else:
					for i,xs in enumerate(expand_xforms([G1,G2],num_sym_expand)):
						xs1,xs2 = xs*ax1,xs*ax2
						if swap: gx.append((xs2,xs1))
						else:    gx.append((xs1,xs2))
						if showme=="ALL" or showme==symlab:
							symaxssubunit1,symorisubunit1,hlxdirsubunit1,hlxregsubunit1 = (xs1.R*v for v in (Uz,Uy,a1,r1))
							symaxssubunit2,symorisubunit2,hlxdirsubunit2,hlxregsubunit2 = (xs2.R*v for v in (Uz,Uy,a2,r2))
							symcensubunit1,hlxcensubunit1                               = (xs1  *v for v in (V0,c1))
							symcensubunit2,hlxcensubunit2                               = (xs2  *v for v in (V0,c2))
							cmd.create("sg1%i"%i,o1)
							cmd.create("sg2%i"%i,o2)
							xform("sg1%i"%i,xs1)
							xform("sg2%i"%i,xs2)
							if i==0:
								showvecfrompoint(symaxssubunit1*20,symcensubunit1,(1,0,0),lbl="sg1%isym" %i)
								#showvecfrompoint(symorisubunit1*10,symcensubunit1,(1,0,0),lbl="sg1%iori" %i)
								showvecfrompoint(hlxdirsubunit1*20,hlxcensubunit1,(1,0,0),lbl="sg1%ihdir"%i)
								#showvecfrompoint(hlxregsubunit1*10,hlxcensubunit1,(1,0,0),lbl="sg1%ihreg"%i)
								showvecfrompoint(symaxssubunit2*20,symcensubunit2,(0,0,1),lbl="sg2%isym" %i)
								#showvecfrompoint(symorisubunit2*10,symcensubunit2,(0,0,1),lbl="sg2%iori" %i)
								showvecfrompoint(hlxdirsubunit2*20,hlxcensubunit2,(0,0,1),lbl="sg2%ihdir"%i)
								#showvecfrompoint(hlxregsubunit2*10,hlxcensubunit2,(0,0,1),lbl="sg2%ihreg"%i)
				xforms.append((symlab,swap,gx))

				if showme=="ALL" or showme==symlab:
					print "SHOWING",symlab,samples
					# print symang,(ax1.R*Uz).angle_degrees(ax1.R*Uz)
					# print nf1,tgtsymaxs1,ax1.R*Uz
					# print nf2,tgtsymaxs2,ax2.R*Uz
					raise NotImplementedError
					cmd.delete(o1)
					cmd.delete(o2)

	return xforms,samples

def find_helix_align_cages(sele1,sele2,nf1,nf2):
	ni1,ni2,nr1,cr1 = helix_termini_rays(sele1)
	ni1,ni2,nr2,cr2 = helix_termini_rays(sele2)
	for i,hits in enumerate(find_helix_align_cage_xforms(nr1,cr1,nr2,cr2,nf1,nf2)):
		nsamp,xs = hits
		for j,x in enumerate(xs):
			sym,swap,x1,x2 = x
			cmd.create("test_%i_%i_A"%(i,j),sele1+" and name n+ca+c+o")
			cmd.create("test_%i_%i_B"%(i,j),sele2+" and name n+ca+c+o")
			xform("test_%i_%i_A"%(i,j),x1)
			xform("test_%i_%i_B"%(i,j),x2)
			r1,r2 = (nr1,cr2) if swap else (cr1,nr2)
			showvecfrompoint(x1.R*r1[0]*9,x1*r1[1],(0,1,0))
			showvecfrompoint(x1.R*r1[2]*4,x1*r1[1],(0,1,0))
			showline(x1.R*Uz*50  ,x1*V0   ,(0,1,0))
			showvecfrompoint(x2.R*r2[0]*9,x2*r2[1],(0,0,1))
			showvecfrompoint(x2.R*r2[2]*4,x2*r2[1],(0,0,1))
			showline(x2.R*Uz*50  ,x2*V0   ,(0,0,1))

def showhit(h,outfile):
	f1,f2,tmp = h
	sym,swap,Xs = tmp
	# print f1,f2
	# print h
	# cmd.delete('all')
	# cmd.load(f1,"comp1")
	# cmd.load(f2,"comp2")
	# cmd.alter("comp1","chain='A'")
	# cmd.alter("comp2","chain='B'")
	for i,x in enumerate(Xs):
		cmd.create("cage_%i_A"%i,"comp1")
		cmd.create("cage_%i_B"%i,"comp2")
		cmd.alter("cage_%i_A"%i,"chain='%s'"%alphabet[2*i+0])
		cmd.alter("cage_%i_B"%i,"chain='%s'"%alphabet[2*i+1])
		# print x
		if swap: x0,x1 = x[1],x[0]
		else:    x0,x1 = x[0],x[1]
		xform("cage_%i_A"%i,x0)
		xform("cage_%i_B"%i,x1)
	# print "saving",outfile
	cmd.save(outfile,"cage_*")

def read_hdat(l):
	# print l
	s = l.split()
	pdb = s[0]
	nres = int(s[1])
	s.pop(1)
	h1b,h1e = int(s[1]),int(s[2])
	h1a = Vec(float(s[ 3]),float(s[ 4]),float(s[ 5]))
	h1c = Vec(float(s[ 6]),float(s[ 7]),float(s[ 8]))
	h1r = Vec(float(s[ 9]),float(s[10]),float(s[11]))
	h2b,h2e = int(s[12]),int(s[13])
	h2a = Vec(float(s[14]),float(s[15]),float(s[16]))
	h2c = Vec(float(s[17]),float(s[18]),float(s[19]))
	h2r = Vec(float(s[20]),float(s[21]),float(s[22]))
	return (pdb,nres,(h1b,h1e),(h1a,h1c,h1r),(h2b,h2e),(h2a,h2c,h2r))

def read_helix_term_data(fn,maxres=999999999):
	with open(fn) as in1:
		L = in1.readlines()
		# random.shuffle(L)
		print fn,"num temini (N and/or C)",len(L)
		hdat = (read_hdat(l) for l in L)
		return [x for x in hdat if x[1] < maxres]

def find_helix_align_cages_from_datafiles(d1,d2,nf1,nf2,N=999999999999,outdir=None,num_sym_expand=0,showme="NONE",maxres=100):
	# PDBs = set(s.strip() for s in open("/Users/sheffler/dropbox/project/pdbcatalog/symcmp_c2_100.txt").readlines())
	# global hdat1,hdat2
	# if not hdat1: hdat1 = read_helix_term_data(d1)
	# if not hdat2: hdat2 = read_helix_term_data(d2)
	hdat1 = read_helix_term_data(d1,maxres)
	hdat2 = read_helix_term_data(d2,maxres)
	hits = []
	count,tries,npdbs = 0,0,0
	for pa,nresa,ba1,ra1,ba2,ra2 in hdat1:
		# if os.path.basename(pa)[:4] not in PDBs: continue
		if max(ba1) is 0: ba1,ra1=None,None
		if max(ba2) is 0: ba2,ra2=None,None
		if outdir or (showme and showme.upper()!="NONE"):
			# print showme
			cmd.delete("comp1")
			cmd.load(pa,"comp1")
			cmd.hide('ev')
			cmd.remove('comp1 and not name n+ca+c+o and not resn A+T+G+C')
			cmd.alter("comp1","chain='A'")
			# if ra1: showvecfrompoint(ra1[0]*10,ra1[1],lbl='raa1')
			# if ra1: showvecfrompoint(ra1[2]*05,ra1[1],lbl='rar1')
			# if ra2: showvecfrompoint(ra2[0]*10,ra2[1],lbl='raa2')
			# if ra2: showvecfrompoint(ra2[2]*05,ra2[1],lbl='rar2')
		for pb,nresb,bb1,rb1,bb2,rb2 in hdat2:
			if max(bb1) is 0: bb1,rb1=None,None
			if max(bb2) is 0: bb2,rb2=None,None
			if (ba1 and bb2) or (bb1 and ba2):
				npdbs += 1
				if outdir or (showme and showme.upper()!="NONE"):
					cmd.delete("comp2")
					cmd.load(pb,"comp2")
					cmd.hide('ev')
					cmd.remove('comp2 and not name n+ca+c+o and not resn A+T+G+C')
					cmd.alter("comp2","chain='B'")
					# raise NotImplementedError
					# if rb1: showvecfrompoint(rb1[0]*10,rb1[1],lbl='rba1')
					# if rb1: showvecfrompoint(rb1[2]*05,rb1[1],lbl='rbr1')
					# if rb2: showvecfrompoint(rb2[0]*10,rb2[1],lbl='rba2')
					# if rb2: showvecfrompoint(rb2[2]*05,rb2[1],lbl='rbr2')
				# print "ra2",ra2
				# print "rb2",rb2
				sx,nsamp1 = find_helix_align_cage_xforms(ra1,ra2,rb1,rb2,nf1,nf2,num_sym_expand=num_sym_expand,showme=showme)
				tries += nsamp1
				for isx,x in enumerate(sx):
					count += 1
					sym,swap,xforms = x
					fn = "FHACFD_%s_%s_%s_%i.pdb"%(os.path.basename(pa).split("_")[0],os.path.basename(pb).split("_")[0],sym,isx)
					print "HIT:",fn,swap,x[0],x[1],len(x[2]),x[2][0][0].R.xx
					hits.append((pa,pb,x))
					# if count > 0:
					if outdir or (showme and showme.upper()!="NONE"):
						print "HIT!"
						showhit(hits[-1],outdir+fn)
						# return
					if count >= N: return count,tries,npdbs,len(hdat1)*len(hdat2)*count/npdbs
	f = len(hdat1)*len(hdat2)*count/npdbs if npdbs else None
	return count,tries,npdbs,f

# collect_helix_termini("/Users/sheffler/project/struts/input/josh_dna/*.pdb","/Users/sheffler/project/struts/input/josh_dna/josh_dna_C3.helixterm")
def collect_helix_termini(pattern,outfile):
	cmd.delete("cht_*")
	with open(outfile,'w') as out:
		for i,fn in enumerate(glob.glob(pattern)):
			print i,fn
			tag = "cht_"+str(random.random())
			cmd.load(fn,tag)
			nres = cmd.select(tag+" and name ca")
			ni,ci,nr,cr = helix_termini_rays(tag,8)
			# print fn,nr,cr
			if nr or cr:
				out.write(fn)
				out.write("%5i"%nres)
				for i,r in ((ni,nr),(ci,cr)):
					if not r: i,r = (0,0),(V0,V0,V0)
					out.write(" %6i %6i "%(int(i[0]),int(i[1])))
					out.write((r"%10.5f "*9)%(r[0].x,r[0].y,r[0].z,r[1].x,r[1].y,r[1].z,r[2].x,r[2].y,r[2].z))
				out.write('\n')

			# cmd.show('car')
			# util.chainbow()
			# if nr: showvecfrompoint(nr[0]*9,nr[1],(0,1,0))
			# if nr: showvecfrompoint(nr[2]*4,nr[1],(0,1,0))
			# if cr: showvecfrompoint(cr[0]*9,cr[1],(0,0,1))
			# if cr: showvecfrompoint(cr[2]*4,cr[1],(0,0,1))
			cmd.delete(tag)

def teststrut(N=10,maxres=99999,showme="NONE",nexpand=1,outdir="/work/sheffler/project/struts/"):
	outdir += "/"
	print "teststrut","N",N,"maxres",maxres,"showme",showme,"nexpand",nexpand,"outdir",outdir

	cmd.delete('all')
	infile1 = "/Users/sheffler/project/struts/input/C2.helixterm"
	# infile2 = "/Users/sheffler/project/struts/input/C3.helixterm"
	infile2 = "/Users/sheffler/project/struts/input/josh_dna/josh_dna_C3.helixterm"
	# infile2 = "/Users/sheffler/project/struts/input/jorge_C4.dat"
	if outdir:
		for f in glob.glob(outdir+"FHACFD_*.pdb"): os.remove(f)
	count,tries,npdbs,esthits = find_helix_align_cages_from_datafiles(infile1,infile2,2,3,N,outdir,num_sym_expand=nexpand,showme=showme,maxres=maxres)
	print "stopped after",count,"hits,",npdbs,"pdbs,",tries,"tries  EST HITS:", str(esthits)

	print """
	set cartoon_cylindrical_helices=1
	run /Users/sheffler/pymol/struts.py
	cProfile.runctx("test()",globals(),locals(),"/work/sheffler/tmp/test.prof")
	p = pstats.Stats('/work/sheffler/tmp/test.prof')
	p.sort_stats("cumulative").print_stats(30)
	os.remove("/work/sheffler/tmp/test.prof")
	"""


	# v = cmd.get_view()
	# print "TEST HELIX ALIGN"
	# cmd.delete('line*')
	# Haxis = Vec(0,0,1)
	# # showline(Haxis*100,Vec(0,0,0))
	# a1 = randnorm()
	# c1 = randvec()*2
	# a2 = randnorm()
	# c2 = randvec()*2
	# showline(a1*100,c1,col=(1,0,1))
	# for tgt in (ATET,AOCT,AICS):
	# 	for ang in rotation_around_dof_to_set_vec_vec_angle(Haxis,tgt,a1,a2):
	# 		rot = rotation_around_degrees(Haxis,ang,Vec(0,0,0))
	# 		a2rot = rot*a2
	# 		c2rot = rot*c2
	# 		atest = angle_degrees(a1,Vec(0,0,0),a2rot)
	# 		if abs(tgt-atest) > 0.000001 and abs(180-tgt-atest) > 0.000001:
	# 			print "ERROR!!!!!!!!!!!!!"
	# 			continue
	# 		d = slide_to_make_lines_intersect(Haxis,a2rot,c2rot,a1,c1)
	# 		anew = rot*a2
	# 		cnew = rot*c2 + d*Haxis
	# 		print "%-07.3f %-07.3f %-07.3f %-07.3f "%(atest,ang,d,atest)
	# 		showline(anew*10,cnew)
	# cmd.set_view(v)

def testhelix(filepattern="/work/sheffler/sym_comp/C2/*.pdb.gz",minlen=10,sele="all"):
	files = glob.glob(filepattern)
	random.shuffle(files)
	for fn in files[:1]:
		print fn
		cmd.delete('all')
	# for fn in glob.glob("/work/sheffler/project/struts/doc/ideal_alpha.pdb"):
		cmd.load(fn,"test")
		cmd.remove('not name n+ca+c+o and not resn A+T+G+C')
		hinfos = list()
		# atoms = cmd.get_model("not byres (all within 5 of not byres ss h)").atom
		atoms = cmd.get_model("(("+sele+") and (byres ss h))").atom
		# for a in atoms: print a.resi,a.name
		halocal = Vec( 0.889193, 0.344166, 0.301472).normalized()
		hclocal = Vec(-0.327422,-0.862082, 2.074854)
		avg = Vec()
		lastir = 999999
		for ir,n,ca,c,o in ([int(atoms[i].resi),]+[Vec(a.coord) for a in atoms[i:i+4]] for i in range(0,len(atoms),4)):
			if ir != lastir+1: hinfos.append(list())
			lastir = ir
			s = stub(ca,o,n,c)
			axs,cen = s.R*halocal,s*hclocal
			cen -= proj(axs,cen-ca)
			reg = (ca-cen).normalized()
			# print ca.distance(cen)
			hinfos[-1].append((ir,axs,cen,reg))
			# showvecfrompoint(s.R*Ux,ca,(1,0,0))
			# showvecfrompoint(s.R*Uy,ca,(0,1,0))
			# showvecfrompoint(s.R*Uz,ca,(0,0,1))
			# showvecfrompoint((cen-ca)*0.9,ca,(1,1,1))
			# print ir,s.tolocal(Vec(-29.978799,65.021254,-26.611610))
			# avg += s.R.transposed()*Vec(-0.428876,0.840627,-0.330773)
		hinfos = [hinfo for hinfo in hinfos if len(hinfo)-2 >= minlen]
		# print avg.normalized()
		# aa = Vec()
		# ac = Vec()
		for ihinfo,hinfo in enumerate(hinfos):
			e1,e2,e3=0,0,0
			for i in range(1,len(hinfo)-1):
				pi,pa,pc,pr = hinfo[i-1]
				ti,ta,tc,tr = hinfo[i+0]
				ni,na,nc,nr = hinfo[i+1]
				assert pi+1==ti and ti==ni-1
				e1single = (1-(nc-pc).normalized().dot( ta  )   **2)*HELIX_D_PER_RES
				e2single = (1-(      ta           .dot(pa+na)/2)**2)*HELIX_D_PER_RES
				e3single = tc.distance_squared((pc+nc)/2)
				e1 += e1single
				e2 += e2single
				e3 += e3single
				# aa += ta
				# ac += tc
				# showsegment(pc,tc)
				# showsegment(tc,nc)
				# showvecfrompoint(ta*HELIX_D_PER_RES,tc)
				showvecfrompoint(ta*+(HELIX_D_PER_RES/2),tc)
				showvecfrompoint(ta*-(HELIX_D_PER_RES/2),tc)
				showvecfrompoint(tr*+(HELIX_CA_RADIUS  ),tc)
				# showvecfrompoint(-ta,tc)
				# showvecfrompoint(hinfo[i][0]*100,hinfo[i][1])
			print ihinfo,"%9.5f"%(sqrt(e1/(len(hinfo)-2)))
			print ihinfo,"%9.5f"%(sqrt(e2/(len(hinfo)-2)))
			print ihinfo,"%9.5f"%(sqrt(e3/(len(hinfo)-2)))
		# ac /= (len(hinfo)-2)
		# print aa.normalized()
		# print ac
		# showvecfrompoint(aa*50,ac)
		# showvecfrompoint(aa*-50,ac)
		# a = (hinfo[-1][1]-hinfo[1][1]).normalized()
		# print a
		# tgt = Vec(atoms[-3].coord) - projperp(a,Vec(atoms[-3].coord)-(hinfo[-1][1]+hinfo[1][1])/2)
		# print tgt
		# showvecfrompoint(-a*100,tgt)




if __name__ == '__main__':
	# xforms = read_xforms("/Users/sheffler/project/symgen_movie/fere.pdb_xforms.dat")
	# print len(xforms)
	# print xforms[-1]

	import doctest
	r = doctest.testmod()
	print r

	# find_helix_align_cages_from_datafiles("/Users/sheffler/dropbox/project/struts/input/C2.helixterm","/Users/sheffler/dropbox/project/struts/input/C3.helixterm",2,3)

