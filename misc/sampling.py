from pymol_util import *
from math import *

def drawtetra():
	cmd.delete("*")
	N=4
	s = 1
	h = sqrt(3.0)/2.0
	count = 1
	P = list()
	E = list()
	with open("/work/sheffler/tmp/hex.pdb",'w') as o:
		for ix,iy,iz in [(float(x),float(y),float(z)) for x in range(N+0) for y in range(N+1) for z in range(N+1)]:
			x = ix+iy%2/2+iz%2/2
			y = iy*h+iz%2*h/3
			z = iz*h
			P.append(Vec(x,y,z))
	for p in P:
		for q in P:
			if q is p: continue
			if q.distance(p) < 1.1:
				E.append((p,q))
	for p in P:
		showsphere(p,r=0.1,col=(1.0-p.z/h/N,((p.z+1)/h/N)%1.0,p.z/h/N))
	for p,q in E:
		showlineabs(p,q)
		# print >>o, "ATOM %6i CUBE CUB %5i       %03.3f   %03.3f   %03.3f  1.00  1.00"%(count,count,x,y,z)
		# count += 1

def drawcells(E,col=(1,1,1),CR=0.53):
	cgo = []
	for p,q in E:
		dis = p.distance(q)
		print dis
		if dis < 0.01: continue
		cen = (p+q)/2.0
		nrm = (p-q).normalized()
		cgo += cgo_cyl( c1 = cen+nrm/500, c2 = cen-nrm/500, r = sqrt(CR*CR-(dis*dis/4)), col=col )
	return cgo

def drawbcc(N=2,scale=1.0,col=(1,1,1)):
	cgo = []
	P,E = list(),list()
	# with open("/work/sheffler/tmp/hex.pdb",'w') as o:
	for ix,iy,iz in [(float(x),float(y),float(z)) for x in range(N+1) for y in range(N+1) for z in range(N+1)]:
		P.append(scale*(Vec(ix+0.0,iy+0.0,iz+0.0)-Vec(N,N,N)/2.0))
		if ix==N or iy==N or iz==N: continue
		P.append(scale*(Vec(ix+0.5,iy+0.5,iz+0.5)-Vec(N,N,N)/2.0))
	# cen = Vec(0,0,0)
	# for p in P: cen += p
	# cen /= len(P)
	# print cen
	# for i,p in enumerate(P): P[i] = p-cen
	# for p in P: showsphere(p,r=0.458*scale,col=col)
	for p in P:
		cgo += cgo_sphere(p,r=0.1,col=col)
	print N, len(P)
	E = [(p,q) for p in P for q in P if p is not q and p.distance(q) < 1.1*scale]
	cgo += drawcells(E ,col=col,CR=0.54*scale)

	# for p in P:
	# 	R = 0.2 if p.z%1==0 else 0.15
	# 	showsphere(p,r=R,col=(1.0-p.z/h/N,((p.z+1)/h/N)%1.0,p.z/h/N)) # cov rad .56
	# 	showsphere(p,r=0.1,col=(1.0,1.0,1.0))
	# 	pass
	# for p,q in E:
	# 	# dis = p.distance(q)
	# 	# if 0.99 < dis < 1.01: showlineabs(p,q,col=(0.3,0.3,0.3))
	# 	cen = (p+q)/2.0
	# 	P2.append(cen)
	# 	# nrm = (p-q).normalized()
	# 	# R = 0.1 if 0.99<dis<1.01 else 0.066
	# 	# showsphere(cen,r=R)
	# 	# print >>o, "ATOM %6i CUBE CUB %5i       %03.3f   %03.3f   %03.3f  1.00  1.00"%(count,count,x,y,z)
	# 	# count += 1
	# E2 = [(p,q) for p in P2 for q in P2 if p is not q and p.distance(q) < 0.55]
	# drawcells(E2,col=(0.0,0.0,1.0),CR=0.27)
	cmd.load_cgo(cgo,"bcc"+str(N))


def drawfcc():
	cmd.delete("*")
	N=2
	for x,y,z in [(float(x),float(y),float(z)) for x in range(N+1) for y in range(N+1) for z in range(N+1)]:
		if x==0: showlineabs(Vec(0,y,z),Vec(N,y,z),col=(0.9,0.9,0.9))
		if y==0: showlineabs(Vec(x,0,z),Vec(x,N,z),col=(0.9,0.9,0.9))
		if z==0: showlineabs(Vec(x,y,0),Vec(x,y,N),col=(0.9,0.9,0.9))
		showsphere(Vec(x,y,z),0.06,col=(1.0,1.0,1.0))
		if x!=N and y!=N:
			showsphere(            Vec(x+0.5,y+0.5,z+0.0),col=(1,0.5,0.5),r=0.06)
			showline(Vec(0,1,+1)/2,Vec(x+0.5,y+0.5,z+0.0),col=(0.2,0.2,0.2))
			showline(Vec(0,1,-1)/2,Vec(x+0.5,y+0.5,z+0.0),col=(0.2,0.2,0.2))
			showline(Vec(1,0,+1)/2,Vec(x+0.5,y+0.5,z+0.0),col=(0.2,0.2,0.2))
			showline(Vec(1,0,-1)/2,Vec(x+0.5,y+0.5,z+0.0),col=(0.2,0.2,0.2))
			showline(Vec(1,+1,0)/2,Vec(x+0.5,y+0.5,z+0.0),col=(0.2,0.2,0.2))
			showline(Vec(1,-1,0)/2,Vec(x+0.5,y+0.5,z+0.0),col=(0.2,0.2,0.2))
		if z!=N and x!=N:
			showsphere(            Vec(x+0.5,y+0.0,z+0.5),col=(1,0.5,0.5),r=0.06)
			showline(Vec(0,1,+1)/2,Vec(x+0.5,y+0.0,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(0,1,-1)/2,Vec(x+0.5,y+0.0,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(1,0,+1)/2,Vec(x+0.5,y+0.0,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(1,0,-1)/2,Vec(x+0.5,y+0.0,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(1,+1,0)/2,Vec(x+0.5,y+0.0,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(1,-1,0)/2,Vec(x+0.5,y+0.0,z+0.5),col=(0.2,0.2,0.2))
		if y!=N and z!=N:
			showsphere(            Vec(x+0.0,y+0.5,z+0.5),col=(1,0.5,0.5),r=0.06)
			showline(Vec(0,1,+1)/2,Vec(x+0.0,y+0.5,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(0,1,-1)/2,Vec(x+0.0,y+0.5,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(1,0,+1)/2,Vec(x+0.0,y+0.5,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(1,0,-1)/2,Vec(x+0.0,y+0.5,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(1,+1,0)/2,Vec(x+0.0,y+0.5,z+0.5),col=(0.2,0.2,0.2))
			showline(Vec(1,-1,0)/2,Vec(x+0.0,y+0.5,z+0.5),col=(0.2,0.2,0.2))


cmd.delete('all')
drawbcc(N= 1,scale=8.0,col=(0,1,0))
drawbcc(N= 2,scale=4.0,col=(0,0,1))
# drawbcc(N= 4,scale=2.0,col=(0,1,1))
#drawbcc(N= 8,scale=1.0,col=(1,1,1))
#drawbcc(N=16,scale=0.625,col=(1,1,0))
#drawbcc(N=32,scale=0.625,col=(1,1,0))
