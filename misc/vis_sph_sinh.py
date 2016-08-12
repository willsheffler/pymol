# https://github.com/cffk/orientation/blob/v1.1/ExpandSet.cpp

from xyzMath import Vec
from math import sinh, asinh

max_sinh = 0.0
max_asinh = 0.0


def myasinh(x):
	#x-1/6x^3+3/(40)x^5-5/(112)x^7+(35)/(1152)x^9+...
	x2 = x*x
	x3 = x2*x
	x5 = x2*x3
	x7 = x2*x5
	x9 = x2*x7
	return x - x3/6.0 + x5*3.0/40.0 - x7*5.0/112.0 + x9*35.0/1152.0

def mysinh(x):
	# return sinh(x)
	x2 = x*x
	x3 = x2*x
	x5 = x2*x3
	x7 = x2*x5
	x9 = x2*x7
	v = x + x3/6.0 #+ x5/120.0 #+ x7/5040.0 + x9/362880.0
	global max_sinh, max_asinh
	if abs(x) > max_sinh:
		max_sinh = abs(x)
	if abs(v) > max_asinh:
		max_asinh = abs(v)
	return v


def hexgrid(n, w, maxang, sigma):
	global max_sinh, max_asinh
	dotthresh = cos(radians(maxang))
	# if sigma != 0.0:
		# w = w / sigma
	X = w*Vec(1,0,0)
	Y = w*Vec(0.5,sqrt(3.0)/2,0)
	points = []
	for i in range(-n,n+1):
		for j in range(-n,n+1):
			p = i*X + j*Y + Vec(0,0,1)
			old = max_sinh, max_asinh
			if sigma != 0.0:
				p.x = mysinh(p.x*sigma)/sigma
				p.y = mysinh(p.y*sigma)/sigma
			p2 = p
			p = p / p.length()
			if p.z >= dotthresh:
				points.append((p,p2))
			else:
				# restore if not accepted
				max_sinh, max_asinh = old
	return points

# for 60°, this looks nice-ish...
# PyMOL>run /Users/sheffler/pymol/misc/vis_sph_sinh.py; testhex(radius=10, spacing=0.1, sigma=3,maxang=60, scale=50, n=15)
def testhex(radius=5.0, spacing=0.5, sigma=0.0, scale=20.0, n=30, maxang=60.0):
	v = cmd.get_view()
	cmd.delete('hextest')
	cgo = []
	for p,p2 in hexgrid(n, spacing, maxang, sigma):
		p = scale * p
		p2 = scale * p2
		cgo.extend(cgo_sphere(p))
		cgo.extend(cgo_segment(p,p2))
		cgo.extend(cgo_sphere(p2))
	cmd.load_cgo(cgo,'hextest')
	cmd.set_view(v)
	print "max used val to sinh:", max_sinh, ", max used val from sinh:",max_asinh

def test_sinh_inv():
	for x in range(1,17):
		x = x/2.0
		print x, asinh(x), myasinh(x)/asinh(x)
