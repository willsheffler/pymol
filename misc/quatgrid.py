import itertools

def cube(lb=Vec(-100.0,-100.0,-100.0),ub=Vec(100.0,100.0,100.0),r=0.5,xform=Xform()):
	v = cmd.get_view()
	lb = xform*lb
	ub = xform*ub
	cmd.load_cgo([
		cgo.CYLINDER, ub.x, ub.y, ub.z   , ub.x, ub.y, lb.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, ub.x, ub.y, lb.z   , ub.x, lb.y, lb.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, ub.x, lb.y, lb.z   , ub.x, lb.y, ub.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, ub.x, lb.y, ub.z   , ub.x, ub.y, ub.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, lb.x, ub.y, ub.z   , lb.x, ub.y, lb.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, lb.x, ub.y, lb.z   , lb.x, lb.y, lb.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, lb.x, lb.y, lb.z   , lb.x, lb.y, ub.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, lb.x, lb.y, ub.z   , lb.x, ub.y, ub.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, lb.x, ub.y, ub.z   , ub.x, ub.y, ub.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, lb.x, ub.y, lb.z   , ub.x, ub.y, lb.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, lb.x, lb.y, ub.z   , ub.x, lb.y, ub.z, r, 1,1,1, 1,1,1,
		cgo.CYLINDER, lb.x, lb.y, lb.z   , ub.x, lb.y, lb.z, r, 1,1,1, 1,1,1,
	],"UNIT_CELL")
	cmd.set_view(v)

def sign(x):
	return 1.0 if x > 0.0 else -1.0

def test(N=10):
	v = cmd.get_view()
	cgo = []
	for i,j,k in itertools.product(range(N),range(N),range(N)):
		p = Vec(i,j,k)*2.0/N-Vec(1.0-1.0/N)
		p.x = 0.8*p.x + 0.2*sign(p.x)*p.x**2
		p.y = 0.8*p.y + 0.2*sign(p.y)*p.y**2
		p.z = 0.8*p.z + 0.2*sign(p.z)*p.z**2
		cgo.extend( cgo_sphere(100*p) )
	cmd.load_cgo(cgo,"test")
	cube()
	cmd.set_view(v)
