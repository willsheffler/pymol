import math,random

class Vec(object):
   def __init__(self,x,y=None,z=None):
      if type(x) is type(self):
         self.x,self.y,self.z = x.x,x.y,x.z
      elif type(x) in (type([]),type((1,))):
         self.x,self.y,self.z = x[0],x[1],x[2]
      elif y is None:
         if type(x) in (type(0),type(0.0)):
            self.x,self.y,self.z = x,x,x
         elif type(x) is type(Vec):
            self.x,self.y,self.z = x.x,x.y,x.z
      else:
         self.x,self.y,self.z = float(x),float(y),float(z)
   def dot(u,v):
      return u.x*v.x+u.y*v.y+u.z*v.z
   def length(u):
      return math.sqrt(u.dot(u))
   def distance(u,v):
      return (u-v).length()
   def distance_squared(u,v):
      return u.x*v.x+u.y*v.y+u.z*v.z
   def d(u,v):
      return (u-v).length()
   def cross(u,v):
      return Vec(u.y*v.z-u.z*v.y,u.z*v.x-u.x*v.z,u.x*v.y-u.y*v.x)
   def __mul__(u,a):
      if type(a) is type(0) or type(a) is type(0.0):
         return Vec(u.x*a,u.y*a,u.z*a)
      elif type(a) is Vec: 
         return u.dot(a)
      else:
         # print type(a)
         assert False
   def __rmul__(u,a):
      return u*a
   def __add__(u,v):
      return Vec(u.x+v.x,u.y+v.y,u.z+v.z)
   def __radd__(u,v):
      return u+v
   def __sub__(u,v):
      return Vec(u.x-v.x,u.y-v.y,u.z-v.z)
   def __rsub__(u,v):
      return u+(-v)
   def __neg__(u):
      return Vec(-u.x,-u.y,-u.z)
   def __div__(u,a):
      return u*(1.0/a)
   def __str__(self):
      return "%f, %f, %f"%(self.x,self.y,self.z)
   def __repr__(self):
      return "Vec( %f, %f, %f )"%(self.x,self.y,self.z)
   def normalize(u):
      l = u.length()
      u.x /= l
      u.y /= l
      u.z /= l
   def normalized(u):
      v = Vec(u)
      u.normalize()
      return u
   def outer(u,v):
		return Mat( u.x*v.x, u.x*v.y, u.x*v.z,
		            u.y*v.x, u.y*v.y, u.y*v.z,
		            u.z*v.x, u.z*v.y, u.z*v.z		 )
      
      
X = Vec(1,0,0)
Y = Vec(0,1,0)
Z = Vec(0,0,1)

def randvec():
	return Vec(random.gauss(0,1),random.gauss(0,1),random.gauss(0,1))

class Mat(object):
   """docstring for Mat"""
   def __init__(self, xx=None, xy=None, xz=None, yx=None, yy=None, yz=None, zx=None, zy=None, zz=None):
      super(Mat, self).__init__()
      if xx is None: # identity default
         self.xx = 1.0
         self.yy = 1.0
         self.zz = 1.0
      if type(xx) in (type(0),type(0.0)):
         self.xx = float(xx)
         self.xy = float(xy)
         self.xz = float(xz)
         self.yx = float(yx)
         self.yy = float(yy)
         self.yz = float(yz)
         self.zx = float(zx)
         self.zy = float(zy)
         self.zz = float(zz)
      elif type(xx) is type(Vec(0)): # cols specified as vectors
         self.xx = xx.x
         self.xy = xy.x
         self.xz = xz.x
         self.yx = xx.y
         self.yy = xy.y
         self.yz = xz.y
         self.zx = xx.z
         self.zy = xy.z
         self.zz = xz.z     
      else:
         assert false
   def row(m,i):
      assert type(i) is type(1)
      if   i is 0: return Vec(m.xx,m.xy,m.xz)
      elif i is 1: return Vec(m.yx,m.yy,m.yz)
      elif i is 2: return Vec(m.zx,m.zy,m.zz)
      else: assert 0 <= i and i <= 2
   def col(m,i):
      assert type(i) is type(1)
      if   i is 0: return Vec(m.xx,m.yx,m.zx)
      elif i is 1: return Vec(m.xy,m.yy,m.zy)
      elif i is 2: return Vec(m.xz,m.yz,m.zz)
      else: assert 0 <= i and i <= 2
   def rowx(m): return m.row(0)
   def rowy(m): return m.row(1)
   def rowz(m): return m.row(2)      
   def colx(m): return m.col(0)
   def coly(m): return m.col(1)
   def colz(m): return m.col(2)      
   def __mul__(m,v):
      if type(v) in(type(0),type(0.0)):
         return Mat( v*m.xx, v*m.xy, v*m.xz, v*m.yx, v*m.yy, v*m.yz, v*m.zx, v*m.zy, v*m.zz )
      elif type(v) is Vec:
         return Vec( m.rowx()*v, m.rowy()*v, m.rowz()*v )
      elif type(v) is Mat:
         return Mat( m.rowx()*v.colx(), m.rowy()*v.colx(), m.rowz()*v.colx(),
                     m.rowx()*v.coly(), m.rowy()*v.coly(), m.rowz()*v.coly(),
                     m.rowx()*v.colz(), m.rowy()*v.colz(), m.rowz()*v.colz() )
      else:
         try:
            return v.__rmul__(m)
         except:
            print type(v)
            raise NotImplementedError
   def __rmul__(m,v):
      if type(v) in(type(0),type(0.0)):
         return Mat( v*m.xx, v*m.xy, v*m.xz, v*m.yx, v*m.yy, v*m.yz, v*m.zx, v*m.zy, v*m.zz )
      elif type(v) is Vec:
         return Vec( m.colx()*v, m.coly()*v, m.colz()*v )
      else:
         try:
            return v.__rmul__(m)
         except:
            print type(v)
            raise NotImplementedError
   def __div__(m,v):
      return m*(1/v)
   def __add__(m,v):
      if type(v) in(type(0),type(0.0)):
         return Mat( v+m.xx, v+m.xy, v+m.xz, v+m.yx, v+m.yy, v+m.yz, v+m.zx, v+m.zy, v+m.zz )
      elif type(v) is Mat:
         return Mat( v.xx+m.xx, v.xy+m.xy, v.xz+m.xz, v.yx+m.yx, v.yy+m.yy, v.yz+m.yz, v.zx+m.zx, v.zy+m.zy, v.zz+m.zz )
      else:
         try:
            return v.__rmul__(m)
         except:
            print type(v)
            raise NotImplementedError
   def __sub__(m,v):
      return m + -v
   def __neg__(m):
      return m * -1
   def __str__(m):
      return "Mat( "+str(m.rowx())+"\n     "+str(m.rowy())+"\n     "+str(m.rowz()) + "  )"
   def transpose(m):
      return Mat( m.xx, m.yx, m.zx, m.xy, m.yy, m.zy, m.xz, m.yz, m.zz )
   def trace(m):
      return m.xx+m.yy+m.zz
   def add_diagonal(m,v):
      return Mat( v.x+m.xx, m.xy, m.xz, m.yx, v.y+m.yy, m.yz, m.zx, m.zy, v.z+m.zz )

class Stub(object):
   """docstring for Stub"""
   def __init__(self, frame, cen):
      super(Stub, self).__init__()
      self.frame = frame
      self.cen = cen
   def from_four_points(s,cen,a,b,c):
      s.cen = cen
      e1 = Vec(a-b).normalized()
      e3 = Vec(e1.cross(c-b)).normalized()
      e2 = Vec(e1.cross(e3)).normalized()
      s.frame = Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z)
      return s
   def to_frame(s,x):
      return s.frame * (x - s.cen)
   def from_frame(s,x):
      return (s.frame.transpose() * x) + s.cen
      

class Jump(object):
   """docstring for Jump"""
   def __init__(self, rot, trans):
      super(Jump, self).__init__()
      assert type(rot) is Mat and type(trans) is Vec
      self.rot = rot
      self.trans = trans

def proj(u,v):
   return projection_matrix(u)*v

def projperp(u,v):
   return v - proj(u,v)



def projection_matrix(v):
   m = Mat( v.x * v.x, v.x * v.y, v.x * v.z, v.y * v.x, v.y * v.y, v.y * v.z, v.z * v.x, v.z * v.y, v.z * v.z )
   return m / v.dot(v)

def rotation_matrix_radians(axis,angle):
   n = axis.normalized()
   sin_theta = math.sin( angle )
   cos_theta = math.cos( angle )
   R = projection_matrix(n)
   R *= 1.0 - cos_theta
   R.xx += cos_theta;       R.xy -= sin_theta * n.z; R.xz += sin_theta * n.y
   R.yx += sin_theta * n.z; R.yy += cos_theta;       R.yz -= sin_theta * n.x
   R.zx -= sin_theta * n.y; R.zy += sin_theta * n.x; R.zz += cos_theta
   return R;

def rotation_matrix(axis,angle):
   return rotation_matrix_radians(axis,angle*math.pi/180.0)


def dihedral(p1,p2,p3,p4):
   a = ( p2 - p1 ).normalized()
   b = ( p3 - p2 ).normalized()
   c = ( p4 - p3 ).normalized()
   x = -a.dot(c) + a.dot(b) * b.dot(c)
   y =  a.dot( b.cross(c) );
   return abs(math.atan2(y,x)) * 180.0 / 3.14159


def angle(p1,p2,p3=None):
	if p3 is None:
   		return math.acos( p1.normalized().dot(p2.normalized()) ) * 180.0 / 3.14159
	else:
   		a = ( p2 - p1 ).normalized()
   		b = ( p2 - p3 ).normalized()
   		return math.acos( a.dot(b) ) * 180.0 / 3.14159


def sin_cos_range(x):
	assert -1.001 < x < 1.001
	return min(1.0,max(-1.0,x))

def rotation_axis(R):
	cos_theta = sin_cos_range((R.trace()-1.0)/2.0);
	tolerance = 0.000000001
	if cos_theta > -1.0+tolerance and cos_theta < 1.0-tolerance:
		x = ( 1.0 if R.zy > R.yz else -1.0 ) * math.sqrt( ( R.xx - cos_theta ) / ( 1.0 - cos_theta ) )
		y = ( 1.0 if R.xz > R.zx else -1.0 ) * math.sqrt( ( R.yy - cos_theta ) / ( 1.0 - cos_theta ) )
		z = ( 1.0 if R.yx > R.xy else -1.0 ) * math.sqrt( ( R.zz - cos_theta ) / ( 1.0 - cos_theta ) )
		theta = math.acos( cos_theta );
		assert abs( x*x + y*y + z*z - 1 ) <= 0.01
		return Vec(x,y,z),theta
	elif cos_theta >= 1.0-tolerance:
		return Vec(1.0,0.0,0.0),0.0
	else:
		nnT = R.add_diagonal(Vec(1.0,1.0,1.0)) / 2.0
		x,y,z = 0.0,0.0,0.0;
		if nnT.xx > tolerance:
			x = math.sqrt( nnT.xx )
			y = nnT.yx / x
			z = nnT.zx / x
		elif nnT.yy > tolerance:
			x = ZERO
			y = math.sqrt(nnT.yy)
			z = nnT.zy / y
		else:
			assert( nnT.zz > tolerance );
			x = ZERO
			y = ZERO
			z = sqrt( nnT.zz )
		assert abs( x*x + y*y + z*z - 1.0 ) <= 0.01
		return Vec( x, y, z ),math.pi

def test_rotation_mat():
	import random
	for i in range(10000):
		a0 = Vec(random.gauss(0.0,1.0),random.gauss(0.0,1.0),random.gauss(0.0,1.0)).normalized()
		t0 = random.uniform(0,math.pi)
		a,t = rotation_axis(rotation_matrix_radians(a0,t0))
		if t == 0.0 and t0 < 0.001:
			continue
		if abs(t-math.pi) < 0.00001:
			if (abs(a.x-a0.x) < 0.001 and abs(a.y-a0.y) < 0.001 and abs(a.z-a0.z) < 0.001) or \
			   (abs(a.x+a0.x) < 0.001 and abs(a.y+a0.y) < 0.001 and abs(a.z+a0.z) < 0.001):
				continue
			else:
				print a0
				print a
				continue
		if not abs(t-t0) < 0.0001 or not (a.normalized()-a0.normalized()).length() < 0.0001:
			print a0.normalized(), t0
			print a.normalized() , t
			print "FAIL"
			return
	print "test_rotation_mat PASS"



def test():
	test_rotation_mat()

if __name__ == '__main__':
	test()

