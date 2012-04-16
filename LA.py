"""
Easy 3D Linear Algebra, like xyz\* in rosetta
"""
import math,random

TOLERANCE = 0.000000001

class Vec(object):
   """a Vector like xyzVector<Real> in rosetta

   >>> v = Vec(1,2,3)
   >>> print v, 10*v
   (1.000000,2.000000,3.000000) (10.000000,20.000000,30.000000)

   multiplication is a dot prod at the moment
   >>> v*v
   14.0
   >>> assert Vec(1,0,-0) == Vec(1,-0,0)
   """
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
      return "(%f,%f,%f)"%(self.x,self.y,self.z)
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
                  u.z*v.x, u.z*v.y, u.z*v.z      )
   def __eq__(self,other):
      return ( abs(self.x-other.x) < TOLERANCE and 
               abs(self.y-other.y) < TOLERANCE and
               abs(self.z-other.z) < TOLERANCE )

      
X = Vec(1,0,0)
Y = Vec(0,1,0)
Z = Vec(0,0,1)

def randvec():
   return Vec(random.gauss(0,1),random.gauss(0,1),random.gauss(0,1))

class Mat(object):
   """docstring for Mat

   >>> m = Mat(2,0,0,0,1,0,0,0,1)
   >>> print m
   Mat[ (2.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ]
   >>> print m*m
   Mat[ (4.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ]
   >>> print Mat(*range(1,10)) * Mat(*range(10,19))
   Mat[ (84.000000,90.000000,96.000000), (201.000000,216.000000,231.000000), (318.000000,342.000000,366.000000) ]
   >>> assert Mat(0.0,1.0,2.0,3,4,5,6,7,8) == Mat(-0,1,2,3,4,5.0,6.0,7.0,8.0)
   >>> print Mat(100,2,3,4,5,6,7,8,9).det()
   -297.0
   >>> m = Mat(100,2,3,4,5,6,7,8,9)
   >>> assert m * ~m == Identity
   """
   def __init__(self, xx=None, xy=None, xz=None, yx=None, yy=None, yz=None, zx=None, zy=None, zz=None):
      super(Mat, self).__init__()
      if xx is None: # identity default
         self.xx = 1.0
         self.xy = 0.0
         self.xz = 0.0         
         self.yx = 0.0         
         self.yy = 1.0
         self.yz = 0.0
         self.zx = 0.0
         self.zy = 0.0
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
      elif type(xx) is type(Vec(0)) and type(xy) is type(Vec(0)) and type(xz) is type(Vec(0)):
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
         print "can't make Mat"
         assert False
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
   def __invert__(m):
      return Mat(   m.zz*m.yy-m.zy*m.yz  , -(m.zz*m.xy-m.zy*m.xz) ,   m.yz*m.xy-m.yy*m.xz  ,
                  -(m.zz*m.yx-m.zx*m.yz) ,   m.zz*m.xx-m.zx*m.xz  , -(m.yz*m.xx-m.yx*m.xz) ,
                    m.zy*m.yx-m.zx*m.yy  , -(m.zy*m.xx-m.zx*m.xy) ,   m.yy*m.xx-m.yx*m.xy  ) / m.det()
   def __mul__(m,v):
      if type(v) in(type(0),type(0.0)):
         return Mat( v*m.xx, v*m.xy, v*m.xz, v*m.yx, v*m.yy, v*m.yz, v*m.zx, v*m.zy, v*m.zz )
      elif type(v) is Vec:
         return Vec( m.rowx()*v, m.rowy()*v, m.rowz()*v )
      elif type(v) is Mat:
         return Mat( m.rowx()*v.colx(), m.rowx()*v.coly(), m.rowx()*v.colz(),
                     m.rowy()*v.colx(), m.rowy()*v.coly(), m.rowy()*v.colz(),
                     m.rowz()*v.colx(), m.rowz()*v.coly(), m.rowz()*v.colz() )
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
      return "Mat[ %s, %s, %s ]" % (str(m.rowx()),str(m.rowy()),str(m.rowz()))
   def transpose(m):
      m = Mat( m.xx, m.yx, m.zx, m.xy, m.yy, m.zy, m.xz, m.yz, m.zz )
   def transposed(m):
      return Mat( m.xx, m.yx, m.zx, m.xy, m.yy, m.zy, m.xz, m.yz, m.zz )
   def det(m):
          #   a11  (a33  a22- a32  a23)- a21 ( a33  a12- a32  a13)+ a31(  a23  a12- a22  a13)
      return m.xx*(m.zz*m.yy-m.zy*m.yz)-m.yx*(m.zz*m.xy-m.zy*m.xz)+m.zx*(m.yz*m.xy-m.yy*m.xz)
   def trace(m):
      return m.xx+m.yy+m.zz
   def add_diagonal(m,v):
      return Mat( v.x+m.xx, m.xy, m.xz, m.yx, v.y+m.yy, m.yz, m.zx, m.zy, v.z+m.zz )
   def __eq__(self,other):
      return ( abs(self.xx-other.xx) < TOLERANCE and 
               abs(self.xy-other.xy) < TOLERANCE and
               abs(self.xz-other.xz) < TOLERANCE and
               abs(self.yx-other.yx) < TOLERANCE and
               abs(self.yy-other.yy) < TOLERANCE and
               abs(self.yz-other.yz) < TOLERANCE and
               abs(self.zx-other.zx) < TOLERANCE and
               abs(self.zy-other.zy) < TOLERANCE and
               abs(self.zz-other.zz) < TOLERANCE )


Identity = Mat(1,0,0,0,1,0,0,0,1)

class RT(object):
   """Coordinate frame like rosetta RT, behaves also as a rosetta Stub

   >>> x = RT(rot=Identity,cen=Z)
   >>> print x
   RT( Mat[ (1.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ], (0.000000,0.000000,1.000000) )
   >>> assert (x*x) == RT(rot=Identity,cen=2*Z)
   >>> x = RT(rot=rotation_matrix(Vec(1,0,0),90.0),cen=Vec(0,0,0))
   >>> print x
   RT( Mat[ (1.000000,0.000000,0.000000), (0.000000,0.000000,-1.000000), (0.000000,1.000000,0.000000) ], (0.000000,0.000000,0.000000) )
   >>> assert x*x*x*x == RTI
   >>> x.cen = X
   >>> assert x*x*x*x == RT(rot=Identity,cen=4*X)
   >>> x.cen = Z
   >>> print x
   RT( Mat[ (1.000000,0.000000,0.000000), (0.000000,0.000000,-1.000000), (0.000000,1.000000,0.000000) ], (0.000000,0.000000,1.000000) )
   >>> assert x               == RT(rot=rotation_matrix(X, 90.0),cen=Vec(0, 0,1))
   >>> assert x*x             == RT(rot=rotation_matrix(X,180.0),cen=Vec(0,-1,1))
   >>> assert x*x*x           == RT(rot=rotation_matrix(X,270.0),cen=Vec(0,-1,0))
   >>> assert x*x*x*x         == RT(rot=rotation_matrix(X,  0.0),cen=Vec(0, 0,0))
   >>> assert x*x*x*x*x       == RT(rot=rotation_matrix(X, 90.0),cen=Vec(0, 0,1))
   >>> assert x*x*x*x*x*x     == RT(rot=rotation_matrix(X,180.0),cen=Vec(0,-1,1))
   >>> assert x*x*x*x*x*x*x   == RT(rot=rotation_matrix(X,270.0),cen=Vec(0,-1,0))
   >>> assert x*x*x*x*x*x*x*x == RT(rot=rotation_matrix(X,  0.0),cen=Vec(0, 0,0))
   >>> x = RT(rotation_matrix(Vec(1,2,3),123),Vec(5,7,9))
   >>> assert ~x *  x == RTI
   >>> assert  x * ~x == RTI

   Frames / RTs are interchangable:

   >>> fr = RT(rotation_matrix(Vec(1,2,3), 65.64),cen=Vec(3,2,1))
   >>> to = RT(rotation_matrix(Vec(7,5,3),105.44),cen=Vec(10,9,8))
   >>> x = RT(frm=fr,to=to)
   >>> assert RT(frm=RTI,to=to) ==  to
   >>> assert RT(frm=fr,to=RTI) == ~fr
   >>> assert (to * ~fr) * fr == to
   >>> assert x * fr == to
   """
   def __init__(self, rot=None, cen=None, a=None, b=None, c=None, frm=None, to=None):
      super(RT, self).__init__()
      if rot is not None:
         self.rot = rot
         self.cen = cen
      elif frm is not None and to is not None:
         tmp = to * ~frm
         self.rot = tmp.rot
         self.cen = tmp.cen
      elif a is not None and b is not None and c is not None:
         if cen is None: cen = a
         self.from_four_points(cen,a,b,c)
      else:
         assert False
   def from_four_points(s,cen,a,b,c):
      s.cen = cen
      e1 = Vec(a-b).normalized()
      e3 = Vec(e1.cross(c-b)).normalized()
      e2 = Vec(e1.cross(e3)).normalized()
      s.rot = Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z)
      return s
   def to_frame(s,x):
      return s.rot * (x - s.cen)
   def from_frame(s,x):
      return (s.rot.transpose() * x) + s.cen
   def __invert__(self):
      rot = ~self.rot
      cen = rot * -self.cen
      return RT(rot,cen)
   def __mul__(self,other):
      rot = self.rot*other.rot
      cen = self.rot*(other.cen) + self.cen
      return RT(rot=rot,cen=cen)
   def __eq__(self,other):
      return self.rot==other.rot and self.cen==other.cen
   def __str__(self):
      return "RT( %s, %s )" % (str(self.rot),str(self.cen))

RTI = RT(rot=Identity,cen=Vec(0,0,0))

class Jump(object):
   """docstring for Jump"""
   def __init__(self, rot, trans):
      super(Jump, self).__init__()
      assert type(rot) is Mat and type(trans) is Vec
      self.rot = rot
      self.trans = trans
def proj(u,v):
   """
   >>> u = Vec(1,0,0); v = Vec(1,1,1)
   >>> proj(u,v)
   Vec( 1.000000, 0.000000, 0.000000 )
   """
   return projection_matrix(u)*v
def projperp(u,v):
   """
   >>> u = Vec(1,0,0); v = Vec(1,1,1)
   >>> projperp(u,v)
   Vec( 0.000000, 1.000000, 1.000000 )
   """
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
   """ get a rotation matrix

   >>> rx180 = rotation_matrix(Vec(1,0,0),180.0)
   >>> rx90  = rotation_matrix(Vec(1,0,0),90.0)
   >>> print rx90*rx90 == rx180
   True
   >>> r = rotation_matrix(Vec(1,0,0),45.0)
   >>> print r
   Mat[ (1.000000,0.000000,0.000000), (0.000000,0.707107,-0.707107), (0.000000,0.707107,0.707107) ]
   >>> assert r*r == rx90
   >>> assert r*r*r*r == rx180
   >>> assert r*r*r*r*r*r*r*r == Identity
   >>> assert ~r == r.transposed()
   """
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
   """
   >>> assert rotation_axis( rotation_matrix_radians(Vec(1,2,3).normalized(),1.23) ) == (Vec(1,2,3).normalized(),1.23)
   """
   cos_theta = sin_cos_range((R.trace()-1.0)/2.0);
   if cos_theta > -1.0+TOLERANCE and cos_theta < 1.0-TOLERANCE:
      x = ( 1.0 if R.zy > R.yz else -1.0 ) * math.sqrt( ( R.xx - cos_theta ) / ( 1.0 - cos_theta ) )
      y = ( 1.0 if R.xz > R.zx else -1.0 ) * math.sqrt( ( R.yy - cos_theta ) / ( 1.0 - cos_theta ) )
      z = ( 1.0 if R.yx > R.xy else -1.0 ) * math.sqrt( ( R.zz - cos_theta ) / ( 1.0 - cos_theta ) )
      theta = math.acos( cos_theta );
      assert abs( x*x + y*y + z*z - 1 ) <= 0.01
      return Vec(x,y,z),theta
   elif cos_theta >= 1.0-TOLERANCE:
      return Vec(1.0,0.0,0.0),0.0
   else:
      nnT = R.add_diagonal(Vec(1.0,1.0,1.0)) / 2.0
      x,y,z = 0.0,0.0,0.0;
      if nnT.xx > TOLERANCE:
         x = math.sqrt( nnT.xx )
         y = nnT.yx / x
         z = nnT.zx / x
      elif nnT.yy > TOLERANCE:
         x = ZERO
         y = math.sqrt(nnT.yy)
         z = nnT.zy / y
      else:
         assert( nnT.zz > TOLERANCE );
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
   import doctest
   print doctest.testmod()

