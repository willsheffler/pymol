"""
Easy 3D Linear Algebra, like xyz\* in rosetta
"""
import math,random
from itertools import chain,product

EPS = 0.000000001
SQRTEPS = math.sqrt(EPS)

def isnum(x): return (type(x) in (int,float))
def isiter(x): return hasattr(x,"__iter__")
def radians(deg): return deg/180.0*math.pi
def degrees(rad): return rad/math.pi*180.0

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
   def __init__(self,x=0.0,y=None,z=None):
      if y is None:
         if isnum(x):
            self.x,self.y,self.z = (float(x),)*3
         elif type(x) is Vec:
            self.x,self.y,self.z = x.x,x.y,x.z
         elif isiter(x): 
            i = iter(x)
            self.x,self.y,self.z = i.next(),i.next(),i.next()
         else: raise NotImplementedError
      elif z is not None:
         assert isnum(x) and isnum(y) and isnum(z)
         self.x,self.y,self.z = float(x),float(y),float(z)
      else: raise NotImplementedError
      assert type(self.x) is float
      assert type(self.y) is float
      assert type(self.z) is float
   def dot(u,v):
      return u.x*v.x+u.y*v.y+u.z*v.z
   def normdot(u,v):
      return min(1.0,max(-1.0,u.dot(v)/u.length()/v.length()))
   def angle(u,v):
      d = u.normdot(v)
      if d > 1.0-EPS: return 0.0;
      if d < EPS-1.0: return math.pi
      return math.acos(d)
   def angle_degrees(u,v):
      return degrees(u.angle(v))
   def lineangle(u,v):
      if u.length() < SQRTEPS or v.length < SQRTEPS: return 0.0
      ang = abs(math.acos( u.normdot(v) ))
      return ang if ang < math.pi/2.0 else math.pi-ang
   def lineangle_degrees(u,v):
      return degrees(u.lineangle(v))
   def length(u):
      return math.sqrt(u.dot(u))
   def distance(u,v):
      return (u-v).length()
   def distance_squared(u,v):
      return u.x*v.x+u.y*v.y+u.z*v.z
   def cross(u,v):
      return Vec(u.y*v.z-u.z*v.y,u.z*v.x-u.x*v.z,u.x*v.y-u.y*v.x)
   def __mul__(u,a):
      if isnum(a):
         return Vec(u.x*a,u.y*a,u.z*a)
      elif type(a) is Vec: 
         return u.dot(a)
      else:
         a.__rmul__(u)
   def __rmul__(u,a):
      return u*a
   def __add__(u,v):
      if type(v) is Vec: return Vec(u.x+v.x,u.y+v.y,u.z+v.z)
      return v.__radd__(u)
   def __sub__(u,v):
      if type(v) is Vec: return Vec(u.x-v.x,u.y-v.y,u.z-v.z)
      return v.__rsub__(u)
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
      v.normalize()
      return v
   def outer(u,v):
      return Mat( u.x*v.x, u.x*v.y, u.x*v.z,
                  u.y*v.x, u.y*v.y, u.y*v.z,
                  u.z*v.x, u.z*v.y, u.z*v.z      )
   def __eq__(self,other):
      return ( abs(self.x-other.x) < EPS and 
               abs(self.y-other.y) < EPS and
               abs(self.z-other.z) < EPS )
   def rounded(self,sd):
      return Vec(round(self.x,sd), round(self.y,sd), round(self.z,sd) )
      
Ux = Vec(1,0,0)
Uy = Vec(0,1,0)
Uz = Vec(0,0,1)
U0 = Vec(0,0,0)

def randvec(n=1):
   if n is 1: return Vec(random.gauss(0,1),random.gauss(0,1),random.gauss(0,1))
   return [Vec(random.gauss(0,1),random.gauss(0,1),random.gauss(0,1)) for i in range(n)]

def randnorm(n=1):
   """
   >>> assert abs(randnorm().length()-1.0) < 0.0000001
   """
   if n is 1: return randvec().normalized()
   return (randvec().normalized() for i in range(n))

def coplanar(x1,x2,x3,x4):
   """
   >>> u,v,w = randvec(3)
   >>> a,b,c = (random.gauss(0,10) for i in range(3))
   >>> assert     coplanar(u, v, w, u + a*(u-v) + b*(v-w) + c*(w-u) )
   >>> assert not coplanar(u, v, w, u + a*(u-v) + b*(v-w) + c*(w-u) + randvec().cross(u-v) )   
   """
   return abs((x3-x1).dot((x2-x1).cross(x4-x3))) < SQRTEPS

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
         self.xx, self.xy, self.xz = 1.0,0.0,0.0
         self.yx, self.yy, self.yz = 0.0,1.0,0.0
         self.zx, self.zy, self.zz = 0.0,0.0,1.0
      elif xy is None and type(xx) is Mat:
         self.xx, self.xy, self.xz = xx.xx, xx.xy, xx.xz
         self.yx, self.yy, self.yz = xx.yx, xx.yy, xx.yz
         self.zx, self.zy, self.zz = xx.zx, xx.zy, xx.zz
      elif yx is None and type(xx) is Vec and type(xy) is Vec and type(xz) is Vec:
         self.xx, self.xy, self.xz = xx.x, xy.x, xz.x
         self.yx, self.yy, self.yz = xx.y, xy.y, xz.y
         self.zx, self.zy, self.zz = xx.z, xy.z, xz.z
      elif type(xx) in (int,float):
         self.xx, self.xy, self.xz = float(xx), float(xy), float(xz)
         self.yx, self.yy, self.yz = float(yx), float(yy), float(yz)
         self.zx, self.zy, self.zz = float(zx), float(zy), float(zz)
      else:
         raise NotImplementedError
      assert type(self.xx) is float and type(self.xy) is float and type(self.xz) is float
      assert type(self.yx) is float and type(self.yy) is float and type(self.yz) is float
      assert type(self.zx) is float and type(self.zy) is float and type(self.zz) is float
   def row(m,i):
      assert type(i) is int
      if   i is 0: return Vec(m.xx,m.xy,m.xz)
      elif i is 1: return Vec(m.yx,m.yy,m.yz)
      elif i is 2: return Vec(m.zx,m.zy,m.zz)
      else: assert 0 <= i and i <= 2
   def col(m,i):
      assert type(i) is int
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
      if type(v) in(int,float):
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
      if type(v) in(int,float):
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
      if type(v) in(int,float):
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
   def is_rotation(m):
      return (m.colx().isnormal() and m.coly().isnormal() and m.colz().isnormal() and
              m.rowx().isnormal() and m.rowy().isnormal() and m.rowz().isnormal()   )
   def __eq__(self,other):
      return ( abs(self.xx-other.xx) < EPS and 
               abs(self.xy-other.xy) < EPS and
               abs(self.xz-other.xz) < EPS and
               abs(self.yx-other.yx) < EPS and
               abs(self.yy-other.yy) < EPS and
               abs(self.yz-other.yz) < EPS and
               abs(self.zx-other.zx) < EPS and
               abs(self.zy-other.zy) < EPS and
               abs(self.zz-other.zz) < EPS )


Identity = Mat(1,0,0,0,1,0,0,0,1)

def projection_matrix(v):
   m = Mat( v.x * v.x, v.x * v.y, v.x * v.z, v.y * v.x, v.y * v.y, v.y * v.z, v.z * v.x, v.z * v.y, v.z * v.z )
   return m / v.dot(v)
def rotation_matrix(axis,angle):
   n = axis.normalized()
   sin_theta = math.sin( angle )
   cos_theta = math.cos( angle )
   R = projection_matrix(n)
   R *= 1.0 - cos_theta
   R.xx += cos_theta;       R.xy -= sin_theta * n.z; R.xz += sin_theta * n.y
   R.yx += sin_theta * n.z; R.yy += cos_theta;       R.yz -= sin_theta * n.x
   R.zx -= sin_theta * n.y; R.zy += sin_theta * n.x; R.zz += cos_theta
   return R;

def rotation_matrix_degrees(axis,angle):
   """ get a rotation matrix

   >>> rx180 = rotation_matrix_degrees(Vec(1,0,0),180.0)
   >>> rx90  = rotation_matrix_degrees(Vec(1,0,0),90.0)
   >>> print rx90*rx90 == rx180
   True
   >>> r = rotation_matrix_degrees(Vec(1,0,0),45.0)
   >>> print r
   Mat[ (1.000000,0.000000,0.000000), (0.000000,0.707107,-0.707107), (0.000000,0.707107,0.707107) ]
   >>> assert r*r == rx90
   >>> assert r*r*r*r == rx180
   >>> assert r*r*r*r*r*r*r*r == Identity
   >>> assert ~r == r.transposed()

   >>> ang = random.random()*360.0-180.0
   >>> v = randvec()
   >>> axs = randnorm()
   >>> while(abs(v.dot(axs))>0.9): axs = randnorm()
   >>> u = rotation_matrix_degrees(projperp(v,axs),ang)*v
   >>> assert abs(u.angle_degrees(v)-abs(ang)) < SQRTEPS
   >>> test_rotation_mat()
   test_rotation_mat PASS
   """
   return rotation_matrix(axis,radians(angle))

def test_rotation_mat():
   import random
   for i in range(100):
      a0 = randnorm()
      t0 = random.uniform(-math.pi,math.pi)
      a,t = rotation_axis(rotation_matrix(a0,t0))
      if t0 < 0.01: continue
      if abs(t-math.pi) < EPS:
         if (abs(a.x-a0.x) < 0.001 and abs(a.y-a0.y) < 0.001 and abs(a.z-a0.z) < 0.001) or \
            (abs(a.x+a0.x) < 0.001 and abs(a.y+a0.y) < 0.001 and abs(a.z+a0.z) < 0.001):
            continue
         else:
            print a0
            print a
            return False
      if not abs(t-t0) < EPS or not (a.normalized()-a0.normalized()).length() < EPS:
         print a0.normalized(), t0
         print a.normalized() , t
         print "FAIL"
         return
   print "test_rotation_mat PASS"


def randrot(n=1):
   if n is 1: return rotation_matrix_degrees(randvec(),random.random()*360)
   return (rotation_matrix_degrees(randvec(),random.random()*360) for i in range(n))

class Xform(object):
   """Coordinate frame like rosetta Xform, behaves also as a rosetta Stub

   >>> x = Xform(R=Identity,t=Uz)
   >>> print x
   Xform( Mat[ (1.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ], (0.000000,0.000000,1.000000) )
   >>> assert (x*x) == Xform(R=Identity,t=2*Uz)
   >>> x = Xform(R=rotation_matrix_degrees(Vec(1,0,0),90.0),t=Vec(0,0,0))
   >>> print x
   Xform( Mat[ (1.000000,0.000000,0.000000), (0.000000,0.000000,-1.000000), (0.000000,1.000000,0.000000) ], (0.000000,0.000000,0.000000) )
   >>> assert x*x*x*x == RTI
   >>> x.t = Ux
   >>> assert x*x*x*x == Xform(R=Identity,t=4*Ux)
   >>> x.t = Uz
   >>> print x
   Xform( Mat[ (1.000000,0.000000,0.000000), (0.000000,0.000000,-1.000000), (0.000000,1.000000,0.000000) ], (0.000000,0.000000,1.000000) )
   >>> assert x               == Xform(R=rotation_matrix_degrees(Ux, 90.0),t=Vec(0, 0,1))
   >>> assert x*x             == Xform(R=rotation_matrix_degrees(Ux,180.0),t=Vec(0,-1,1))
   >>> assert x*x*x           == Xform(R=rotation_matrix_degrees(Ux,270.0),t=Vec(0,-1,0))
   >>> assert x*x*x*x         == Xform(R=rotation_matrix_degrees(Ux,  0.0),t=Vec(0, 0,0))
   >>> assert x*x*x*x*x       == Xform(R=rotation_matrix_degrees(Ux, 90.0),t=Vec(0, 0,1))
   >>> assert x*x*x*x*x*x     == Xform(R=rotation_matrix_degrees(Ux,180.0),t=Vec(0,-1,1))
   >>> assert x*x*x*x*x*x*x   == Xform(R=rotation_matrix_degrees(Ux,270.0),t=Vec(0,-1,0))
   >>> assert x*x*x*x*x*x*x*x == Xform(R=rotation_matrix_degrees(Ux,  0.0),t=Vec(0, 0,0))
   >>> x = Xform(rotation_matrix_degrees(Vec(1,2,3),123),Vec(5,7,9))
   >>> assert ~x *  x == RTI
   >>> assert  x * ~x == RTI

   Frames / RTs are interchangable:

   >>> fr = Xform(rotation_matrix_degrees(Vec(1,2,3), 65.64),t=Vec(3,2,1))
   >>> to = Xform(rotation_matrix_degrees(Vec(7,5,3),105.44),t=Vec(10,9,8))
   >>> x = to/fr
   >>> assert to/RTI ==  to
   >>> assert RTI/fr == ~fr
   >>> assert (to * ~fr) * fr == to
   >>> assert x * fr == to

   >>> a1 = randnorm()
   >>> b1 = randnorm()
   >>> ang = random.random()*360.0-180.0
   >>> a2 = rotation_matrix_degrees(a1.cross(randnorm()),ang) * a1
   >>> b2 = rotation_matrix_degrees(b1.cross(randnorm()),ang) * b1
   >>> assert abs(angle(a1,a2) - angle(b1,b2)) < EPS
   >>> xa = Xform().from_two_vecs(a1,a2)
   >>> xb = Xform().from_two_vecs(b1,b2)
   >>> assert xa.tolocal(a1).distance(xb.tolocal(b1)) < EPS
   >>> assert xa.tolocal(a2).distance(xb.tolocal(b2)) < EPS
   >>> assert (~xa*a1).distance(~xb*b1) < EPS
   >>> assert (~xa*a2).distance(~xb*b2) < EPS   
   >>> assert b1.distance((xb/xa)*a1) < EPS
   >>> assert b2.distance((xb/xa)*a2) < EPS
   >>> x = randxform()
   >>> u,v = randvec(2)
   >>> assert ((x*v)+u).distance((u+x)*v) < EPS
   >>> assert (x*(v+u)).distance((x+u)*v) < EPS
   >>> assert (Xform(u)*x*v).distance((u+x)*v) < EPS
   >>> assert (x*Xform(u)*v).distance((x+u)*v) < EPS
   >>> assert ((x*v)-u).distance((u-x)*v) < EPS
   >>> assert (x*(v-u)).distance((x-u)*v) < EPS
   >>> R = randrot()
   >>> assert ((R*x)*u).distance(R*(x*u)) < EPS
   >>> assert ((x*R)*u).distance(x*(R*u)) < EPS
   >>> assert ((Xform(R)*x)*u).distance(Xform(R,U0)*(x*u)) < EPS
   >>> assert ((x*Xform(R))*u).distance(x*(Xform(R,U0)*u)) < EPS
   >>> assert ((x/x)*v).distance(v) < EPS
   >>> y = randxform()
   >>> assert ((x/y)*v).distance(x*~y*v) < EPS
   """
   def __init__(self, R=None, t=None):
      super(Xform, self).__init__()
      if type(R) is Vec and t is None: R,t = Identity,R
      self.R = R if R else Identity
      self.t = t if t else U0
      assert type(self.R) is Mat
      assert type(self.t) is Vec
   def from_four_points(s,cen,a,b,c):
      s.t = cen
      e1 = (a-b).normalized()
      e3 = e1.cross(c-b).normalized()
      e2 = e1.cross(e3).normalized()
      # print "from_four_points"
      # print e1
      # print e2
      # print e3            
      s.R = Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z)
      return s
   def from_two_vecs(s,a,b):
      e1 = a.normalized()
      e2 = projperp(a,b).normalized()
      e3 = e1.cross(e2)
      return Xform( Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z),U0)
   def tolocal(s,x):
      return s.R.transposed() * (x - s.t)
   def toglobal(s,x):
      return (s.R * x) + s.t
   def __invert__(self):
      R = ~self.R
      t = R * -self.t
      return Xform(R,t)
   def __mul__(X,o):
      if(type(o) is Vec):
         return X.R * o + X.t
      elif type(o) is Xform:
         return Xform(X.R*o.R,X.R*(o.t) + X.t)
      elif type(o) is Mat:
         return Xform(X.R*o,X.t)
      else:
         return o.__rmul__(X)
   def __rmul__(X,o):
      if type(o) is Mat:
         return Xform(o*X.R,o*X.t)
      raise NotImplementedError      
   def __div__(X,o):
      if type(o) is Xform: return X*~o
      return o.__rdiv__(X)
   def __add__(X,v):
      if type(v) is Vec: return Xform( X.R, X.t + X.R*v )
      return v.__radd__(X)
   def __radd__(X,v):
      if type(v) is Vec: return Xform( X.R, X.t + v )
      raise NotImplementedError
   def __sub__(X,v):
      if type(v) is Vec: return Xform( X.R, X.t - X.R*v )
      return v.__rsub__(X)
   def __rsub__(X,v):
      if type(v) is Vec: return Xform( X.R, X.t - v )
      raise NotImplementedError
   def __eq__(self,other):
      return self.R==other.R and self.t==other.t
   def __repr__(self):
      return "Xform( %s, %s )" % (str(self.R),str(self.t))
   def pretty(self):
      a,r = rotation_axis(self.R)
      if self.t.length() > EPS:
         return "Xform( axis=%s, ang=%f, dir=%s, dis=%f )"%(str(a),math.degrees(r),str(self.t.normalized()),self.t.length())
      else:
         return "Xform( axis=%s, ang=%f, dir=%s, dis=%f )"%(str(a),math.degrees(r),str(U0),0)         

RTI = Xform(R=Identity,t=Vec(0,0,0))
def stub(cen=None, a=None, b=None, c=None):
   return Xform().from_four_points(cen,a,b,c)

def randxform(n=1):
   if n is 1: return Xform(randrot(),randvec())
   return (Xform(randrot(),randvec()) for i in range(n))

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

def rotation_around(axs,ang,cen):
   R = rotation_matrix(axs,ang)
   return Xform(R,R*-cen+cen)

def rotation_around_degrees(axs,ang,cen):
   return rotation_around(axs,radians(ang),cen)

def test():
   test_rotation_mat()


def dihedral(p1,p2,p3,p4):
   """
   >>> dihedral_degrees(Ux,Uy,U0,Uz)
   90.0
   >>> dihedral_degrees(Ux,U0,Uy,Uz)
   -90.0
   """
   a = ( p2 - p1 ).normalized()
   b = ( p3 - p2 ).normalized()
   c = ( p4 - p3 ).normalized()
   x = -a.dot(c) + a.dot(b) * b.dot(c)
   y =  a.dot( b.cross(c) );
   return math.atan2(y,x)

def dihedral_degrees(p1,p2,p3,p4):
   return degrees(dihedral(p1,p2,p3,p4))

def angle(p1,p2,p3=None):
   if p3 is None:
         return math.acos( p1.normalized().dot(p2.normalized()) )
   else:
         a = ( p2 - p1 ).normalized()
         b = ( p2 - p3 ).normalized()
         return math.acos( a.dot(b) )


def sin_cos_range(x):
   assert -1.001 < x < 1.001
   return min(1.0,max(-1.0,x))

def rotation_axis(R):
   """
   >>> assert rotation_axis( rotation_matrix(Vec(1,2,3).normalized(),1.23) ) == (Vec(1,2,3).normalized(),1.23)
   """
   cos_theta = sin_cos_range((R.trace()-1.0)/2.0);
   if cos_theta > -1.0+EPS and cos_theta < 1.0-EPS:
      x = ( 1.0 if R.zy > R.yz else -1.0 ) * math.sqrt( ( R.xx - cos_theta ) / ( 1.0 - cos_theta ) )
      y = ( 1.0 if R.xz > R.zx else -1.0 ) * math.sqrt( ( R.yy - cos_theta ) / ( 1.0 - cos_theta ) )
      z = ( 1.0 if R.yx > R.xy else -1.0 ) * math.sqrt( ( R.zz - cos_theta ) / ( 1.0 - cos_theta ) )
      theta = math.acos( cos_theta );
      assert abs( x*x + y*y + z*z - 1 ) <= 0.01
      return Vec(x,y,z),theta
   elif cos_theta >= 1.0-EPS:
      return Vec(1.0,0.0,0.0),0.0
   else:
      nnT = R.add_diagonal(Vec(1.0,1.0,1.0)) / 2.0
      x,y,z = 0.0,0.0,0.0;
      if nnT.xx > EPS:
         x = math.sqrt( nnT.xx )
         y = nnT.yx / x
         z = nnT.zx / x
      elif nnT.yy > EPS:
         x = 0
         y = math.sqrt(nnT.yy)
         z = nnT.zy / y
      else:
         assert( nnT.zz > EPS );
         x = 0
         y = 0
         z = math.sqrt( nnT.zz )
      assert abs( x*x + y*y + z*z - 1.0 ) <= 0.01
      return Vec( x, y, z ),math.pi

def point_line_distance(p,a,c):
   """
   >>> point_line_distance(U0,Uy,U0)
   0.0
   >>> round(point_line_distance(U0,Uy,Ux+Uz),8)
   1.41421356
   >>> round(point_line_distance(Ux,Ux,Vec(3,2,1)) , 8)
   2.23606798
   """
   return projperp(a,p-c).length()

def line_line_distance(a1,c1,a2,c2):
   """
   >>> line_line_distance(Ux,U0,Uy,U0)
   0.0
   >>> round(line_line_distance(Ux,Vec(0,1,2),Ux,Vec(3,2,1)) , 8)
   1.41421356
   >>> line_line_distance(Ux,10*Uy,Uz,99.0*Ux)
   10.0
   >>> X = randxform()
   >>> round(line_line_distance(X.R*Ux,X*Vec(0,1,2),X.R*Ux,X*Vec(3,2,1)) , 8)
   1.41421356
   """
   a1 = a1.normalized()
   a2 = a2.normalized()   
   if abs(a1.dot(a2)) > 0.9999:
      return projperp(a1,c2-c1).length()
   a = a1
   b = a2
   c = c2-c1
   n = abs(c.dot(a.cross(b)))
   d = a.cross(b).length()
   if abs(d) < EPS: return 0
   return n/d

def alignvector(a,b):
   """
   >>> u = randvec()
   >>> v = randvec()
   >>> assert v.angle(alignvector(u,v)*u) < EPS
   """
   return rotation_around(a.normalized()+b.normalized(),math.pi,U0)

def alignaroundaxis(axis,u,v):
   """
   >>> axis = randnorm()
   >>> u = randvec()
   >>> angle = (random.random()-0.5)*2*math.pi
   >>> v = rotation_matrix(axis,angle)*u
   >>> uprime = alignaroundaxis(axis,u,v)*u
   >>> assert v.angle(uprime) < EPS
   >>> v = randvec()
   >>> uprime = alignaroundaxis(axis,u,v)*u
   >>> assert coplanar(U0,axis,v,uprime)
   """
   return rotation_around(axis, -dihedral(u,axis,U0,v), U0 )

def alignvectors_minangle(a1,a2,b1,b2):
   """
   exact alignment:
   >>> ang = random.random()*360.0-180.0
   >>> a1 = randvec()
   >>> b1 = randnorm()*a1.length()
   >>> l2 = random.gauss(0,1)
   >>> a2 = rotation_matrix_degrees(a1.cross(randnorm()),ang) * a1 * l2
   >>> b2 = rotation_matrix_degrees(b1.cross(randnorm()),ang) * b1 * l2
   >>> assert abs(angle(a1,a2) - angle(b1,b2)) < EPS
   >>> Xa2b = alignvectors_minangle(a1,a2,b1,b2)
   >>> assert Xa2b.t.length() < EPS
   >>> assert (Xa2b*a1).distance(b1) < EPS
   >>> assert (Xa2b*a2).distance(b2) < EPS

   if angle(a1,a2) != angle(b1,2b), minimize deviation
   >>> a1,a2,b1,b2 = randvec(4)
   >>> Xa2b = alignvectors_minangle(a1,a2,b1,b2)
   >>> assert coplanar(b1,b2,Xa2b*a1,Xa2b*a2)
   >>> assert (b1.angle(a1)+b2.angle(a2)) > (b1.angle(Xa2b*a1)+b2.angle(Xa2b*a2))
   """
   aaxis = (a1.normalized()+a2.normalized())/2.0
   baxis = (b1.normalized()+b2.normalized())/2.0
   Xmiddle = alignvector(aaxis,baxis)
   assert (baxis).angle(Xmiddle*(aaxis)) < SQRTEPS
   Xaround = alignaroundaxis(baxis, Xmiddle*a1, b1 )# 
   return Xaround * Xmiddle
   # not so good if angles don't match:
   # xa = Xform().from_two_vecs(a2,a1)
   # xb = Xform().from_two_vecs(b2,b1)
   # return xb/xa

def alignvectors(a1,a2,b1,b2):
   return alignvectors_minangle(a1,a2,b1,b2)

# def alignvectors_kindamindis(a1,a2,b1,b2):
#    """
#    >>> ang = random.random()*360.0-180.0
#    >>> a1 = randvec()
#    >>> b1 = randnorm()*a1.length()
#    >>> l2 = random.gauss(0,1)
#    >>> a2 = rotation_matrix_degrees(a1.cross(randnorm()),ang) * a1 * l2
#    >>> b2 = rotation_matrix_degrees(b1.cross(randnorm()),ang) * b1 * l2
#    >>> assert abs(angle(a1,a2) - angle(b1,b2)) < EPS
#    >>> Xa2b = alignvectors(a1,a2,b1,b2)
#    >>> assert Xa2b.t.length() < EPS
#    >>> assert (Xa2b*a1).distance(b1) < EPS
#    >>> assert (Xa2b*a2).distance(b2) < EPS

#    >>> a1 = randvec()
#    >>> b1 = randvec()
#    >>> a2 = randvec()
#    >>> b2 = randvec()
#    >>> Xa2b = alignvectors(a1,a2,b1,b2)
#    >>> assert coplanar(b1,b2,Xa2b*a1,Xa2b*a2)
#    >>> if not (b1.distance(a1)+b2.distance(a2)) > (b1.distance(Xa2b*a1)+b2.distance(Xa2b*a2)):
#    ...   print b1
#    ...   print b2
#    ...   print a1
#    ...   print a2
#    ...   print Xa2b*a1
#    ...   print Xa2b*a2               
#    """
#    Xmiddle = alignvector(a1+a2,b1+b2)
#    assert (b1+b2).angle(Xmiddle*(a1+a2)) < SQRTEPS
#    assert (b1+b2).angle(Xmiddle*a1+Xmiddle*a2) < SQRTEPS   
#    Xaround = alignaroundaxis(b1+b2, Xmiddle*a1, b1 )# 
#    return Xaround * Xmiddle
#    # xa = Xform().from_two_vecs(a2,a1)
#    # xb = Xform().from_two_vecs(b2,b1)
#    return xb/xa

def get_test_generators1():
   x1 = rotation_around_degrees(Vec(0,0,1),180,Vec(0,0,0))
   x2 = rotation_around_degrees(Vec(1,1,1),120,Vec(1,0,0))
   return x1,x2

def expand_xforms(G,N=3,c=Vec(1,3,10)):
   """
   >>> G = get_test_generators1()
   >>> for x in expand_xforms(G): print x*Ux
   (-1.000000,0.000000,0.000000)
   (1.000000,0.000000,0.000000)
   (1.000000,-0.000000,0.000000)
   (-1.000000,0.000000,0.000000)
   (1.000000,-2.000000,0.000000)
   (1.000000,0.000000,0.000000)
   (-1.000000,2.000000,0.000000)
   (-1.000000,0.000000,0.000000)
   (1.000000,-2.000000,0.000000)
   (1.000000,-0.000000,-2.000000)
   """
   seenit = set()
   for Xs in chain(G,*(product(G,repeat=n) for n in range(2,N+1))):
      X = Xs if isinstance(Xs,Xform) else reduce(Xform.__mul__,Xs)
      v = X*c
      key = (round(v.x,3),round(v.y,3),round(v.z,3))
      if key not in seenit:
         seenit.add(key)
         yield X

def find_identities(G,n=6,c=Vec(1,3,10)):
   """
   >>> G = get_test_generators1()
   >>> for I in find_identities(G): print I.t
   (0.000000,0.000000,0.000000)
   (-2.000000,2.000000,2.000000)
   (2.000000,-2.000000,2.000000)
   """
   for x in expand_xforms(G,n,c):
      if (abs(x.R.xx-1.0) < 0.0000001 and
          abs(x.R.yy-1.0) < 0.0000001 and
          abs(x.R.zz-1.0) < 0.0000001 ):
         yield x

def get_cell_bounds_orthogonal_only(G,n=6,c=Vec(1,3,10)):
   """
   very slow... need to speed up
   # >>> G = get_test_generators1()
   # >>> get_cell_bounds_orthogonal_only(G[:2],12)
   # (4.0, 4.0, 4.0)
   """
   mnx,mny,mnz = 9e9,9e9,9e9
   for i in (I.t for I in find_identities(G,n)):
      if abs(i.x) > SQRTEPS and abs(i.y) < SQRTEPS and abs(i.z) < SQRTEPS: mnx = min(mnx,abs(i.x))
      if abs(i.x) < SQRTEPS and abs(i.y) > SQRTEPS and abs(i.z) < SQRTEPS: mny = min(mny,abs(i.y))
      if abs(i.x) < SQRTEPS and abs(i.y) < SQRTEPS and abs(i.z) > SQRTEPS: mnz = min(mnz,abs(i.z))
   return round(mnx,3),round(mny,3),round(mnz,3)




























if __name__ == '__main__':
   import doctest
   for i in range(10):
      r = doctest.testmod()
      print r
      if r[0] is not 0: break
