"""
Easy 3D Linear Algebra, like xyz\* in rosetta
"""
from random    import gauss,uniform
from math      import pi,sqrt,sin,cos,acos,asin,atan2,degrees,radians,copysign
from itertools import chain,product,izip
import operator as op
import re

EPS = 0.000000001
SQRTEPS = sqrt(EPS)

def isint   (x): return type(x) is int
def isfloat (x): return type(x) is float
def isnum   (x): return isint(x) or isfloat(x)
def ispoint (x): return type(x) is Point
def ispoints(x): return type(x) is Points
def isvec   (x): return type(x) is Vec
def isvecs  (x): return type(x) is Vecs
def isvorpt (x): return isvec(x) or ispoint(x)
def isline  (x): return type(x) is Line
def isplane (x): return type(x) is Plane
def ismat   (x): return type(x) is Mat
def isxform (x): return type(x) is Xform
def islist  (x): return type(x) is list
def istuple (x): return type(x) is tuple
def isiter  (x): return hasattr(x,"__iter__")

def sametype(x,y): return type(x) is type(y)

def allints  (*X): return reduce(op.and_,(type(x) is int for x in X),True)
def allfloats(*X): return reduce(op.and_,(type(x) is float for x in X),True)
def allnums  (*X): return reduce(op.and_,(isint(x) or isfloat(x) for x in X),True)
def allpoints(*X): return reduce(op.and_,(type(x) is Point for x in X),True)
def allvecs  (*X): return reduce(op.and_,(type(x) is Vec for x in X),True)
def allvorpts(*X): return reduce(op.and_,(isvec(x) or ispoint(x) for x in X),True)
def alllines (*X): return reduce(op.and_,(type(x) is Line for x in X),True)
def allplanes(*X): return reduce(op.and_,(type(x) is Plane for x in X),True)
def allmats  (*X): return reduce(op.and_,(type(x) is Mat for x in X),True)
def allxforms(*X): return reduce(op.and_,(type(x) is Xform for x in X),True)
def alllists (*X): return reduce(op.and_,(type(x) is list for x in X),True)
def alltuples(*X): return reduce(op.and_,(type(x) is tuple for x in X),True)
def alliters (*X): return reduce(op.and_,(hasattr(x,"__iter__") for x in X),True)

def anyints  (*X): return reduce(op.or_,(type(x) is int for x in X),False)
def anyfloats(*X): return reduce(op.or_,(type(x) is float for x in X),False)
def anynums  (*X): return reduce(op.or_,(isint(x) or isfloat(x) for x in X),False)
def anypoints(*X): return reduce(op.or_,(type(x) is Point for x in X),False)
def anyvecs  (*X): return reduce(op.or_,(type(x) is Vec for x in X),False)
def anyvorpts(*X): return reduce(op.or_,(isvec(x) or ispoint(x) for x in X),False)
def anylines (*X): return reduce(op.or_,(type(x) is Line for x in X),False)
def anyplanes(*X): return reduce(op.or_,(type(x) is Plane for x in X),False)
def anymats  (*X): return reduce(op.or_,(type(x) is Mat for x in X),False)
def anyxforms(*X): return reduce(op.or_,(type(x) is Xform for x in X),False)
def anylists (*X): return reduce(op.or_,(type(x) is list for x in X),False)
def anytuples(*X): return reduce(op.or_,(type(x) is tuple for x in X),False)
def anyiters (*X): return reduce(op.or_,(hasattr(x,"__iter__") for x in X),False)

def typeerror(o,t1,t2): 
   raise TypeError("unsupported operand type(s) for "+o+": '"+type(t1).__name__+"' and '"+type(t2).__name__+"'")

def stripfloats(s):
   """
   >>> stripfloats(" 1.10 100.00 0. 1.230000 3.34534500 (0.00000) ")
   ' 1.1 100 0 1.23 3.345345 (0) '
   """
   s = re.sub(r"(\b\d+[.]\d*?)0+\b",r"\1",s)
   s = re.sub(r"(\b\d+)[.]([ ,\s\)$])",r"\1\2",s)
   return s

def sin_cos_range(x):
   assert -1.001 < x < 1.001
   return min(1.0,max(-1.0,x))

class Point(object):
   """a Point like xyzVector<Real> in rosetta

   >>> p = Point(1,2,3)
   >>> p
   P(1.000000,2.000000,3.000000)
   >>> print p
   P(1,2,3)
   >>> print ispoint(p),isvec(p)
   True False
   >>> 10+p
   Traceback (most recent call last):
   TypeError: unsupported operand type(s) for +: 'int' and 'Point'
   >>> print 10*p
   P(10,20,30)
   >>> p.key()
   (1.0, 2.0, 3.0)

   elementwise mult
   >>> p*p
   Traceback (most recent call last):
   TypeError: unsupported operand type(s) for *: 'Point' and 'Point'
   >>> assert Point(1,0,-0) == Point(1,-0,0)

   >>> round(Point(1,2,3).distance(Point(3,2,1)),6)
   2.828427
   >>> r = randpoint()
   >>> assert r.distance(r) < EPS

   >>> p.angle(p)
   Traceback (most recent call last):
   AttributeError: 'Point' object has no attribute 'angle'
   >>> p.dot(p)
   Traceback (most recent call last):
   AttributeError: 'Point' object has no attribute 'dot'
   >>> p.length()
   Traceback (most recent call last):
   AttributeError: 'Point' object has no attribute 'length'
   """
   def __init__(self,x=0.0,y=None,z=None):
      if y is None:
         if isnum(x):
            self.x,self.y,self.z = (float(x),)*3
         elif isvec(x) | ispoint(x):
            self.x,self.y,self.z = x.x,x.y,x.z
         elif isiter(x): 
            i = iter(x)
            self.x,self.y,self.z = i.next(),i.next(),i.next()
         else: raise TypeError
      elif z is not None:
         assert isnum(x) and isnum(y) and isnum(z)
         self.x,self.y,self.z = float(x),float(y),float(z)
      else: raise TypeError
      assert allfloats(self.x,self.y,self.z)
   def distance_squared(p,q): 
      if not allpoints(p,q): raise TypeError("distance between vecs / ints doesn't make sense")
      return reduce(op.add,((f-g)**2 for f,g in zip(p,q)))
   def distance(u,v): return sqrt(u.distance_squared(v))
   def __sub__(u,r):
      if allpoints(u,r): return Vec(u.x-r.x,u.y-r.y,u.z-r.z)
      return u + -r
   def __rsub__(u,l):
      return l+-u
   def __eq__(self,other): 
      return ( type(self) is type(other) and
               abs(self.x-other.x) < EPS and 
               abs(self.y-other.y) < EPS and
               abs(self.z-other.z) < EPS )
   def rounded(self,sd):
      return Vec(round(self.x,sd), round(self.y,sd), round(self.z,sd) )
   def __len__(v):
      return 3
   def abs(v):
      return Vec(abs(v.x),abs(v.y),abs(v.z))
   def __getitem__(v,i):
      if i is 0: return v.x
      if i is 1: return v.x
      if i is 2: return v.x
      raise IndexError
   def tuple(v):
      return (v.x,v.y,v.z)
   def key(v):
      return v.rounded(6).tuple()
   def __iter__(p):
      yield p.x
      yield p.y
      yield p.z
   def __mul__(p,a):
      if isnum(a): return Point(a*p.x,a*p.y,a*p.z)
      typeerror('*',p,a)
   def __rmul__(p,a):
      if isnum(a): return Point(a*p.x,a*p.y,a*p.z)
      typeerror('*',a,p)
   def __repr__(self):
      return "P(%f,%f,%f)"%(self.x,self.y,self.z)
   def __str__(self):
      return stripfloats(repr(self))

class Points(list):
   pass

class Vec(Point):
   """a 3D direction

   >>> p = Point(1,1,1)
   >>> v = Vec(1,2,3)
   >>> p*v
   Traceback (most recent call last):
   TypeError: unsupported operand type(s) for *: 'Point' and 'Vec'
   >>> v*p
   Traceback (most recent call last):
   TypeError: unsupported operand type(s) for *: 'Vec' and 'Point'
   >>> v/p
   Traceback (most recent call last):
   TypeError: unsupported operand type(s) for /: 'Vec' and 'Point'
   >>> print p+v
   P(2,3,4)
   >>> print p-v
   P(0,-1,-2)
   >>> print v-p
   Traceback (most recent call last):
   TypeError: bad operand type for unary -: 'Point'
   >>> assert p-(p+v) == -v and p-p+v == v

   >>> v.distance_squared(v)
   Traceback (most recent call last):
   TypeError: distance between vecs / ints doesn't make sense

   pairwise +/-/*/div on vecs ok
   >>> print v+v
   V(2,4,6)
   >>> print v-v
   V(0,0,0)
   >>> print v*v
   V(1,4,9)
   >>> print v/v
   V(1,1,1)

   """
   def __init__(self,*args,**kwargs):
      super(Vec,self).__init__(*args,**kwargs)
   def dot(u,v):
      assert isvec(v)
      return u.x*v.x+u.y*v.y+u.z*v.z
   def cross(u,v):
      assert isvec(v)
      return Vec(u.y*v.z-u.z*v.y,u.z*v.x-u.x*v.z,u.x*v.y-u.y*v.x)
   __and__ = dot
   __or__ = cross
   def normdot(u,v): 
      assert isvec(v)
      return min(1.0,max(-1.0,u.dot(v)/u.length()/v.length()))
   def angle(u,v):
      assert isvec(v)
      d = u.normdot(v)
      if d > 1.0-EPS: return 0.0;
      if d < EPS-1.0: return pi
      return acos(d)
   def angle_degrees(u,v): return degrees(u.angle(v))
   def lineangle(u,v):
      assert isinstance(v,Vec); 
      if u.length() < SQRTEPS or v.length < SQRTEPS: return 0.0
      ang = abs(acos( u.normdot(v) ))
      return ang if ang < pi/2.0 else pi-ang
   def linemaxangle(u,v):
      return math.pi - u.lineangle(v)
   def lineangle_degrees(u,v):    return degrees(lineangle   (u,v))
   def linemaxangle_degrees(u,v): return degrees(linemaxangle(u,v))
   def length(u):         return sqrt(u.dot(u))
   def length_squared(u): return      u.dot(u)
   def unit(v):
      if   abs(v.x) > SQRTEPS: return v/v.x
      elif abs(v.y) > SQRTEPS: return v/v.y
      elif abs(v.z) > SQRTEPS: return v/v.z
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
      assert isvec(v)
      return Mat( u.x*v.x, u.x*v.y, u.x*v.z,
                  u.y*v.x, u.y*v.y, u.y*v.z,
                  u.z*v.x, u.z*v.y, u.z*v.z      )
   def __add__(v,r):
      assert isvorpt(v)
      if isnum(r): return type(v)(v.x+r,v.y+r,v.z+r)
      elif isvorpt(r):
         if isvec(v) and isvec(r): return Vec(v.x+r.x,v.y+r.y,v.z+r.z)
         if isvec(v) or  isvec(r): return Point(v.x+r.x,v.y+r.y,v.z+r.z)
         raise TypeError
      return v.__radd__(v)
   def __radd__(v,r):
      return v + r
   def __mul__(p,a):
      assert isvorpt(p)
      if isnum(a):     return type(p)(f*a for f   in p)
      if allvecs(p,a): return type(p)(f*g for f,g in izip(p,a))
      if isvorpt(a): typeerror('*',p,a)
      else: return a.__rmul__(p)
   def __rmul__(u,a):
      if anypoints(u,a): typeerror('*',a,u)
      return u*a
   def __neg__(u):
      return Vec(-u.x,-u.y,-u.z)
   def __div__(u,a):
      if anypoints(u,a): typeerror('/',u,a)
      if isnum(a): return Vec(u.x/a  ,u.y/a,  u.z/a)
      if isvec(a): return Vec(u.x/a.x,u.y/a.y,u.z/a.z)
      return a.__rdiv__(u)
   def __rdiv__(u,a):
      return a/u
   def __repr__(self):
      return "V(%f,%f,%f)"%(self.x,self.y,self.z)
   def __str__(self):
      return stripfloats(repr(self))
   def proj(v,u):
      """
      >>> print Vec(1,1,1).proj(Vec(abs(gauss(0,10)),0,0))
      V(1,0,0)
      >>> print Vec(2,2,2).proj(Vec(abs(gauss(0,10)),0,0))
      V(2,0,0)
      >>> u,v = randvec(2)
      >>> puv = v.proj(u).normalized()
      >>> assert abs(abs(puv.dot(u.normalized()))-1.0) < EPS
      """
      return u.dot(v)/u.dot(u)*u
   def perp(v,u):
      """
      >>> u = Vec(1,0,0); v = Vec(1,1,1)
      >>> print v.perp(u)
      V(0,1,1)
      >>> u,v = randvec(2)
      >>> assert abs(v.perp(u).dot(u)) < EPS
      >>> assert abs( u.dot( v.perp(u) ) ) < EPS
      """
      return v - v.proj(u)


Ux = Vec(1,0,0)
Uy = Vec(0,1,0)
Uz = Vec(0,0,1)
V0 = Vec(0,0,0)
Px = Point(1,0,0)
Py = Point(0,1,0)
Pz = Point(0,0,1)
P0 = Point(0,0,0)

class Vecs(list):
   pass


def randpoint(n=1):
   if n is 1: return Point(gauss(0,1),gauss(0,1),gauss(0,1))
   return Points(Point(gauss(0,1),gauss(0,1),gauss(0,1)) for i in range(n))

def randvec(n=1):
   if n is 1: return Vec(gauss(0,1),gauss(0,1),gauss(0,1))
   return Vecs(Vec(gauss(0,1),gauss(0,1),gauss(0,1)) for i in range(n))

def randnorm(n=1):
   """
   >>> assert abs(randnorm().length()-1.0) < 0.0000001
   """
   if n is 1: return randvec().normalized()
   return Vecs(randvec().normalized() for i in range(n))

def coplanar(x1,x2,x3,x4):
   """
   >>> u,v,w = randpoint(3)
   >>> a,b,c = (gauss(0,10) for i in range(3))
   >>> assert     coplanar(u, v, w, u + a*(u-v) + b*(v-w) + c*(w-u) )
   >>> assert not coplanar(u, v, w, u + a*(u-v) + b*(v-w) + c*(w-u) + randvec().cross(u-v) )   
   """
   if allpoints(x1,x2,x3,x4): return abs((x3-x1).dot((x2-x1).cross(x4-x3))) < SQRTEPS
   raise NotImplementedError

def rmsd(l,m):
   """
   >>> l,m = randpoint(6),randpoint(6)
   >>> rmsd(l,l)
   0.0
   """
   assert ispoints(l)
   assert ispoints(m)   
   rmsd = 0.0
   for u,v in izip(l,m): 
      rmsd += u.distance_squared(v)
   return sqrt(rmsd)

def dihedral(p1,p2,p3,p4=None):
   """
   3 Vecs or 4 points
   >>> dihedral_degrees(Px,Py,P0,Pz)
   90.0
   >>> dihedral_degrees(Px,P0,Py,Pz)
   -90.0
   >>> dihedral_degrees(Uy,Uz,Ux)
   90.0
   >>> dihedral_degrees(Uy,Uz,Ux)
   90.0
   """
   if allpoints(p1,p2,p3,p4):
      a = ( p2 - p1 ).normalized()
      b = ( p3 - p2 ).normalized()
      c = ( p4 - p3 ).normalized()
      x = -a.dot(c) + a.dot(b) * b.dot(c)
      y =  a.dot( b.cross(c) );
      return atan2(y,x)
   if allvecs(p1,p2,p3) and p4 is None:
      return dihedral(P0,P0+p1,P0+p1+p2,P0+p1+p2+p3)
   if anypoints(p1,p2,p3,p4):
      raise NotImplementedError

def dihedral_degrees(p1,p2,p3,p4=None): 
   return degrees(dihedral(p1,p2,p3,p4))

def angle(p1,p2,p3=None):
   if allvecs(p1,p2) and p3 is None:
      return p1.angle(p2)
   elif allpoints(p1,p2,p3):
         a = ( p2 - p1 ).normalized()
         b = ( p2 - p3 ).normalized()
         return acos( a.dot(b) )


class Line(object):
   """
   from a direction and a point
   >>> print Line(Ux,P0)
   Line( P(0,0,0) + r * V(1,0,0) )

   from two points: 
   >>> print Line(P0,Px)
   Line( P(1,0,0) + r * V(1,0,0) )
   >>> print Line(P0,P0)
   Traceback (most recent call last):
      assert direction.length_squared() > SQRTEPS
   AssertionError

   >>> assert Line(Ux,P0) == Line(Ux,Px)
   >>> assert Line(Ux,P0) == Line(-Ux,Px)
   >>> assert Line(Ux,P0) != Line(Ux,Py)   
   """
   def __init__(self,direction,position):
      assert ispoint(position)
      if ispoint(direction): 
         direction = position-direction
      assert direction.length_squared() > SQRTEPS
      self.d = direction.normalized()
      self.p = position
   def __eq__(l1,l2):
      return (l1.d==l2.d or l1.d==-l2.d) and l1.d.lineangle(l1.p-l2.p) < EPS
   def __str__(l):
      return "Line( %s + r * %s )"%(str(l.p),str(l.d))
   def __repr__(l):
      return "Line(%s,%s)"%(repr(p.d),repr(p.p))
   def distance(l,r):
      """
      >>> l = Line(Uy,P0)
      >>> l.distance(P0)
      0.0
      >>> round(l.distance(Px+Uz),8)
      1.41421356
      >>> round(Line(Ux,Px).distance(Point(3,2,1)) , 8)
      2.23606798

      >>> Line(Ux,P0).distance(Line(Uy,P0))
      0.0
      >>> l1 = Line(Ux,Point(0,1,2))
      >>> l2 = Line(Ux,Point(3,2,1))
      >>> round(l1.distance(l2) , 8)
      1.41421356
      >>> l3 = Line(Uz,99.0*Px)
      >>> Line(Ux,10*Py).distance(l3)
      10.0

      # >>> X = randxform()
      # >>> round(Line(X.R*Ux,X*Point(0,1,2)).distance(Line(X.R*Ux,X*Point(3,2,1))) , 8)
      """
      if ispoint(r): return (r-l.p).perp(l.d).length()
      if isvec(r): raise TypeError("Line distance to Vec not defined")
      if isline(r):
         a1 = l.d.normalized()
         a2 = r.d.normalized()   
         if abs(a1.dot(a2)) > 0.9999: return (r.p-l.p).perp(a1).length()
         a = a1
         b = a2
         c = r.p-l.p
         n = abs(c.dot(a.cross(b)))
         d = a.cross(b).length()
         if abs(d) < EPS: return 0
         return n/d

# def line_line_distance(a1,c1,a2,c2):
#    """
#    >>> line_line_distance(Ux,V0,Uy,V0)
#    0.0
#    >>> round(line_line_distance(Ux,Vec(0,1,2),Ux,Vec(3,2,1)) , 8)
#    1.41421356
#    >>> line_line_distance(Ux,10*Uy,Uz,99.0*Ux)
#    10.0

#    # >>> X = randxform()
#    # >>> round(line_line_distance(X.R*Ux,X*Vec(0,1,2),X.R*Ux,X*Vec(3,2,1)) , 8)
#    # 1.41421356
#    """
#    a1 = a1.normalized()
#    a2 = a2.normalized()   
#    if abs(a1.dot(a2)) > 0.9999: return (c1-c2).perp(a1).length()
#    a = a1
#    b = a2
#    c = c2-c1
#    n = abs(c.dot(a.cross(b)))
#    d = a.cross(b).length()
#    if abs(d) < EPS: return 0
#    return n/d


def line_plane_intersection(l,l0,n,p0):
   """
   >>> l  = Ux
   >>> l0 = randvec()
   >>> n  = Ux
   >>> p0 = V0
   >>> assert line_plane_intersection(l,l0,n,p0)[1] == Vec(0,l0.y,l0.z)
   >>> n = randnorm()
   >>> p0 = randvec().cross(n)
   >>> l = randvec()
   >>> l0 = p0+l*gauss(0,10)
   >>> assert line_plane_intersection(l,l0,n,p0)[1] == p0
   """
   n = n.normalized()
   d = (p0-l0).dot(n) / l.dot(n)
   return d,d*l+l0

def slide_to_make_lines_intersect(dof,l,l0,m,m0):
   """
   >>> v = randvec()
   >>> assert abs(slide_to_make_lines_intersect(Ux,Uy,v,Uz,V0) + v.x ) < EPS
   >>> dof,l,l0,m,m0 = randvec(5)
   >>> d = slide_to_make_lines_intersect(dof,l,l0,m,m0)
   >>> l0 = l0 + d*dof
   >>> assert abs(Line(l,P0+l0).distance(Line(m,P0+m0))) < EPS
   """
   n  = l.cross(m)
   p0 = m0
   d,i = line_plane_intersection(dof,l0,n,p0)   
   assert ( (i-l0).normalized().dot(dof.normalized()) - 1.0 ) < EPS
   assert i-l0 == dof*d
   return d


# def slide_to_make_lines_intersect(dof,l,l0,m,m0):
#    """
#    >>> v = randvec()
#    >>> assert abs(slide_to_make_lines_intersect(Ux,Uy,v,Uz,V0) + v.x ) < EPS
#    >>> dof,l,l0,m,m0 = randvec(5)
#    >>> d = slide_to_make_lines_intersect(dof,l,l0,m,m0)
#    >>> l0 = l0 + d*dof
#    >>> assert abs(line_line_distance(l,l0,m,m0)) < EPS
#    """
#    l0 = Point(l0)
#    m0 = Point(m0)
#    n  = l.cross(m)
#    d,i = Plane(n,m0).intersection(Line(dof,l0))
#    assert ( (i-l0).normalized().dot(dof.normalized()) - 1.0 ) < EPS
#    assert i == dof*d+l0
#    return d

class Plane(object):
   """
   from normal and center
   >>> print Plane(Ux,P0)
   Plane(norm=V(1,0,0),p0=P(0,0,0))

   from 3 points:
   >>> print Plane(P0,Py,Pz)
   Plane(norm=V(1,0,0),p0=P(0,0,0))

   from line and point
   >>> print Plane( Line(Uy,Pz), P0)
   Plane(norm=V(1,0,0),p0=P(0,0,0))

   from line and vec
   >>> print Plane( Line(Uy,P0), Uz)
   Plane(norm=V(1,0,0),p0=P(0,0,0))

   >>> assert Plane( Line(Uy,P0), Uz) == Plane(-Ux,P0)
   >>> assert Plane( Line(Uy,P0), Uz) != Plane(-Ux,P0+Vec(0.0001) )
   """
   def __init__(self,a,b=None,c=None):
      if isplane(a):  a,b = a.n,a.p
      elif isline(a) and ispoint(b) and c is None: a,b = a.d.cross(a.p-b),b
      elif isline(a) and isvec(b)   and c is None: a,b = a.d.cross(b),a.p
      if allpoints(a,b,c): a,b = (a-b).cross(a-c),a
      assert isvec(a) and ispoint(b)
      assert a.length_squared() > SQRTEPS
      self.n = a.normalized()
      self.p = b
   def __eq__(p1,p2):
      return (p1.n==p2.n or p1.n==-p2.n) and abs(p1.n.dot(p1.p-p2.p)) < EPS
   def __str__(p):
      return "Plane(norm=%s,p0=%s)"%(str(p.n),str(p.p))
   def __repr__(p):
      return "Plane(%s,%s)"%(repr(p.n),repr(p.p))
   def intersection(p,l):
      """
      >>> l  = Ux
      >>> l0 = randpoint()
      >>> n  = Ux
      >>> p0 = P0
      >>> assert Plane(n,p0).intersection(Line(l,l0))[1] == Point(0,l0.y,l0.z)
      >>> n = randnorm()
      >>> p0 = P0 + randvec().cross(n)
      >>> l = randvec()
      >>> l0 = p0+l*gauss(0,10)
      >>> assert Plane(n,p0).intersection(Line(l,l0))[1] == p0
      """
      n = p.n.normalized()
      d = (p.p-l.p).dot(n) / l.d.dot(n)
      return d, d*l.d + l.p




# class Mat(object):
#    """docstring for Mat

#    >>> m = Mat(2,0,0,0,1,0,0,0,1)
#    >>> print m
#    Mat[ (2.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ]
#    >>> print m*m
#    Mat[ (4.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ]
#    >>> print Mat(*range(1,10)) * Mat(*range(10,19))
#    Mat[ (84.000000,90.000000,96.000000), (201.000000,216.000000,231.000000), (318.000000,342.000000,366.000000) ]
#    >>> assert Mat(0.0,1.0,2.0,3,4,5,6,7,8) == Mat(-0,1,2,3,4,5.0,6.0,7.0,8.0)
#    >>> print Mat(100,2,3,4,5,6,7,8,9).det()
#    -297.0
#    >>> m = Mat(100,2,3,4,5,6,7,8,9)
#    >>> assert m * ~m == Imat
#    """
#    def __init__(self, xx=None, xy=None, xz=None, yx=None, yy=None, yz=None, zx=None, zy=None, zz=None):
#       super(Mat, self).__init__()
#       if xx is None: # identity default
#          self.xx, self.xy, self.xz = 1.0,0.0,0.0
#          self.yx, self.yy, self.yz = 0.0,1.0,0.0
#          self.zx, self.zy, self.zz = 0.0,0.0,1.0
#       elif xy is None and ismat(xx):
#          self.xx, self.xy, self.xz = xx.xx, xx.xy, xx.xz
#          self.yx, self.yy, self.yz = xx.yx, xx.yy, xx.yz
#          self.zx, self.zy, self.zz = xx.zx, xx.zy, xx.zz
#       elif yx is None and isvec(xx) and isvec(xy) and isvec(xz):
#          self.xx, self.xy, self.xz = xx.x, xy.x, xz.x
#          self.yx, self.yy, self.yz = xx.y, xy.y, xz.y
#          self.zx, self.zy, self.zz = xx.z, xy.z, xz.z
#       elif isnum(xx):
#          self.xx, self.xy, self.xz = float(xx), float(xy), float(xz)
#          self.yx, self.yy, self.yz = float(yx), float(yy), float(yz)
#          self.zx, self.zy, self.zz = float(zx), float(zy), float(zz)
#       else:
#          assert not isnum(xx)
#          assert not ismat(xx)
#          assert not isvec(xx)
#          raise TypeError
#       assert isfloat(self.xx) and isfloat(self.xy) and isfloat(self.xz)
#       assert isfloat(self.yx) and isfloat(self.yy) and isfloat(self.yz)
#       assert isfloat(self.zx) and isfloat(self.zy) and isfloat(self.zz)
#    def row(m,i):
#       assert isint(i)
#       if   i is 0: return Vec(m.xx,m.xy,m.xz)
#       elif i is 1: return Vec(m.yx,m.yy,m.yz)
#       elif i is 2: return Vec(m.zx,m.zy,m.zz)
#       else: assert 0 <= i and i <= 2
#    def col(m,i):
#       assert isint(i)
#       if   i is 0: return Vec(m.xx,m.yx,m.zx)
#       elif i is 1: return Vec(m.xy,m.yy,m.zy)
#       elif i is 2: return Vec(m.xz,m.yz,m.zz)
#       else: assert 0 <= i and i <= 2
#    def rowx(m): return m.row(0)
#    def rowy(m): return m.row(1)
#    def rowz(m): return m.row(2)      
#    def colx(m): return m.col(0)
#    def coly(m): return m.col(1)
#    def colz(m): return m.col(2)
#    def __invert__(m): return Mat(   m.zz*m.yy-m.zy*m.yz  , -(m.zz*m.xy-m.zy*m.xz) ,   m.yz*m.xy-m.yy*m.xz  ,
#                   -(m.zz*m.yx-m.zx*m.yz) ,   m.zz*m.xx-m.zx*m.xz  , -(m.yz*m.xx-m.yx*m.xz) ,
#                     m.zy*m.yx-m.zx*m.yy  , -(m.zy*m.xx-m.zx*m.xy) ,   m.yy*m.xx-m.yx*m.xy  ) / m.det()
#    def __mul__(m,r):
#       if   isnum(r): return Mat( r*m.xx, r*m.xy, r*m.xz, r*m.yx, r*m.yy, r*m.yz, r*m.zx, r*m.zy, r*m.zz )
#       elif isvec(r): return Vec( m.rowx()*r, m.rowy()*r, m.rowz()*r )
#       elif ismat(r): return Mat( m.rowx()&r.colx(), m.rowx()&r.coly(), m.rowx()&r.colz(),
#                                  m.rowy()&r.colx(), m.rowy()&r.coly(), m.rowy()&r.colz(),
#                                  m.rowz()&r.colx(), m.rowz()&r.coly(), m.rowz()&r.colz() )
#       else: return r.__rmul__(m)
#    def __rmul__(m,v):
#       if   isnum(v): return m*v
#       elif isvec(v): return Vec( m.colx()*v, m.coly()*v, m.colz()*v )
#    def __div__(m,v): return m*(1/v)
#    def __add__(m,v):
#       if   isnum(v): return Mat(v   +m.xx,v   +m.xy,v   +m.xz,v   +m.yx,v   +m.yy,v   +m.yz,v   +m.zx,v   +m.zy,v   +m.zz)
#       elif ismat(v): return Mat(v.xx+m.xx,v.xy+m.xy,v.xz+m.xz,v.yx+m.yx,v.yy+m.yy,v.yz+m.yz,v.zx+m.zx,v.zy+m.zy,v.zz+m.zz)
#       else: return v.__radd__(m)
#    def __sub__(m,v): return m + -v
#    def __neg__(m): return m * -1
#    def __str__(m): return "Mat[ %s, %s, %s ]" % (str(m.rowx()),str(m.rowy()),str(m.rowz()))
#    def transpose(m):
#       m = Mat( m.xx, m.yx, m.zx, m.xy, m.yy, m.zy, m.xz, m.yz, m.zz )
#    def transposed(m): return Mat( m.xx, m.yx, m.zx, m.xy, m.yy, m.zy, m.xz, m.yz, m.zz )
#    def det(m):
#           #   a11  (a33  a22- a32  a23)- a21 ( a33  a12- a32  a13)+ a31(  a23  a12- a22  a13)
#       return m.xx*(m.zz*m.yy-m.zy*m.yz)-m.yx*(m.zz*m.xy-m.zy*m.xz)+m.zx*(m.yz*m.xy-m.yy*m.xz)
#    def trace(m): 
#       return m.xx+m.yy+m.zz
#    def add_diagonal(m,v): 
#       return Mat( v.x+m.xx, m.xy, m.xz, m.yx, v.y+m.yy, m.yz, m.zx, m.zy, v.z+m.zz )
#    def is_rotation(m): 
#       return (m.colx().isnormal() and m.coly().isnormal() and m.colz().isnormal() and
#               m.rowx().isnormal() and m.rowy().isnormal() and m.rowz().isnormal()   )
#    def __eq__(self,other): return ( abs(self.xx-other.xx) < EPS and 
#                abs(self.xy-other.xy) < EPS and
#                abs(self.xz-other.xz) < EPS and
#                abs(self.yx-other.yx) < EPS and
#                abs(self.yy-other.yy) < EPS and
#                abs(self.yz-other.yz) < EPS and
#                abs(self.zx-other.zx) < EPS and
#                abs(self.zy-other.zy) < EPS and
#                abs(self.zz-other.zz) < EPS )
#    def rotation_axis(R):
#       """
#       >>> axis ,ang  = randnorm(),uniform(-pi,pi)
#       >>> axis2,ang2 = rotation_matrix(axis,ang).rotation_axis()
#       >>> assert abs( abs(ang) - abs(ang2) ) < EPS
#       >>> assert axis == axis2 * copysign(1,ang*ang2)
#       """
#       cos_theta = sin_cos_range((R.trace()-1.0)/2.0);
#       if cos_theta > -1.0+EPS and cos_theta < 1.0-EPS:
#          x = ( 1.0 if R.zy > R.yz else -1.0 ) * sqrt( max(0.0, ( R.xx - cos_theta ) / ( 1.0 - cos_theta ) ) )
#          y = ( 1.0 if R.xz > R.zx else -1.0 ) * sqrt( max(0.0, ( R.yy - cos_theta ) / ( 1.0 - cos_theta ) ) )
#          z = ( 1.0 if R.yx > R.xy else -1.0 ) * sqrt( max(0.0, ( R.zz - cos_theta ) / ( 1.0 - cos_theta ) ) )
#          theta = acos( cos_theta );
#          assert abs( x*x + y*y + z*z - 1 ) <= 0.01
#          return Vec(x,y,z),theta
#       elif cos_theta >= 1.0-EPS: return Vec(1.0,0.0,0.0),0.0
#       else:
#          nnT = (R+Imat)/2.0
#          x,y,z = 0.0,0.0,0.0;
#          if nnT.xx > EPS:
#             x = sqrt( nnT.xx )
#             y = nnT.yx / x
#             z = nnT.zx / x
#          elif nnT.yy > EPS:
#             x = 0
#             y = sqrt(nnT.yy)
#             z = nnT.zy / y
#          else:
#             assert( nnT.zz > EPS );
#             x = 0
#             y = 0
#             z = sqrt( nnT.zz )
#          assert abs( x*x + y*y + z*z - 1.0 ) <= 0.01
#          return Vec( x, y, z ),pi


# Imat = Mat(1,0,0,0,1,0,0,0,1)

# def projection_matrix(v):
#    m = Mat( v.x * v.x, v.x * v.y, v.x * v.z, v.y * v.x, v.y * v.y, v.y * v.z, v.z * v.x, v.z * v.y, v.z * v.z )
#    return m / v.dot(v)


# def rotation_matrix(axis,angle):
#    n = axis.normalized()
#    sin_theta = sin( angle )
#    cos_theta = cos( angle )
#    R = projection_matrix(n)
#    R *= 1.0 - cos_theta
#    R.xx += cos_theta;       R.xy -= sin_theta * n.z; R.xz += sin_theta * n.y
#    R.yx += sin_theta * n.z; R.yy += cos_theta;       R.yz -= sin_theta * n.x
#    R.zx -= sin_theta * n.y; R.zy += sin_theta * n.x; R.zz += cos_theta
#    return R;

# def rotation_matrix_degrees(axis,angle):
#    """ get a rotation matrix

#    >>> rx180 = rotation_matrix_degrees(Vec(1,0,0),180.0)
#    >>> rx90  = rotation_matrix_degrees(Vec(1,0,0),90.0)
#    >>> print rx90*rx90 == rx180
#    True
#    >>> r = rotation_matrix_degrees(Vec(1,0,0),45.0)
#    >>> print r
#    Mat[ (1.000000,0.000000,0.000000), (0.000000,0.707107,-0.707107), (0.000000,0.707107,0.707107) ]
#    >>> assert r*r == rx90
#    >>> assert r*r*r*r == rx180
#    >>> assert r*r*r*r*r*r*r*r == Imat
#    >>> assert ~r == r.transposed()

#    >>> ang = uniform(0,1)*360.0-180.0
#    >>> v = randvec()
#    >>> axs = randnorm()
#    >>> while(abs(v.dot(axs))>0.9): axs = randnorm()
#    >>> u = rotation_matrix_degrees(projperp(v,axs),ang)*v
#    >>> assert abs(u.angle_degrees(v)-abs(ang)) < SQRTEPS
#    >>> test_rotation_mat()
#    test_rotation_mat PASS
#    """
#    return rotation_matrix(axis,radians(angle))

# def test_rotation_mat():
#    import random
#    for i in range(100):
#       a0 = randnorm()
#       t0 = uniform(-pi,pi)
#       a,t = rotation_matrix(a0,t0).rotation_axis()
#       if t0 < 0.01: continue
#       if abs(t-pi) < EPS:
#          if (abs(a.x-a0.x) < 0.001 and abs(a.y-a0.y) < 0.001 and abs(a.z-a0.z) < 0.001) or \
#             (abs(a.x+a0.x) < 0.001 and abs(a.y+a0.y) < 0.001 and abs(a.z+a0.z) < 0.001):
#             continue
#          else:
#             print a0
#             print a
#             return False
#       if not abs(t-t0) < EPS or not (a.normalized()-a0.normalized()).length() < EPS:
#          print a0.normalized(), t0
#          print a.normalized() , t
#          print "FAIL"
#          return
#    print "test_rotation_mat PASS"


# def randrot(n=1):
#    if n is 1: return rotation_matrix_degrees(randvec(),uniform(0,1)*360)
#    return (rotation_matrix_degrees(randvec(),uniform(0,1)*360) for i in range(n))

# class Xform(object):
#    """Coordinate frame like rosetta Xform, behaves also as a rosetta Stub

#    >>> x = Xform(R=Imat,t=Uz)
#    >>> print x
#    Xform( Mat[ (1.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ], (0.000000,0.000000,1.000000) )
#    >>> assert (x*x) == Xform(R=Imat,t=2*Uz)
#    >>> x = Xform(R=rotation_matrix_degrees(Vec(1,0,0),90.0),t=Vec(0,0,0))
#    >>> print x
#    Xform( Mat[ (1.000000,0.000000,0.000000), (0.000000,0.000000,-1.000000), (0.000000,1.000000,0.000000) ], (0.000000,0.000000,0.000000) )
#    >>> assert x*x*x*x == Ixform
#    >>> x.t = Ux
#    >>> assert x*x*x*x == Xform(R=Imat,t=4*Ux)
#    >>> x.t = Uz
#    >>> print x
#    Xform( Mat[ (1.000000,0.000000,0.000000), (0.000000,0.000000,-1.000000), (0.000000,1.000000,0.000000) ], (0.000000,0.000000,1.000000) )
#    >>> assert x               == Xform(R=rotation_matrix_degrees(Ux, 90.0),t=Vec(0, 0,1))
#    >>> assert x*x             == Xform(R=rotation_matrix_degrees(Ux,180.0),t=Vec(0,-1,1))
#    >>> assert x*x*x           == Xform(R=rotation_matrix_degrees(Ux,270.0),t=Vec(0,-1,0))
#    >>> assert x*x*x*x         == Xform(R=rotation_matrix_degrees(Ux,  0.0),t=Vec(0, 0,0))
#    >>> assert x*x*x*x*x       == Xform(R=rotation_matrix_degrees(Ux, 90.0),t=Vec(0, 0,1))
#    >>> assert x*x*x*x*x*x     == Xform(R=rotation_matrix_degrees(Ux,180.0),t=Vec(0,-1,1))
#    >>> assert x*x*x*x*x*x*x   == Xform(R=rotation_matrix_degrees(Ux,270.0),t=Vec(0,-1,0))
#    >>> assert x*x*x*x*x*x*x*x == Xform(R=rotation_matrix_degrees(Ux,  0.0),t=Vec(0, 0,0))
#    >>> x = Xform(rotation_matrix_degrees(Vec(1,2,3),123),Vec(5,7,9))
#    >>> assert ~x *  x == Ixform
#    >>> assert  x * ~x == Ixform

#    Frames / RTs are interchangable:

#    >>> fr = Xform(rotation_matrix_degrees(Vec(1,2,3), 65.64),t=Vec(3,2,1))
#    >>> to = Xform(rotation_matrix_degrees(Vec(7,5,3),105.44),t=Vec(10,9,8))
#    >>> x = to/fr
#    >>> assert to/Ixform ==  to
#    >>> assert Ixform/fr == ~fr
#    >>> assert (to * ~fr) * fr == to
#    >>> assert x * fr == to

#    >>> a1 = randnorm()
#    >>> b1 = randnorm()
#    >>> ang = uniform(0,1)*360.0-180.0
#    >>> a2 = rotation_matrix_degrees(a1.cross(randnorm()),ang) * a1
#    >>> b2 = rotation_matrix_degrees(b1.cross(randnorm()),ang) * b1
#    >>> assert abs(angle(a1,a2) - angle(b1,b2)) < EPS
#    >>> xa = Xform().from_two_vecs(a1,a2)
#    >>> xb = Xform().from_two_vecs(b1,b2)
#    >>> assert xa.tolocal(a1) == xb.tolocal(b1)
#    >>> assert xa.tolocal(a2) == xb.tolocal(b2)
#    >>> assert ~xa*a1 == ~xb*b1
#    >>> assert ~xa*a2 == ~xb*b2
#    >>> assert xb/xa*a1 == b1
#    >>> assert xb/xa*a2 == b2

#    add/sub with Vecs:
#    >>> X = randxform()
#    >>> u,v = randvec(2)
#    >>> assert isxform(u+X) and isxform(X+u) and isxform(u-X) and isxform(X-u)
#    >>> assert X*v+u == (u+X)*v
#    >>> assert X*(v+u) == (X+u)*v
#    >>> assert Xform(u)*X*v == (u+X)*v
#    >>> assert X*Xform(u)*v == (X+u)*v
#    >>> assert X*v-u == (u-X)*v
#    >>> assert X*(v-u) == (X-u)*v
   
#    mul,div with Mats:
#    >>> R = randrot()
#    >>> assert isxform(R*X) and isxform(X*R)
#    >>> assert R*X*u == (R*X)*u == R*(X*u)
#    >>> assert X*R*u == (X*R)*u == X*(R*u)
#    >>> assert Xform(R)*X*u == Xform(R)*(X*u)
#    >>> assert X*Xform(R)*u == X*(Xform(R,V0)*u)
#    >>> assert X/X*v == v

#    mul/div Xforms:
#    >>> Y = randxform()
#    >>> assert isxform(X/Y) and isxform(X*Y)
#    >>> assert X/Y*v == X*~Y*v

#    # >>> axis,ang,cen = randnorm(),uniform(-pi,pi),randvec()
#    # >>> X = rotation_around(axis,ang,cen)
#    # >>> axis2,ang2,cen2 = X.rotation_center()
#    # >>> assert abs( abs(ang) - abs(ang2) ) < EPS
#    # >>> assert axis == axis2 * copysign(1,ang*ang2)
#    # >>> print cen
#    # >>> print cen2
#    """
#    def __init__(self, R=None, t=None):
#       super(Xform, self).__init__()
#       if isvec(R) and t is None: R,t = Imat,R
#       self.R = R if R else Imat
#       self.t = t if t else V0
#       assert ismat(self.R) and isvec(self.t)
#    # def rotation_center(X):
#    #    axis,ang = X.rotation_axis()
#    #    cen = -(X.R-Imat).transposed()*X.t
#    #    return axis,ang,cen
#    def from_four_points(s,cen,a,b,c):
#       s.t = cen
#       e1 = (a-b).normalized()
#       e3 = e1.cross(c-b).normalized()
#       e2 = e1.cross(e3).normalized()
#       # print "from_four_points"
#       # print e1
#       # print e2
#       # print e3            
#       s.R = Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z)
#       return s
#    def from_two_vecs(s,a,b):
#       e1 = a.normalized()
#       e2 = projperp(a,b).normalized()
#       e3 = e1.cross(e2)
#       return Xform( Mat(e1.x,e2.x,e3.x,e1.y,e2.y,e3.y,e1.z,e2.z,e3.z),V0)
#    def tolocal(s,x): return s.R.transposed() * (x - s.t)
#    def toglobal(s,x): return (s.R * x) + s.t
#    def __invert__(self):
#       R = ~self.R
#       t = R * -self.t
#       return Xform(R,t)
#    def __mul__(X,o):
#       if   isvec(o):   return X.R * o + X.t
#       elif isxform(o): return Xform(X.R*o.R,X.R*(o.t) + X.t)
#       elif ismat(o):   return Xform(X.R*o,X.t)
#       elif islist(o):  return [X*x for x in o]
#       elif istuple(o): return tuple([X*x for x in o])
#       elif isiter(o):  return (X*x for x in o)
#       else:            return o.__rmul__(X)
#    def __rmul__(X,o):
#       if ismat(o): return Xform(o*X.R,o*X.t)
#       raise TypeError      
#    def __div__(X,o):
#       if isxform(o): return X*~o
#       return o.__rdiv__(X)
#    def __add__(X,v):
#       if isvec(v): return Xform( X.R, X.t + X.R*v )
#       return v.__radd__(X)
#    def __radd__(X,v):
#       if isvec(v): return Xform( X.R, X.t + v )
#       raise TypeError
#    def __sub__(X,v):
#       if isvec(v): return Xform( X.R, X.t - X.R*v )
#       return v.__rsub__(X)
#    def __rsub__(X,v):
#       if isvec(v): return Xform( X.R, X.t - v )
#       raise TypeError
#    def __eq__(self,other): return self.R==other.R and self.t==other.t
#    def __repr__(self): return "Xform( %s, %s )" % (str(self.R),str(self.t))
#    def __eq__(X,Y):
#       assert isxform(Y)
#       return X.R == Y.R and X.t == Y.t
#    def rotation_axis(X): return X.R.rotation_axis()
#    def pretty(self):
#       a,r = self.rotation_axis()
#       if self.t.length() > EPS: return "Xform( axis=%s, ang=%f, dir=%s, dis=%f )"%(str(a),degrees(r),str(self.t.normalized()),self.t.length())
#       else: return "Xform( axis=%s, ang=%f, dir=%s, dis=%f )"%(str(a),degrees(r),str(V0),0)         

# Ixform = Xform(Imat,V0)

# def stub(cen=None, a=None, b=None, c=None): return Xform().from_four_points(cen,a,b,c)

# def randxform(n=1):
#    if n is 1: return Xform(randrot(),randvec())
#    return (Xform(randrot(),randvec()) for i in range(n))

# def rotation_around(axs,ang,cen):
#    """
#    >>> 
#    """
#    R = rotation_matrix(axs,ang)
#    return Xform(R,R*-cen+cen)

# def rotation_around_degrees(axs,ang,cen): return rotation_around(axs,radians(ang),cen)

# def test():
#    test_rotation_mat()



# def alignvector(a,b):
#    """
#    >>> u = randvec()
#    >>> v = randvec()
#    >>> assert v.angle(alignvector(u,v)*u) < EPS
#    """
#    return rotation_around(a.normalized()+b.normalized(),pi,V0)

# def alignaroundaxis(axis,u,v):
#    """
#    >>> axis = randnorm()
#    >>> u = randvec()
#    >>> angle = uniform(-pi,pi)
#    >>> v = rotation_matrix(axis,angle)*u
#    >>> uprime = alignaroundaxis(axis,u,v)*u
#    >>> assert v.angle(uprime) < EPS
#    >>> v = randvec()
#    >>> uprime = alignaroundaxis(axis,u,v)*u
#    >>> assert coplanar(V0,axis,v,uprime)
#    """
#    return rotation_around(axis, -dihedral(u,axis,V0,v), V0 )

# def alignvectors_minangle(a1,a2,b1,b2):
#    """
#    exact alignment:
#    >>> angdeg = uniform(-180,180)
#    >>> a1 = randvec()
#    >>> b1 = randnorm()*a1.length()
#    >>> l2 = gauss(0,1)
#    >>> a2 = rotation_matrix_degrees(a1.cross(randnorm()),angdeg) * a1 * l2
#    >>> b2 = rotation_matrix_degrees(b1.cross(randnorm()),angdeg) * b1 * l2
#    >>> assert abs(angle(a1,a2) - angle(b1,b2)) < EPS
#    >>> Xa2b = alignvectors_minangle(a1,a2,b1,b2)
#    >>> assert Xa2b.t.length() < EPS
#    >>> assert (Xa2b*a1).distance(b1) < EPS
#    >>> assert (Xa2b*a2).distance(b2) < EPS

#    if angle(a1,a2) != angle(b1,2b), minimize deviation
#    >>> a1,a2,b1,b2 = randvec(4)
#    >>> Xa2b = alignvectors_minangle(a1,a2,b1,b2)
#    >>> assert coplanar(b1,b2,Xa2b*a1,Xa2b*a2)
#    >>> assert (b1.angle(a1)+b2.angle(a2)) > (b1.angle(Xa2b*a1)+b2.angle(Xa2b*a2))
#    """
#    aaxis = (a1.normalized()+a2.normalized())/2.0
#    baxis = (b1.normalized()+b2.normalized())/2.0
#    Xmiddle = alignvector(aaxis,baxis)
#    assert (baxis).angle(Xmiddle*(aaxis)) < SQRTEPS
#    Xaround = alignaroundaxis(baxis, Xmiddle*a1, b1 )# 
#    X = Xaround * Xmiddle
#    assert (b1.angle(a1)+b2.angle(a2)) > (b1.angle(X*a1)+b2.angle(X*a2))
#    return X
#    # not so good if angles don't match:
#    # xa = Xform().from_two_vecs(a2,a1)
#    # xb = Xform().from_two_vecs(b2,b1)
#    # return xb/xa

# def alignvectors(a1,a2,b1,b2): return alignvectors_minangle(a1,a2,b1,b2)

# # def alignvectors_kindamindis(a1,a2,b1,b2):
# #    """
# #    >>> ang = uniform(0,1)*360.0-180.0
# #    >>> a1 = randvec()
# #    >>> b1 = randnorm()*a1.length()
# #    >>> l2 = gauss(0,1)
# #    >>> a2 = rotation_matrix_degrees(a1.cross(randnorm()),ang) * a1 * l2
# #    >>> b2 = rotation_matrix_degrees(b1.cross(randnorm()),ang) * b1 * l2
# #    >>> assert abs(angle(a1,a2) - angle(b1,b2)) < EPS
# #    >>> Xa2b = alignvectors(a1,a2,b1,b2)
# #    >>> assert Xa2b.t.length() < EPS
# #    >>> assert (Xa2b*a1).distance(b1) < EPS
# #    >>> assert (Xa2b*a2).distance(b2) < EPS

# #    >>> a1 = randvec()
# #    >>> b1 = randvec()
# #    >>> a2 = randvec()
# #    >>> b2 = randvec()
# #    >>> Xa2b = alignvectors(a1,a2,b1,b2)
# #    >>> assert coplanar(b1,b2,Xa2b*a1,Xa2b*a2)
# #    >>> if not (b1.distance(a1)+b2.distance(a2)) > (b1.distance(Xa2b*a1)+b2.distance(Xa2b*a2)):
# #    ...   print b1
# #    ...   print b2
# #    ...   print a1
# #    ...   print a2
# #    ...   print Xa2b*a1
# #    ...   print Xa2b*a2               
# #    """
# #    Xmiddle = alignvector(a1+a2,b1+b2)
# #    assert (b1+b2).angle(Xmiddle*(a1+a2)) < SQRTEPS
# #    assert (b1+b2).angle(Xmiddle*a1+Xmiddle*a2) < SQRTEPS   
# #    Xaround = alignaroundaxis(b1+b2, Xmiddle*a1, b1 )# 
# #    return Xaround * Xmiddle
# #    # xa = Xform().from_two_vecs(a2,a1)
# #    # xb = Xform().from_two_vecs(b2,b1)
# #    return xb/xa

# def get_test_generators1():
#    x1 = rotation_around_degrees(Vec(0,0,1),180,Vec(0,0,0))
#    x2 = rotation_around_degrees(Vec(1,1,1),120,Vec(1,0,0))
#    return x1,x2

# def expand_xforms(G,N=3,c=Vec(1,3,10)):
#    """
#    >>> G = get_test_generators1()
#    >>> for x in expand_xforms(G): print x*Ux
#    (-1.000000,0.000000,0.000000)
#    (1.000000,0.000000,0.000000)
#    (1.000000,-0.000000,0.000000)
#    (-1.000000,0.000000,0.000000)
#    (1.000000,-2.000000,0.000000)
#    (1.000000,0.000000,0.000000)
#    (-1.000000,2.000000,0.000000)
#    (-1.000000,0.000000,0.000000)
#    (1.000000,-2.000000,0.000000)
#    (1.000000,-0.000000,-2.000000)
#    """
#    seenit = set()
#    for Xs in chain(G,*(product(G,repeat=n) for n in range(2,N+1))):
#       X = Xs if isinstance(Xs,Xform) else reduce(Xform.__mul__,Xs)
#       v = X*c
#       key = (round(v.x,3),round(v.y,3),round(v.z,3))
#       if key not in seenit:
#          seenit.add(key)
#          yield X

# def find_identities(G,n=6,c=Vec(1,3,10)):
#    """
#    >>> G = get_test_generators1()
#    >>> for I in find_identities(G): print I.t
#    (0.000000,0.000000,0.000000)
#    (-2.000000,2.000000,2.000000)
#    (2.000000,-2.000000,2.000000)
#    """
#    for x in expand_xforms(G,n,c):
#       if (abs(x.R.xx-1.0) < 0.0000001 and
#           abs(x.R.yy-1.0) < 0.0000001 and
#           abs(x.R.zz-1.0) < 0.0000001 ):
#          yield x

# def get_cell_bounds_orthogonal_only(G,n=6,c=Vec(1,3,10)):
#    """
#    very slow... need to speed up
#    # >>> G = get_test_generators1()
#    # >>> get_cell_bounds_orthogonal_only(G[:2],12)
#    # (4.0, 4.0, 4.0)
#    """
#    mnx,mny,mnz = 9e9,9e9,9e9
#    for i in (I.t for I in find_identities(G,n)):
#       if abs(i.x) > SQRTEPS and abs(i.y) < SQRTEPS and abs(i.z) < SQRTEPS: mnx = min(mnx,abs(i.x))
#       if abs(i.x) < SQRTEPS and abs(i.y) > SQRTEPS and abs(i.z) < SQRTEPS: mny = min(mny,abs(i.y))
#       if abs(i.x) < SQRTEPS and abs(i.y) < SQRTEPS and abs(i.z) > SQRTEPS: mnz = min(mnz,abs(i.z))
#    return round(mnx,3),round(mny,3),round(mnz,3)



























if __name__ == '__main__':
   import doctest
   for i in range(10):
      r = doctest.testmod()
      print r
      if r[0] is not 0: break
