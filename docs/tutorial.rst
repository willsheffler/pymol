1Tutorial
============

.. testsetup::

    from xyzMath import *

.. toctree::

    tutorial.rst


Compatibility with pyrosetta Vectors/Matrices
-------------------------------------------------------

use v.to_rosetta() and Vec(rosetta_vec)

to/from rosetta.numeric.xyzVector_double_t
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
>>> import rosetta
>>> vros = rosetta.numeric.xyzVector_double_t(1,2,3)
>>> # convert from rosetta xyzVector to Vec
>>> v = Vec(vros)
>>> print(type(v))
<class 'xyzMath.Vec'>
>>> print(v)
(1.000000,2.000000,3.000000)
>>> # convert to rosetta xyzVector from Vec
>>> u = v.to_rosetta()
>>> print(type(u))
<class 'rosetta.numeric.xyzVector_double_t'>
>>> print(u)
      1.000000000000000       2.000000000000000       3.000000000000000


to/from rosetta.numeric.xyzMatrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
>>> import rosetta
>>> mros = rosetta.numeric.xyzMatrix_double_t(0)
>>> mros.xx(1); mros.yy(1); mros.zz(1)
>>> # convert from rosetta xyzMatrix to Mat
>>> m = Mat(mros)
>>> print(type(m))
<class 'xyzMath.Mat'>
>>> print(m)
Mat[ (1.000000,0.000000,0.000000), (0.000000,1.000000,0.000000), (0.000000,0.000000,1.000000) ]
>>> # convert to rosetta xyzMatrix from Mat
>>> n = m.to_rosetta()
>>> print(type(n))
<class 'rosetta.numeric.xyzMatrix_double_t'>
>>> print(n)
      1.000000000000000       0.000000000000000       0.000000000000000
      0.000000000000000       1.000000000000000       0.000000000000000
      0.000000000000000       0.000000000000000       1.000000000000000
<BLANKLINE>




Getting helical parameters from a transform
---------------------------------------------


Get the transform you are interested in
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

first, get the transform you're interested in. if you are in pymol and want to get a transform between selection, you can do this::

    from pymol_util import getrelframe
    xform = getrelframe( 'resi 1 and name n+ca+c', 'resi 4 and name n+ca+c' )

if you have atomic coordinates (say from pyrosetta) and want to use those to define stubs, you can do something like this

>>> some_xform = rotation_around_degrees(axs=Vec(1,0,0), ang=180.0, cen=Vec(0,1,0))
>>> some_xform.t.x = 10 # translation along axis
>>> N_1  = Vec(1,0,0)
>>> CA_1 = Vec(0,1,0)
>>> C_1  = Vec(0,0,1)
>>> N_2  = some_xform * N_1
>>> CA_2 = some_xform * CA_1
>>> C_2  = some_xform * C_1
>>> print N_2, CA_2, C_2
(11.000000,2.000000,-0.000000) (10.000000,1.000000,0.000000) (10.000000,2.000000,-1.000000)


then use the stub function do get stubs (coordinate frames) from the coords:

>>> stub1 = stub(N_1, CA_1, C_1)
>>> stub2 = stub(N_2, CA_2, C_2)

now the transform that takes stub1 to stub2 is the following (approx equal to some_xform in this example):

>>> xform = stub2 * ~stub1
>>> assert xform == some_xform

Get the helical parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

axis of rotation, rotation magnitude and center of rotation can be computed with the
:func:`rotation_axis_center <xyzMath.Xform.rotation_axis_center>` method of class :class:`.Xform`

>>> axis, ang, cen = xform.rotation_axis_center()

translation along the rotation axis can be obtained:

>>> translation_along_axis_of_rotation = axis.dot(xform.t)
>>> print translation_along_axis_of_rotation
10.0

helical radius is relative to the coordinates. can be obtained with :func:`.projperp`

*Note: rounding is only so these code examples test correctly*

>>> radius_from_N_1  = projperp(axis, cen - N_1 ).length()
>>> radius_from_CA_1 = projperp(axis, cen - CA_1).length()
>>> radius_from_C_1  = projperp(axis, cen - C_1 ).length()
>>> round(radius_from_N_1, 6)
1.0
>>> round(radius_from_CA_1, 6)
0.0
>>> round(radius_from_C_1, 6)
1.414214

sanity check: transform shouldn't change the radius

>>> radius_from_N_2  = projperp(axis, cen - N_2 ).length()
>>> radius_from_CA_2 = projperp(axis, cen - CA_2).length()
>>> radius_from_C_2  = projperp(axis, cen - C_2 ).length()
>>> assert round(radius_from_N_1 , 6) == round(radius_from_N_2 , 6)
>>> assert round(radius_from_CA_1, 6) == round(radius_from_CA_2, 6)
>>> assert round(radius_from_C_1 , 6) == round(radius_from_C_2 , 6)

