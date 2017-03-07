Tutorial
============

.. testsetup::

    from xyzMath import *

.. toctree::

    tutorial.rst



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

now the transform that takes sub1 to stub2 is the following (approx equal to some_xform in this example):

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

