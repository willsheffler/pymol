some docs

>>> print 1
1

some more docs

.. testsetup:: helixparam
    import xyzMath

.. testcode:: helixparams

   1+1        # this will give no output!
   print 2+2  # this will give output

.. testoutput:: helixparams

   4

.. doctest::

    >>> import math
    >>> print math.sqrt(2.)
    1.41421356237
