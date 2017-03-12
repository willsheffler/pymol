import doctest
import unittest

import minimock
import xyzMath
import pymol_util
import sym_util
import sym_comp
from misc import symgen
from misc import truss


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(xyzMath))
    tests.addTests(doctest.DocTestSuite(pymol_util))
    tests.addTests(doctest.DocTestSuite(sym_util))
    tests.addTests(doctest.DocTestSuite(sym_comp))
    tests.addTests(doctest.DocTestSuite(truss))
    return tests


if __name__ == '__main__':
    import doctest
    print "testing xyzMath......", doctest.testmod(xyzMath)
    print "testing pymol_util...", doctest.testmod(pymol_util)
    print "testing sym_util.....", doctest.testmod(sym_util)
    print "testing sym_comp.....", doctest.testmod(sym_comp)
    print "testing truss.py......", doctest.testmod(symgen)
    print "testing truss.py......", doctest.testmod(truss)
    print "testing test.py......", doctest.testmod()
