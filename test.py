import doctest,unittest,xyzMath,pymol_util,sym_util,sym_comp

def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(xyzMath))
    tests.addTests(doctest.DocTestSuite(pymol_util))
    tests.addTests(doctest.DocTestSuite(sym_util))
    tests.addTests(doctest.DocTestSuite(sym_comp))
    return tests

if __name__ == '__main__':
	import doctest
	print "testing xyzMath......", doctest.testmod(xyzMath)
	print "testing pymol_util...", doctest.testmod(pymol_util)
	print "testing sym_util.....", doctest.testmod(sym_util)
	print "testing sym_comp.....", doctest.testmod(sym_comp)
	print "testing test.py......", doctest.testmod()

