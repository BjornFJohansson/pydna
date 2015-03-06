#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import doctest
import os
import sys
import os.path
import pkgutil
import pydna

os.environ["PYTHONDONTWRITEBYTECODE"] = "True"
sys.dont_write_bytecode = True
os.environ["pydna_cache"]  = "nocache"

def load_tests(loader, tests, pattern):
    return load_my_tests()

def load_my_tests():
    suite = unittest.TestSuite()
    for all_test_suite in unittest.defaultTestLoader.discover("."):
        for test_suite in all_test_suite:
            suite.addTests(test_suite)
    finder = doctest.DocTestFinder()
    pkgpath = os.path.dirname(pydna.__file__)
    for _,name,_ in [n for n in pkgutil.iter_modules((pkgpath,)) if n[1]!="gel"]:
        n="pydna.{}".format(name)
        my_module_with_doctests = __import__(n, None, None, n.split("."))
        for tst in finder.find(my_module_with_doctests):
            suite.addTest(doctest.DocTestCase(tst))
        del my_module_with_doctests
    return suite

def main():
    cwd = os.getcwd()
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(os.path.join(dname,"tests"))
    print "Python version:", sys.version
    print "Operating system:", os.name, sys.platform
    print "running unit tests and doctests"
    print
    runner = unittest.TextTestRunner(verbosity = 3)
    unittest.main(testRunner=runner)
    os.chdir(cwd)

if __name__ == '__main__':
    main()
