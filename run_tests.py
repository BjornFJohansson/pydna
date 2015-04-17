#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import nose

os.environ["PYTHONDONTWRITEBYTECODE"] = "True"
sys.dont_write_bytecode = True
os.environ["pydna_cache"]  = "nocache"

def main():
    cwd = os.getcwd()
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(os.path.join(dname,"tests"))
    print "Python version:", sys.version
    print "Operating system:", os.name, sys.platform
    print "running unit tests"
    print
    nose.run( ["verbosity=3", "--nocapture"] )
    os.chdir(os.path.join(dname,"pydna"))
    print "doctests"
    nose.run( ["verbosity=3","--nocapture", "with-doctest=1"])
    os.chdir(cwd)

if __name__ == '__main__':
    main()
