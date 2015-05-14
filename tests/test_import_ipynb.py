#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nose

from pydna import ipynb_importer


def test_ipynb_import():

    import sys
    
    #sys.path.append(".")
    import os

    print "#"*500
    print  sys.path   
        
    sys.path.append(os.path.dirname(os.path.realpath(__file__)))

    print sys.path

    import mynotebook

    assert mynotebook.foo() == "bar"

if __name__ == '__main__':
    nose.runmodule()

