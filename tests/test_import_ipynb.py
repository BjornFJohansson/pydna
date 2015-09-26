#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nose, sys

def test_ipynb_import():

    import mymodule

    assert mymodule.foo2() == "bar2"

    try:
        import IPython
    except ImportError:
        print "IPython not installed"
        assert True
    else:
        from pydna import ipynb_importer
        import mynotebook
        assert mynotebook.foo() == "bar"

if __name__ == '__main__':
    nose.runmodule(argv=[sys.argv[0], '--nocapture'])