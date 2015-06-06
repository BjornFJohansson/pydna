#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nose

def test_ipynb_import():

    import mymodule

    assert mymodule.foo2() == "bar2"

    from pydna import ipynb_importer

    import mynotebook

    assert mynotebook.foo() == "bar"    

if __name__ == '__main__':
    nose.runmodule()

