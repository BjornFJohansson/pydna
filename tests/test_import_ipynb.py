#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import sys

def test_ipynb_import():

    import mymodule

    assert mymodule.foo2() == "bar2"

    try:
        import IPython
    except ImportError:
        print("*** IPython not installed Jupyter notebook import not tested ***")
        assert True
    else:
        from pydna.ipynb_importer import NotebookFinder, NotebookLoader
        import mynotebook
        assert mynotebook.foo() == "bar"

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
