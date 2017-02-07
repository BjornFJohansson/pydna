#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Importing IPython Notebooks as Modules

# It is a common problem that people want to import code from IPython Notebooks.
# This is made difficult by the fact that Notebooks are not plain Python files,
# and thus cannot be imported by the regular Python machinery.
#
# Fortunately, Python provides some fairly sophisticated [hooks]
# (http://www.python.org/dev/peps/pep-0302/) into the import machinery,
# so we can actually make IPython notebooks importable without much difficulty,
# and only using public APIs.

import io     as _io
import sys    as _sys
import types  as _types
import os     as _os


try:
    import IPython  as _IPython
except ImportError:
    pass

try:
    import nbformat  as _nbformat
except ImportError:
    try:
        from IPython import nbformat as _nbformat
    except ImportError:
        pass

def _find_notebook(fullname, path=None):
    """find a notebook, given its fully qualified name and an optional path

    This turns "foo.bar" into "foo/bar.ipynb"
    and tries turning "Foo_Bar" into "Foo Bar" if Foo_Bar
    does not exist.
    """

    name = fullname.rsplit('.', 1)[-1]
    if not path:
         path = _sys.path #['']
    for d in path:
        nb_path = _os.path.join(d, name + ".ipynb")
        if _os.path.isfile(nb_path):
            return nb_path
        # let import Notebook_Name find "Notebook Name.ipynb"
        nb_path = nb_path.replace("_", " ")
        if _os.path.isfile(nb_path):
            return nb_path


class NotebookLoader(object):
    """Module Loader for IPython Notebooks"""
    def __init__(self, path=None):
        self.shell = _IPython.core.interactiveshell.InteractiveShell.instance()
        self.path = path

    def load_module(self, fullname):
        """import a notebook as a module"""
        path = _find_notebook(fullname, self.path)

        #print ("importing IPython notebook from %s" % path)

        # load the notebook object
        with _io.open(path, 'r', encoding='utf-8') as f:
            nb = _nbformat.read(f, 4)

        #print type(nb)
        #print dir(nb)


        # nbformat_minor cells nbformat metadata

        # create the module and add it to _sys.modules
        # if name in _sys.modules:
        #    return _sys.modules[name]
        mod = _types.ModuleType(fullname)
        mod.__file__ = path
        mod.__loader__ = self
        _sys.modules[fullname] = mod

        # extra work to ensure that magics that would affect the user_ns
        # actually affect the notebook module's ns
        save_user_ns = self.shell.user_ns
        self.shell.user_ns = mod.__dict__

        try:
          for cell in nb.cells:
            #if cell.cell_type == 'code' and cell.language == 'python':
            if cell.cell_type == 'code':
                #print cell.source
                #raw_input("!!!")
                # transform the input to executable Python
                code = self.shell.input_transformer_manager.transform_cell(cell.source)
                # run the code in themodule
                exec(code, mod.__dict__)
        finally:
            self.shell.user_ns = save_user_ns
        return mod

class NotebookFinder(object):
    """Module finder that locates IPython Notebooks"""
    def __init__(self):
        self.loaders = {}


    def find_module(self, fullname, path=None):
        nb_path = _find_notebook(fullname, path)
        if not nb_path:
            return
        key = path
        if path:
            # lists aren't hashable
            key = _os.path.sep.join(path)
        if key not in self.loaders:
            self.loaders[key] = NotebookLoader(path)
        return self.loaders[key]

_sys.meta_path.append(NotebookFinder())

if __name__=="__main__":
    cache = _os.getenv("pydna_cache", "nocache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
