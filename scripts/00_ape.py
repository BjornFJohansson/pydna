#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This code can be put in any Python module, it does not require IPython
# itself to be running already.  It only creates the magics subclass but
# doesn't instantiate it yet.


from __future__ import print_function
from IPython.core.magic import (Magics, magics_class, line_magic,
                                cell_magic, line_cell_magic)

try:
    from pydna.editor import Editor as _Ape
except ImportError:
    pass
else:
    # The class MUST call this class decorator at creation time
    @magics_class
    class MyMagics(Magics):

        _apeloader = _Ape("tclsh /home/bjorn/.ApE/apeextractor/ApE.vfs/lib/app-AppMain/AppMain.tcl")

        @line_magic
        def ape(self, line):
            import pydna
            seq=''
            #print self.shell.user_ns[line]
            try:
                seq = self.shell.user_ns[line]
            except KeyError:
                pass
            try:
                seq = pydna.read(line)
            except ValueError:
                pass
            if seq:
                seq.description = line   # new
                MyMagics._apeloader.open(seq)  #(*args,**kwargs)
            return

        @cell_magic
        def cmagic(self, line, cell):
            "my cell magic"
            return line, cell

        @line_cell_magic
        def lcmagic(self, line, cell=None):
            "Magic that works both as %lcmagic and as %%lcmagic"
            if cell is None:
                print("Called as line magic")
                return line
            else:
                print("Called as cell magic")
                return line, cell



# In order to actually use these magics, you must register them with a
# running IPython.  This code must be placed in a file that is loaded once
# IPython is up and running:
ip = get_ipython()
# You can register the class itself without instantiating it.  IPython will
# call the default constructor on it.
ip.register_magics(MyMagics)





