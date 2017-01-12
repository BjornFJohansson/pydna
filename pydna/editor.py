#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''This module provides a class for opening a sequence using an editor
that accepts a file as a command line argument.

ApE - A plasmid Editor [#]_ is and excellent editor for this purpose.

References
----------

.. [#] http://biologylabs.utah.edu/jorgensen/wayned/ape/

'''

import time        as _time
import tempfile    as _tempfile
import os          as _os
import subprocess  as _subprocess
import operator    as _operator

class Editor:
    '''
    The Editor class needs to be instantiated before use.

    Parameters
    ----------

    shell_command_for_editor : str
        String containing the path to the editor
        
    tmpdir : str, optional
        String containing path to the temprary directory where sequence
        files are stored before opening.

    Examples
    --------

    >>> import pydna
    >>> #ape = pydna.Editor("tclsh8.6 /home/bjorn/.ApE/apeextractor/ApE.vfs/lib/app-AppMain/AppMain.tcl")
    >>> #ape.open("aaa") # This command opens the sequence in the ApE editor

    '''
    def __init__(self, shell_command_for_editor, tmpdir=None):
        self.path_to_editor = shell_command_for_editor
        self.tmpdir = tmpdir or _os.path.join(_tempfile.gettempdir(),"ApE")
        try:
            _os.makedirs(self. tmpdir)
        except OSError:
            pass

    def open(self, seq):
        '''Open a sequence for editing in an external (DNA) editor.

        Parameters
        ----------
        args : SeqRecord or Dseqrecord object
        
        '''
        for feature in seq.features:
            qf = feature.qualifiers
            if not "label" in qf:
                try:
                    qf["label"] = qf["note"]
                except KeyError:
                    qf["label"] = "feat{}".format(len(feature))
            if not "ApEinfo_fwdcolor" in qf:
                qf["ApEinfo_fwdcolor"]="cyan"
            if not "ApEinfo_revcolor" in qf:
                qf["ApEinfo_revcolor"]="red"
        seq.features.sort(key = _operator.attrgetter("location.start"))
        
        name = seq.name.strip(".?").replace(" ","_")
        tdir = _tempfile.mkdtemp(dir=self.tmpdir)
        tpth = _os.path.join(tdir, name+".gb")
        
        with open(tpth, "w") as f:
            f.write(seq.format("gb"))

        _subprocess.Popen("{} {}".format(self.path_to_editor, tpth),
                         shell=True,
                         stdout = _tempfile.TemporaryFile(),
                         stderr = _tempfile.TemporaryFile()).pid
        _time.sleep(0.5)

apeloader = Editor( _os.getenv("pydna_ape") )

def ape(*args,**kwargs):
    return apeloader.open(*args,**kwargs)


if __name__=="__main__":
    cache = _os.getenv("pydna_cache", "nocache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
