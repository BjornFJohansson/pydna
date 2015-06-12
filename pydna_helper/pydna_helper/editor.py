#!/usr/bin/env python
# -*- coding: utf-8 -*-

from settings import apecommand

from pydna.editor import Editor as _Ape

_apeloader = _Ape(apecommand)

def ape(*args,**kwargs):
    return _apeloader.open(*args,**kwargs)

if __name__=="__main__":
    from pydna import read

    sr1 = read("../../tests/pUC19.gb","gb")
    sr2 = read("../../tests/pCAPs.gb","gb")
    sr3 = read(">abc\naaac")
    #print ">>", (sr1, sr2, sr3)
    #ape(sr1)
    #ape(sr2)
    #ape(sr3)
    
    ape(sr1, sr2, sr3)
    print "Done!"
    #/opt/ActiveTcl-8.6/bin/
