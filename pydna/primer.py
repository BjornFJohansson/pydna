#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''This module provide the Primer class that is a subclass of the biopython SeqRecord.'''

from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA  as _IUPACAmbiguousDNA
from Bio.Seq            import Seq                as _Seq
from pydna.seqrecord    import SeqRecord          as _SeqRecord
from pydna.tm                  import tm_default   as _tm_default
from pydna.tm                  import tm_dbd       as _tm_dbd

class Primer(_SeqRecord):
    '''This class can hold information about a primer and its position on a template 
    footprint and tail.'''

    def __init__(self, record, 
                 *args,
                 amplicon  = None,
                 position  = None, 
                 footprint = 0,
                 tm_func =_tm_default,
                 tm_func_dbd =_tm_dbd,
                 **kwargs):
        
        if hasattr(record, "features"):
            for key, value in record.__dict__.items():
                setattr(self, key, value )
        elif hasattr(record, "alphabet"):
            super().__init__(record, *args, **kwargs)            
        else:        
            super().__init__(_Seq(record, _IUPACAmbiguousDNA()), *args, **kwargs)
          
        self.position      = position
        self._fp           = footprint or len(record)
        self.tm_func       = tm_func
        self.tm_func_dbd   = tm_func_dbd 
        
    @property
    def footprint(self):
        return self.seq[-self._fp:] if self._fp else ""


    @property
    def tail(self):
        return self.seq[:-self._fp] if self._fp else ""
    
    
    def __repr__(self):
        s = min( (self.seq,"{}..{}".format(self.seq[:15], self.seq[-3:])), key=len)
        return "{id} {len}-mer:5'-{seq}-3'".format(id=self.id,len=len(self),seq=s)
    

    def __radd__(self, other):
        new = super().__radd__(other)
        return Primer(new, amplicon = self.amplicon, position=self.position, footprint=self._fp)
    

    def __getitem__(self, index):
        result = super().__getitem__(index)
        if hasattr(index, "indices"): # index is a slice
            i1,i2,i3 = index.indices(len(self))
            j1,j2,j3 = slice(-(self._fp or 0), None).indices(len(self))
            result._fp = self._fp - (i1-j1>0)*abs(i1-j1)
        return result


    def tm(self):
        return self.tm_func(self.footprint)


    def tm_dbd(self):
        return self.tm_func_dbd(self.footprint)


if __name__=="__main__":
    import os as _os
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"]=""
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"]=cached