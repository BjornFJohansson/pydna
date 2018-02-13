#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2013, 2014 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''
This module contain functions for primer design.

'''
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA  as _IUPACAmbiguousDNA
from Bio.Seq import Seq                           as _Seq
#from Bio.SeqRecord import SeqRecord               as _SeqRecord
from pydna.seqrecord import SeqRecord             as _SeqRecord
from pydna.tm import tmbresluc                    as _tmbresluc

class Primer(_SeqRecord):
    '''This class can hold information about a primer and its position on a template 
       footprint and tail.   
    
    '''
    def __init__(self, record, 
                 *args,
                 template  = None,
                 position  = None, 
                 footprint = None,
                 concentration = 1000.0,   # nM (= 1µM)
                 **kwargs):
        if hasattr(record, "features"):
            for key, value in record.__dict__.items():
                setattr(self, key, value )
        elif hasattr(record, "alphabet"):
            super().__init__(record, *args, **kwargs)            
        else:        
            super().__init__(_Seq(record, _IUPACAmbiguousDNA()), *args, **kwargs)
        self.concentration = concentration           
        self.position      = position
        self._fp           = footprint
        self.template      = template
    
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
        return Primer(new, template = self.template, position=self.position, footprint=self._fp)
    
    def __getitem__(self, index):
        result = super().__getitem__(index)
        if hasattr(index, "indices"): # index is a slice
            i1,i2,i3 = index.indices(len(self))
            j1,j2,j3 = slice(-self._fp, None).indices(len(self))
            result._fp = self._fp - (i1-j1>0)*abs(i1-j1)
        return result
    
    def tm(self, saltc=50.0, formula=_tmbresluc):
        return formula( str(self.seq).upper(), primerc=self.concentration, saltc=saltc )

if __name__=="__main__":
    import os as _os
    cache = _os.getenv("pydna_cache", "nocache")
    _os.environ["pydna_cache"]="nocache"
    
    x=Primer("gtatcatatctatctatcta", footprint=12)
    print(str(x.seq).rjust(20))
    print(str(x.tail))
    print(str(x.footprint).rjust(20))
    
    w=x[:15]
    print(str(w.seq).ljust(20))
    print(str(w.tail))
    print(str(w.footprint).rjust(15))
    
#    y=x[-13:]
#    print("---")
#    print(str(y.seq).rjust(20))
#    print(str(y.tail))   
#    print(str(y.footprint).rjust(20))
#    
#    z=x[7:20]
#    print("===")
#    print(str(z.seq).rjust(20))
#    print(str(z.tail))   
#    print(str(z.footprint).rjust(20))
    
    
    #import doctest
    #doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
    
    
    
    
