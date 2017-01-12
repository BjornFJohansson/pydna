#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2013-1016 by Björn Johansson. All rights reserved.
# This code is part of the pydna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

from pydna.dseqrecord import Dseqrecord  as _Dseqrecord
from pydna._pretty import pretty_str as _ps

class GenbankRecord(_Dseqrecord):

    def __init__(self, record, *args, item="accession", start=None, stop=None, strand=1,**kwargs):
        super().__init__(record, *args, **kwargs)
        self.item = item
        self.start = start
        self.stop = stop
        self.strand = strand
        self._repr = item
        if self.start != None and self.stop != None:
            self._repr += " {}-{}".format(self.start, self.stop)
        link = "<a href='https://www.ncbi.nlm.nih.gov/nuccore/{item}?from={start}&to={stop}&strand={strand}' target='_blank'>{text}</a>"
        self.hyperlink = _ps(link.format(item=self.item, start=self.start or "", stop=self.stop or "", strand=self.strand, text=self._repr))

    def __repr__(self):
        '''returns a short string representation of the object'''
        return "Gbank({})({}{})".format(self._repr, {True:"-", False:"o"}[self.linear],len(self))
        
    def _repr_pretty_(self, p, cycle):
        '''returns a short string representation of the object'''
        p.text("Gbank({})({}{})".format(self._repr, {True:"-", False:"o"}[self.linear],len(self)))
            
    def _repr_html_(self):        
        return self.hyperlink

if __name__=="__main__":
    import os as _os
    cache = _os.getenv("pydna_cache", "nocache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
