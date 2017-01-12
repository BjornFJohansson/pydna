#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2013 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

from pydna.dseqrecord import Dseqrecord as _Dseqrecord

class GenbankFile(_Dseqrecord):

    def __init__(self, record, *args, path=None, **kwargs):
        super().__init__(record, *args, **kwargs)
        self.path=path

    def __repr__(self):
        '''returns a short string representation of the object'''
        return "File({})({}{})".format(self.id, {True:"-", False:"o"}[self.linear],len(self))
        
    def _repr_pretty_(self, p, cycle):
        '''returns a short string representation of the object'''
        p.text("File({})({}{})".format(self.id, {True:"-", False:"o"}[self.linear],len(self)))
            
    def _repr_html_(self):
        return "<a href='{path}' target='_blank'>{path}</a><br>".format(path=self.path)

if __name__=="__main__":
    import os as _os
    cache = _os.getenv("pydna_cache", "nocache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache     
