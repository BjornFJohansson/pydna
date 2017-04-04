#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''This module provides three ways to access a primer list specified 
   in the pydna.ini file.

'''
import os as _os

from pydna.parsers import parse_primers as _parse_primers

list_primers = _parse_primers( _os.environ["pydna_primers"] )[::-1]

dict_primers = dict((p.id, p) for p in list_primers)

for _i, _p in enumerate(list_primers):
    globals()["p{:03d}".format(_i)] = _p
           

def append_primer_list(primers):
    pass


if __name__=="__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
