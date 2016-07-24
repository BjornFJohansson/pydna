#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''This module provides three ways to access a primer list specified 
   in the pydna.ini file.

'''
import os as _os

from .dsdna import parse as _parse

primer_list = _parse( _os.environ["pydna_primers"] , ds=False)[::-1]

primer_dict = dict((p.id, p) for p in primer_list)

for _i, _p in enumerate(primer_list):
    globals()["p{:03d}".format(_i)] = _p

if __name__=="__main__":
    import doctest
    doctest.testmod()
