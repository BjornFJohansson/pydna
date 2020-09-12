#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""This module provides three ways to access a primer list specified 
in the primers entry in the pydna.ini file or in the pydna_primers
environment variable.

The list has to have sequences in FASTA, Genbank or EMBL formats.
The primer list can have the format below for example:

::
    
    >first_primer
    tgatcgtcatgctgactatactat
    >second_primer
    ctaggatcgtagatctagctg
    ...

list_primers is a list :class:`pydna.primer.Primer` objects
dict_primers is a dict where the key is the id of the object

each primer is instantiated as p001, p002, ... for all primers"""
import os as _os

from pydna.parsers import parse_primers as _parse_primers

list_primers = _parse_primers(_os.environ["pydna_primers"])[::-1]

# TODO: warn if no primers found!
# TODO: use logger

dict_primers = dict((p.id, p) for p in list_primers)


for _i, _p in enumerate(list_primers):
    globals()["p{:03d}".format(_i)] = _p


def append_primer_list(primers: list):
    # TODO: implement!
    raise NotImplementedError("still todo")


if __name__ == "__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"] = cache
