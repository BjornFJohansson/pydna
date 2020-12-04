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


def primerlist():
    """docstring."""
    return _parse_primers(_os.environ["pydna_primers"])[::-1]


def primerdict():
    """docstring."""
    return dict((p.id, p) for p in primerlist())


# def myprimers():
#     """docstring."""
#     obj = 1
#     for _i, _p in enumerate(primerlist()):
#         name_space["p{:03d}".format(_i)] = _p
#     return obj


def prepend_primerlist(primerlist):
    """docstring."""

    primers = primerlist()

    no = len(primers)

    newprimers = []

    for i,p in zip(range(no+len(primerlist)-1, no-1, -1), primerlist):
        suff = p.id.split(str(i))[-1]
        suff.lstrip("_")
        newprimer = _copy.copy(p)
        newprimer.id = f"{i}_{suff}"
        newprimers.append(newprimer)

    return _pretty_str("\n".join([p.format("fasta") for p in newprimers]))


def check_primer_list(primerlist=primerlist):

    unique = set(str(p.seq).lower() for p in primerlist())

    print("number of primers", len(primerlist))
    print("unique primers", len(unique))
    print("no seq (n)", len([p for p in primerlist if set(p.seq.lower()) == set("n")]))

    for i, p in enumerate(primerlist):
        if not p.name.startswith(str(i)):
            print(i, p.format("tab"))

    print("names checked for correct primer number for ", i + 1, "primers")

    from collections import defaultdict

    dct = defaultdict(list)

    for u in unique:
        for p in primerlist:
            if set(u) == set("n") or set(p.seq.lower()) == set("n"):
                continue
            if u == str(p.seq).lower():
                dct[u].append(p.name)


    for seq, names in dct.items():
        if len(names) > 1:
            print("\n".join(names))
            print(seq)
            print()



if __name__ == "__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"] = cache
