#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Provides three ways to access a list of primer sequences.

The path of a text file can be specified in the pydna.ini file or by the in
´pydna_primers´ environment variable.

The file is expected to contain sequences in FASTA, Genbank or EMBL formats.
The primer list can have the format below for example:

::

    >1_second_primer
    tgatcgtcatgctgactatactat
    >0_first_primer
    ctaggatcgtagatctagctg
    ...

The primerlist funtion returns a list of :class:`pydna.primer.Primer` objects
primerdict returns a dict where the key is the id of the object.
"""

import os as _os
import copy as _copy
from pydna.parsers import parse_primers as _parse_primers
from pydna._pretty import pretty_str as _pretty_str
from collections import defaultdict as _defaultdict


def primerlist():
    """docstring."""
    lines = []
    with open(_os.environ["pydna_primers"]) as f:
        for line in f.readlines():
            if not line.startswith("#"):
                lines.append(line)
    return _parse_primers("\n".join(lines))[::-1]


def primerdict():
    """docstring."""
    return dict((p.id, p) for p in primerlist())


def prepend_primerlist(newprimers: list, oldprimers: list):
    """docstring."""
    new = []
    found = []
    no = len(oldprimers)
    oldstrs = [str(p.seq).upper() for p in oldprimers]
    for p in newprimers[::-1]:
        try:
            i = oldstrs.index(str(p.seq).upper())
        except ValueError:
            i = no + len(new)
            suff = p.id.split(str(i))[-1]
            suff.lstrip("_")
            newprimer = _copy.copy(p)
            newprimer.id = f"{i}_{suff}"
            new.append(newprimer)
        else:
            found.append(oldprimers[i])
    new = new[::-1]
    return found[::-1] or _pretty_str("\n".join([p.format("fasta") for p in new]))


def check_primer_list(primerlist: list):
    """docstring."""
    pl = primerlist

    unique_seqs = set(str(p.seq).lower() for p in pl)

    msg = f"{len(pl)} primers, {len(unique_seqs)} unique primer sequences\n"

    defined = [p for p in pl if set(p.seq.lower()) != set("n")]

    msg += f"{len(pl) - len(defined)} primer(s) without sequence (N)\n"

    for i, p in enumerate(pl):
        if not p.name.startswith(str(i)):
            msg += f"\nWrong number: {i} {p.format('tab')}"

    dct = _defaultdict(list)

    for u in unique_seqs:
        for p in defined:
            if u == str(p.seq).lower():
                dct[u].append(p.name)

    for seq, names in dct.items():
        if len(names) > 1:
            msg += " ".join(names)
            msg += f" {seq}\n"

    return _pretty_str(msg.strip())


if __name__ == "__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"] = cache
