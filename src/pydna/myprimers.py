#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""
Provides a practical way to access a list of primer sequences in a text file.

The path of a text file can be specified in the pydna.ini file or by the
´pydna_primers´ environment variable.

The file is expected to contain sequences in FASTA, Genbank or EMBL formats or
any format readable by the parse_primers function.

The primer list is expected to follow the convension below. The primer name is
expected to begin with the number.

can have the format below for example:

::

    >2_third_primer
    tgagtagtcgtagtcgtcgtat

    >1_second_primer
    tgatcgtcatgctgactatactat

    >0_first_primer
    ctaggatcgtagatctagctg
    ...

The primerlist funtion returns a list of :class:`pydna.primer.Primer` objects
primerdict returns a dict where the key is the id of the object.
"""

import os as _os
import re as _re
from typing import Iterable
from pathlib import Path
import copy as _copy
from keyword import iskeyword as _iskeyword
from pydna.parsers import parse_primers as _parse_primers
from pydna._pretty import pretty_str as _pretty_str
from collections import UserList as _UserList
from pydna.utils import open_folder as _open_folder
from builtins import __dict__ as _kw


class PrimerList(_UserList):
    """Read a text file with primers.

    The primers can be of any format readable by the parse_primers
    function. Lines beginning with # are ignored. Path defaults to
    the path given by the pydna_primers environment variable.

    The primer list does not accept new primers. Use the
    assign_numbers_to_new_primers method and paste the new
    primers at the top of the list.

    The primer list remembers the numbers of accessed primers.
    The indices of accessed primers are stored in the .accessed
    property.
    """

    def __init__(self,
                 initlist: Iterable = None,
                 path: (str, Path) = None,
                 *args,
                 identifier: str = "p",
                 **kwargs):
        if initlist:
            self.data = initlist
            self.path = None
        else:
            lines = []
            self.path = Path(path or _os.environ["pydna_primers"])
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.path, "a+") as f:
                f.seek(0)
                for line in f.readlines():
                    if not line.startswith("#"):
                        lines.append(line)
            self.data = _parse_primers("\n".join(lines))[::-1]
        # super().__init__(*args, **kwargs)
        self.accessed = []
        if (identifier.isidentifier() and not _iskeyword(identifier) and identifier not in _kw):
            self.identifier = identifier
        else:
            raise ValueError(f"{identifier} is not a valid identifier.")

    def __getitem__(self, i):
        """Save indices of accessed items."""
        if isinstance(i, slice):
            result = self.__class__(self.data[i])
            for ind in range(i.start, i.stop, i.step or 1):
                if ind not in self.accessed:
                    self.accessed.append(ind)
        else:
            try:
                result = self.data[i]
            except IndexError as e:
                raise(e)
            else:
                if i not in self.accessed:
                    self.accessed.append(i)
        return result

    def __setitem__(self, i, item):
        """Items already present can be set to the same value."""
        if abs(i) > len(self):
            raise IndexError(f"index {i} does not exist.")
        else:
            if str(item.seq).lower() != str(self.data[i].seq).lower():
                raise ValueError("Cannot change existing primer.")
        if i not in self.accessed:
            self.accessed.append(i)

    def assign_numbers_to_new_primers(self, lst: list):
        """Find new primers in lst.

        Returns a string containing new primers with their assigned
        numbers. This string can be copied and pasted to the primer
        text file.
        """
        new = []
        found = []
        oldstrs = [str(p.seq).upper() for p in self.data]
        no = len(oldstrs)
        for p in lst[::-1]:
            try:
                i = oldstrs.index(str(p.seq).upper())
            except ValueError:
                i = no + len(new)
                suffix = p.id.split(str(i))[-1]
                suffix.lstrip("_")
                newprimer = _copy.copy(p)
                newprimer.id = f"{i}_{suffix}"
                new.append(newprimer)
            else:
                found.append(self[i])
        new = new[::-1]
        return _pretty_str("\n".join([p.format("fasta") for p in new]))

    def pydna_code_from_list(self, lst: list):
        """Pydna code for a list of primer objects."""
        indices = []
        prstrs = [str(p.seq).upper() for p in self.data]
        err = None
        for p in lst:
            try:
                i = prstrs.index(str(p.seq).upper())
            except ValueError as e:
                print(f"{p.format('fasta')}")
                err = e
            else:
                indices.append(i)
        if err:
            raise ValueError("At least one primer not in list.")
        return self.pydna_code_from_indices(indices)

    def pydna_code_from_indices(self, indices: list = None):
        """Pydna code for a list of primer indices."""
        curly = "{}"
        msg = f"{self.identifier} = {curly}\n\n"
        msg += ", ".join(f"{self.identifier}[{i}]" for i in indices)
        msg += " = parse_primers('''\n\n"
        msg += "\n".join(self[i].format("fasta") for i in indices)
        msg += "\n''')"
        return _pretty_str(msg)

    def pydna_code_from_accessed(self):
        """Pydna code for acessed primers."""
        return self.pydna_code_from_indices(self.accessed)

    def open_folder(self):
        """Open folder where primer file is located."""
        _open_folder(self.path.parent)

    code = pydna_code_from_indices


def check_primer_numbers(pl: list = None):
    """Find primers whose number do not match position in list."""
    if not pl:
        pl = PrimerList()
    primers_with_wrong_number = []
    for i, p in enumerate(pl):
        if not p.name.startswith(str(i)):
            primers_with_wrong_number.append(p)
    return primers_with_wrong_number


def undefined_sequence(pl: list = None):
    """Primers in list with N or n instead of a sequence."""
    if not pl:
        pl = PrimerList()
    return [p for p in pl if _re.match("N+", str(p.seq.upper()))]


def find_duplicate_primers(pl: list = None):
    """Find a list of lists with duplicated primer sequences."""
    if not pl:
        pl = PrimerList()
    pg = {}
    for p in pl:
        pg.setdefault(str(p.seq).upper(), []).append(p)
    return [pl for ps, pl in pg.items() if len(pl) > 1]


if __name__ == "__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"] = cache
