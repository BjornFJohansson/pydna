#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s
"""
class PrimerDict(dict):
    """Special dictionary for Primers.

    This dictionary is not meant to be used directly, but is returned
    from the primerdict_number_key function. It keeps track of accessed
    keys through the accessed_keys property.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.accessed_keys = []

    def __getitem__(self, key):
        """docstring."""
        if key not in self.accessed_keys:
            self.accessed_keys.append(key)
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        """docstring."""
        try:
            new = super().__getitem__(key)
        except KeyError:
            raise ValueError(f"key {key} does not exist.")
        else:
            if new.seq != self[key].seq:
                raise ValueError
        if key not in self.accessed_keys:
            self.accessed_keys.append(key)

    def assign_numbers_to_new_primers(self, lst):
        """docstring."""
        new = []
        found = []
        oldstrs = [str(p.seq).upper() for p in self.values()]
        no = len(oldstrs)
        for p in lst[::-1]:
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
                found.append(self[i])
        new = new[::-1]
        return _pretty_str("\n".join([p.format("fasta") for p in new]))

    def pydna_code_from_list(self, lst):
        """docstring."""
        keys = []
        prstrs = [str(p.seq).upper() for p in self.values()]
        for p in lst:
            try:
                i = prstrs.index(str(p.seq).upper())
            except ValueError:
                pass
            else:
                keys.append(i)
        return self.pydna_code_from_keys(keys)

    def pydna_code_from_keys(self, keys=None):
        """docstring."""
        keys = keys or list(dict.fromkeys(self.accessed_keys))
        msg = ", ".join(f"p[{i}]" for i in keys)
        msg += " = parse_primers('''\n\n"
        msg += "\n".join(self[i].format("fasta") for i in keys)
        msg += "\n''')"
        return _pretty_str(msg)


def primerdict_number_key():
    """docstring."""
    return PrimerDict((i, p) for i, p in enumerate(primerlist()))


def primerdict_sequence_key():
    """docstring."""
    return dict((str(p.seq).upper(), p) for p in primerlist())


def pick_from_list(primers: list) -> list:
    """docstring."""
    found = []
    pd = primerdict_sequence_key()
    for p in primers[::-1]:
        key = (p.seq).upper()
        result = pd.get(key)
        if result:
            found.append(result)
        else:
            raise ValueError(f"{key} not found in list.")
    return found


def pick_from_list_code_gen(primers: list) -> _pretty_str:
    """docstring."""
    pl = tuple(str(p.seq).upper() for p in primerlist())
    result = "from pydna.myprimers import primerlist\n\n"
    result += "p = primerlist()\n\n"
    msg = ""
    for p in primers:
        no = pl.index(str(p.seq).upper())
        msg += f"p[{no}], "
    msg = msg.rstrip(", ")
    result += f"primerlist = [{msg}]\n"
    return _pretty_str(result)


def code_for_literal_primers(primers: list) -> _pretty_str:
    """docstring."""
    pl = tuple(str(p.seq).upper() for p in primerlist())
    result = "from pydna.parsers import parse_primers\n\n"
    result += "if not \"p\" in locals(): p = {}\n\n"
    result += "new_primers = "
    msg = ""

    result = "from pydna.parsers import parse_primers\n\n"
    result += "p = {}\n\n"
    for name in names:
        result += f"{name}, "

    result += "= parse_primers(\"\"\"\n\n"

    for p in primerlist:
        result += p.format("fasta")+"\n"

    result += "\"\"\")"

    return result


def prepend_primerlist(newprimers: list,
                       oldprimers: list = None) -> _pretty_str:
    """docstring."""
    new = []
    found = []
    if not oldprimers:
        oldprimers = primerlist()
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
    return _pretty_str("\n".join([p.format("fasta") for p in new]))