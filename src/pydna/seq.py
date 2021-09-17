#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A subclass of the Biopython SeqRecord class.

Has a number of extra methods and uses
the :class:`pydna._pretty_str.pretty_str` class instread of str for a
nicer output in the IPython shell.
"""

from pydna.codon import weights as _weights
from pydna.codon import rare_codons as _rare_codons
from pydna.codon import start as _start
from pydna.codon import stop as _stop
from pydna.codon import n_end as _n_end

from Bio.SeqUtils import seq3 as _seq3
from Bio.SeqUtils import GC as _GC
import re as _re
from Bio.Seq import Seq as _Seq
from pydna._pretty import PrettyTable as _PrettyTable

import logging as _logging

_module_logger = _logging.getLogger("pydna." + __name__)


class Seq(_Seq):
    """docstring."""

    def __init__(self, data, length=None):
        super().__init__(data, length=None)

    def gc(self):
        """Return GC content."""
        return round(_GC(self._data.upper().decode("ASCII"))/100.0, 3)

    def cai(self, organism="sce"):
        """docstring."""
        from CAI import CAI as _CAI
        return round(_CAI(self._data.upper().decode("ASCII"),
                          weights=_weights[organism]), 3)

    def rarecodons(self, organism="sce"):
        """docstring."""
        rare = _rare_codons[organism]
        s = self._data.upper().decode("ASCII")
        slices = []
        for i in range(0, len(self)//3):
            x, y = i*3, i*3+3
            trip = s[x:y]
            if trip in rare:
                slices.append(slice(x, y, 1))
        return slices

    def startcodon(self, organism="sce"):
        """docstring."""
        return _start[organism].get(self._data.upper().decode("ASCII")[:3])

    def stopcodon(self, organism="sce"):
        """docstring."""
        return _stop[organism].get(self._data.upper().decode("ASCII")[-3:])

    def express(self, organism="sce"):
        """docstring."""
        x = _PrettyTable(["cds", "len", "cai", "gc", "sta", "stp",
                          "n-end"]+_rare_codons[organism]+["rare"])
        val = []

        val.append(f"{self._data.upper().decode('ASCII')[:3]}..."
                   f"{self._data.upper().decode('ASCII')[-3:]}")
        val.append(len(self)/3)
        val.append(self.cai(organism))
        val.append(self.gc())
        val.append(self.startcodon())
        val.append(self.stopcodon())
        val.append(_n_end[organism].get(_seq3(self[3:6].translate())))
        s = self._data.upper().decode("ASCII")
        trps = [s[i*3:i*3+3] for i in range(0, len(s)//3)]
        tot = 0
        for cdn in _rare_codons[organism]:
            cnt = trps.count(cdn)
            tot += cnt
            val.append(cnt)
        val.append(round(tot/len(trps), 3))
        x.add_row(val)
        return x

    def orfs(self, minsize=30):
        """docstring."""
        orf = _re.compile(f"ATG(?:...){{{minsize},}}?(?:TAG|TAA|TGA)",
                          flags=_re.IGNORECASE)
        start = 0
        matches = []
        s = self._data.decode("ASCII")
        while True:
            match = orf.search(s, pos=start)
            if match:
                matches.append(slice(match.start(), match.end()))
                start = start + match.start() + 1
            else:
                break
        return sorted([self[sl] for sl in matches], key=len, reverse=True)
