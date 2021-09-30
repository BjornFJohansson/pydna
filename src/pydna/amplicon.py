#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
# doctest: +NORMALIZE_WHITESPACE
# doctest: +SKIP
"""This module provides the :class:`Amplicon` class for PCR simulation.
This class is not meant to be use directly but is
used by the :mod:`amplify` module"""

from pydna.tm import dbd_program as _dbd_program
from pydna.tm import program as _program
from pydna.primer import Primer as _Primer
from pydna._pretty import pretty_str as _pretty_str
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from Bio.SeqRecord import SeqRecord as _SeqRecord
import textwrap as _textwrap
import copy as _copy
import logging as _logging


_module_logger = _logging.getLogger("pydna." + __name__)


class Amplicon(_Dseqrecord):
    """The Amplicon class holds information about a PCR reaction involving two
    primers and one template. This class is used by the Anneal class and is not
    meant to be instantiated directly.

    Parameters
    ----------
    forward_primer : SeqRecord(Biopython)
        SeqRecord object holding the forward (sense) primer

    reverse_primer : SeqRecord(Biopython)
        SeqRecord object holding the reverse (antisense) primer

    template : Dseqrecord
        Dseqrecord object holding the template (circular or linear)


    """

    def __init__(
        self,
        record,
        *args,
        template=None,
        forward_primer=None,
        reverse_primer=None,
        **kwargs
    ):

        super().__init__(record, *args)
        self.template = template
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.__dict__.update(kwargs)

        # https://medium.com/@chipiga86/circular-references-without-memory-
        # leaks-and-destruction-of-objects-in-python-43da57915b8d

    @classmethod
    def from_SeqRecord(cls, record, *args, path=None, **kwargs):
        obj = super().from_SeqRecord(record, *args, **kwargs)
        obj.path = path
        return obj

    def __getitem__(self, sl):
        answer = _copy.copy(self)
        answer.seq = answer.seq.__getitem__(sl)
        # answer.seq.alphabet = self.seq.alphabet
        sr = _SeqRecord("n" * len(self))
        sr.features = self.features
        answer.features = _SeqRecord.__getitem__(sr, sl).features
        return answer

    def __repr__(self):
        """returns a short string representation of the object"""
        return "Amplicon({})".format(self.__len__())

    def _repr_pretty_(self, p, cycle):
        p.text("Amplicon({})".format(self.__len__()))

    def _repr_html_(self):
        return "Amplicon({})".format(self.__len__())

    def reverse_complement(self):
        answer = type(self)(super().reverse_complement())
        answer.template = self.template.rc()
        answer.forward_primer = self.reverse_primer
        answer.reverse_primer = self.forward_primer
        return answer

    rc = reverse_complement

    def figure(self):
        """
        This method returns a simple figure of the two primers binding
        to a part of the template.

        ::

         5tacactcaccgtctatcattatc...cgactgtatcatctgatagcac3
                                    ||||||||||||||||||||||
                                   3gctgacatagtagactatcgtg5
         5tacactcaccgtctatcattatc3
          |||||||||||||||||||||||
         3atgtgagtggcagatagtaatag...gctgacatagtagactatcgtg5



        Returns
        -------
        figure:string
             A string containing a text representation of the primers
             annealing on the template (see example above).

        """

        f = """
            {sp1}5{faz}...{raz}3
             {sp3}{rap}
            {sp3}3{rp}5
            5{fp}3
             {fap:>{fplength}}
            {sp2}3{fzc}...{rzc}5
            """.format(
            fp=self.forward_primer.seq,
            fap="|" * len(self.forward_primer.footprint),
            fplength=len(self.forward_primer.seq),
            rp=self.reverse_primer.seq[::-1],
            rap="|" * len(self.reverse_primer.footprint),
            faz=self.forward_primer.footprint,
            raz=self.reverse_primer.footprint.reverse_complement(),
            fzc=self.forward_primer.footprint.complement(),
            rzc=self.reverse_primer.footprint[::-1],
            sp1=" "
            * (len(self.forward_primer.seq) - len(self.forward_primer.footprint)),
            sp2=" "
            * (len(self.forward_primer.seq) - len(self.forward_primer.footprint)),
            sp3=" " * (3 + len(self.forward_primer.seq)),
        )
        return _pretty_str(_textwrap.dedent(f).strip("\n"))

    def set_forward_primer_footprint(self, length):
        self.forward_primer = _Primer(
            self.forward_primer.tail + self.seq[:length], footprint=length
        )

    def set_reverse_primer_footprint(self, length):
        self.reverse_primer = _Primer(
            self.reverse_primer.tail + self.seq[:length], footprint=length
        )

    def program(self):
        return _program(self)

    def dbd_program(self):
        return _dbd_program(self)


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
