#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Provides the Dseq class for handling double stranded DNA sequences.
Dseq is a subclass of :class:`Bio.Seq.Seq`. The Dseq class
is mostly useful as a part of the :class:`pydna.dseqrecord.Dseqrecord` class which
can hold more meta data.

The Dseq class support the notion of circular and linear DNA topology."""

from collections.abc import Iterable as _Iterable

import copy as _copy
import itertools as _itertools

import sys as _sys
import math as _math

from Bio.Seq import Seq as _Seq

# from Bio.Alphabet import generic_dna as _generic_dna
from pydna._pretty import pretty_str as _pretty_str
from pydna.utils import seguid as _seg
from pydna.utils import rc as _rc
from pydna.utils import complement as _complement
from pydna.utils import flatten as _flatten
from pydna.common_sub_strings import common_sub_strings as _common_sub_strings

from Bio.Restriction import RestrictionBatch as _RestrictionBatch
from Bio.Restriction import CommOnly


class Dseq(_Seq):
    """Dseq is a class designed to hold information for a double stranded
    DNA fragment. Dseq also holds information describing the topology of
    the DNA fragment (linear or circular).

    Parameters
    ----------

    watson : str
        a string representing the watson (sense) DNA strand.

    crick : str, optional
        a string representing the crick (antisense) DNA strand.

    ovhg : int, optional
        A positive or negative number to describe the stagger between the watson and crick strands.
        see below for a detailed explanation.

    linear : bool, optional
        True indicates that sequence is linear, False that it is circular.

    circular : bool, optional
        True indicates that sequence is circular, False that it is linear.


    Examples
    --------

    Dseq is a subclass of the Biopython Seq object. It stores two
    strings representing the watson (sense) and crick(antisense) strands.
    two properties called linear and circular, and a numeric value ovhg
    (overhang) describing the stagger for the watson and crick strand
    in the 5' end of the fragment.

    The most common usage is probably to create a Dseq object as a
    part of a Dseqrecord object (see :class:`pydna.dseqrecord.Dseqrecord`).

    There are three ways of creating a Dseq object directly:

    Only one argument (string):

    >>> from pydna.dseq import Dseq
    >>> Dseq("aaa")
    Dseq(-3)
    aaa
    ttt

    The given string will be interpreted as the watson strand of a
    blunt, linear double stranded sequence object. The crick strand
    is created automatically from the watson strand.

    Two arguments (string, string):

    >>> from pydna.dseq import Dseq
    >>> Dseq("gggaaat","ttt")
    Dseq(-7)
    gggaaat
       ttt

    If both watson and crick are given, but not ovhg an attempt
    will be made to find the best annealing between the strands.
    There are limitations to this! For long fragments it is quite
    slow. The length of the annealing sequences have to be at least
    half the length of the shortest of the strands.

    Three arguments (string, string, ovhg=int):

    The ovhg parameter controls the stagger at the five prime end::

        ovhg = -2

        XXXXX
          XXX

        ovhg = -1

        XXXXX
         XXXX

        ovhg = 0

        XXXXX
        XXXXX

        ovhg = 1

         XXXX
        XXXXX

        ovhg = 2

          XXX
        XXXXX

    Example of creating Dseq objects with different amounts of stagger:

    >>> Dseq(watson="agt",crick="actta",ovhg=-2)
    Dseq(-7)
    agt
      attca
    >>> Dseq(watson="agt",crick="actta",ovhg=-1)
    Dseq(-6)
    agt
     attca
    >>> Dseq(watson="agt",crick="actta",ovhg=0)
    Dseq(-5)
    agt
    attca
    >>> Dseq(watson="agt",crick="actta",ovhg=1)
    Dseq(-5)
     agt
    attca
    >>> Dseq(watson="agt",crick="actta",ovhg=2)
    Dseq(-5)
      agt
    attca

    If the ovhg parameter is specified a crick strand also needs to be supplied,
    otherwise an exception is raised.

    >>> Dseq(watson="agt",ovhg=2)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/usr/local/lib/python2.7/dist-packages/pydna_/dsdna.py", line 169, in __init__
        else:
    ValueError: ovhg defined without crick strand!


    The shape of the fragment is set by either:

    1. linear   = False, True
    2. circular = True, False

    Note that both ends of the DNA fragment has to be compatible to set
    circular = True (or linear = False).


    >>> Dseq("aaa","ttt")
    Dseq(-3)
    aaa
    ttt
    >>> Dseq("aaa","ttt",ovhg=0)
    Dseq(-3)
    aaa
    ttt
    >>> Dseq("aaa","ttt",ovhg=1)
    Dseq(-4)
     aaa
    ttt
    >>> Dseq("aaa","ttt",ovhg=-1)
    Dseq(-4)
    aaa
     ttt
    >>> Dseq("aaa", "ttt", linear = False ,ovhg=0)
    Dseq(o3)
    aaa
    ttt
    >>> Dseq("aaa", "ttt", circular = True , ovhg=0)
    Dseq(o3)
    aaa
    ttt

    >>> a=Dseq("tttcccc","aaacccc")
    >>> a
    Dseq(-11)
        tttcccc
    ccccaaa
    >>> a.ovhg
    4

    >>> b=Dseq("ccccttt","ccccaaa")
    >>> b
    Dseq(-11)
    ccccttt
        aaacccc
    >>> b.ovhg
    -4
    >>>

    Coercing to string

    >>> str(a)
    'ggggtttcccc'

    A Dseq object can be longer that either the watson or crick strands.

    ::

        <-- length -->
        GATCCTTT
             AAAGCCTAG

        <-- length -->
              GATCCTTT
        AAAGCCCTA

    The slicing of a linear Dseq object works mostly as it does for a string.

    >>> s="ggatcc"
    >>> s[2:3]
    'a'
    >>> s[2:4]
    'at'
    >>> s[2:4:-1]
    ''
    >>> s[::2]
    'gac'
    >>> from pydna.dseq import Dseq
    >>> d=Dseq(s, linear=True)
    >>> d[2:3]
    Dseq(-1)
    a
    t
    >>> d[2:4]
    Dseq(-2)
    at
    ta
    >>> d[2:4:-1]
    Dseq(-0)
    <BLANKLINE>
    <BLANKLINE>
    >>> d[::2]
    Dseq(-3)
    gac
    ctg


    The slicing of a circular Dseq object has a slightly different meaning.


    >>> s="ggAtCc"
    >>> d=Dseq(s, circular=True)
    >>> d
    Dseq(o6)
    ggAtCc
    ccTaGg
    >>> d[4:3]
    Dseq(-5)
    CcggA
    GgccT


    The slice [X:X] produces an empty slice for a string, while this will return
    the linearized sequence starting at X:

    >>> s="ggatcc"
    >>> d=Dseq(s, circular=True)
    >>> d
    Dseq(o6)
    ggatcc
    cctagg
    >>> d[3:3]
    Dseq(-6)
    tccgga
    aggcct
    >>>


    See also
    --------
    pydna.dseqrecord.Dseqrecord

    """

    def __init__(
        self, watson, crick=None, ovhg=None, linear=None, circular=None, pos=0
    ):

        if crick is None:
            if ovhg is None:
                crick = _rc(watson)
                ovhg = 0
                self._data = watson
            else:  # ovhg given, but no crick strand
                raise ValueError("ovhg defined without crick strand!")
        else:  # crick strand given
            if ovhg is None:  # ovhg not given

                olaps = _common_sub_strings(
                    str(watson).lower(),
                    str(_rc(crick).lower()),
                    int(_math.log(len(watson)) / _math.log(4)),
                )
                try:
                    F, T, L = olaps[0]
                except IndexError:
                    raise ValueError(
                        "Could not anneal the two strands. Please provide ovhg value"
                    )
                ovhgs = [ol[1] - ol[0] for ol in olaps if ol[2] == L]
                if len(ovhgs) > 1:
                    raise ValueError(
                        "More than one way of annealing the strands. Please provide ovhg value"
                    )
                ovhg = T - F

                sns = (ovhg * " ") + _pretty_str(watson)
                asn = (-ovhg * " ") + _pretty_str(_rc(crick))

                self._data = "".join(
                    [
                        a.strip() or b.strip()
                        for a, b in _itertools.zip_longest(sns, asn, fillvalue=" ")
                    ]
                )

            else:  # ovhg given
                if ovhg == 0:
                    if len(watson) == len(crick):
                        self._data = watson
                    elif len(watson) > len(crick):
                        self._data = watson
                    else:
                        self._data = watson + _rc(crick[: len(crick) - len(watson)])
                elif ovhg > 0:
                    if ovhg + len(watson) > len(crick):
                        self._data = _rc(crick[-ovhg:]) + watson
                    else:
                        self._data = (
                            _rc(crick[-ovhg:])
                            + watson
                            + _rc(crick[: len(crick) - ovhg - len(watson)])
                        )
                else:  # ovhg < 0
                    if -ovhg + len(crick) > len(watson):
                        self._data = watson + _rc(
                            crick[: -ovhg + len(crick) - len(watson)]
                        )
                    else:
                        self._data = watson

        self._circular = (
            bool(circular)
            and bool(linear) ^ bool(circular)
            or linear == False
            and circular is None
        )
        self._linear = not self._circular
        self.watson = _pretty_str(watson)
        self.crick = _pretty_str(crick)
        # self.length = max(len(watson)+max(0,ovhg), len(crick)+max(0,-ovhg))
        self.length = len(self._data)
        self._ovhg = ovhg
        self.pos = pos
        self._data = self._data
        # self.alphabet = _generic_dna

    @classmethod
    def quick(cls, watson: str, crick: str, ovhg=0, linear=True, circular=False, pos=0):
        obj = cls.__new__(cls)  # Does not call __init__
        obj.watson = _pretty_str(watson)
        obj.crick = _pretty_str(crick)
        obj._ovhg = ovhg
        obj._circular = circular
        obj._linear = linear
        obj.length = max(len(watson) + max(0, ovhg), len(crick) + max(0, -ovhg))
        obj.pos = pos
        obj._data = (
            _rc(crick[-max(0, ovhg) or len(crick) :])
            + watson
            + _rc(crick[: max(0, len(crick) - ovhg - len(watson))])
        )
        # obj.alphabet = _generic_dna
        return obj

    @classmethod
    def from_string(cls, dna: str, *args, linear=True, circular=False, **kwargs):
        obj = cls.__new__(cls)  # Does not call __init__
        obj.watson = _pretty_str(dna)
        obj.crick = _pretty_str(_rc(dna))
        obj._ovhg = 0
        obj._circular = circular
        obj._linear = linear
        obj.length = len(dna)
        obj.pos = 0
        obj._data = dna
        # obj.alphabet = _generic_dna
        return obj

    @property
    def ovhg(self):
        """The ovhg property. This cannot be set directly, but is a
        consequence of how the watson and crick strands anneal to
        each other"""
        return self._ovhg

    @property
    def linear(self):
        """The linear property can not be set directly.
        Use an empty slice [:] to create a linear object."""
        return self._linear

    @property
    def circular(self):
        """The circular property can not be set directly.
        Use :meth:`looped` to create a circular Dseq object"""
        return self._circular

    def mw(self):
        """This method returns the molecular weight of the DNA molecule
        in g/mol. The following formula is used::

               MW = (A x 313.2) + (T x 304.2) +
                    (C x 289.2) + (G x 329.2) +
                    (N x 308.9) + 79.0
        """
        nts = (self.watson + self.crick).lower()

        return (
            sum(
                [313.2 for n in nts if n == "a"]
                + [304.2 for n in nts if n == "t"]
                + [289.2 for n in nts if n == "c"]
                + [329.2 for n in nts if n == "g"]
                + [308.9 for n in nts if n == "n"]
            )
            + 79
        )

    def upper(self):
        """Return an upper case copy of the sequence.

        >>> from pydna.dseq import Dseq
        >>> my_seq = Dseq("aAa")
        >>> my_seq
        Dseq(-3)
        aAa
        tTt
        >>> my_seq.upper()
        Dseq(-3)
        AAA
        TTT

        Returns
        -------
        Dseq
            Dseq object in uppercase

        See also
        --------
        pydna.dseq.Dseq.lower

        """
        return self.quick(
            self.watson.upper(),
            self.crick.upper(),
            ovhg=self.ovhg,
            linear=self.linear,
            circular=self.circular,
            pos=self.pos,
        )

    def lower(self):
        """Return a lower case copy of the sequence.

        >>> from pydna.dseq import Dseq
        >>> my_seq = Dseq("aAa")
        >>> my_seq
        Dseq(-3)
        aAa
        tTt
        >>> my_seq.lower()
        Dseq(-3)
        aaa
        ttt

        Returns
        -------
        Dseq
            Dseq object in lowercase

        See also
        --------
        pydna.dseq.Dseq.upper
        """
        return self.quick(
            self.watson.lower(),
            self.crick.lower(),
            ovhg=self.ovhg,
            linear=self.linear,
            circular=self.circular,
            pos=self.pos,
        )

    def find(self, sub, start=0, end=_sys.maxsize):
        """This method behaves like the python string method of the same name.

        Returns an integer, the index of the first occurrence of substring
        argument sub in the (sub)sequence given by [start:end].

        Returns -1 if the subsequence is NOT found.

        Parameters
        ----------

        sub : string or Seq object
            a string or another Seq object to look for.

        start : int, optional
            slice start.

        end : int, optional
            slice end.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> seq = Dseq("atcgactgacgtgtt")
        >>> seq
        Dseq(-15)
        atcgactgacgtgtt
        tagctgactgcacaa
        >>> seq.find("gac")
        3
        >>> seq = Dseq(watson="agt",crick="actta",ovhg=-2)
        >>> seq
        Dseq(-7)
        agt
          attca
        >>> seq.find("taa")
        2
        """

        if self.linear:
            return _Seq.find(self, sub, start, end)

        return (_pretty_str(self) + _pretty_str(self)).find(sub, start, end)

    def __getitem__(self, sl):
        """Returns a subsequence. This method is used by the slice notation"""

        if self.linear:

            x = len(self.crick) - self._ovhg - len(self.watson)

            sns = (self._ovhg * " " + self.watson + x * " ")[sl]
            asn = (-self._ovhg * " " + self.crick[::-1] + -x * " ")[sl]

            ovhg = max(
                (len(sns) - len(sns.lstrip()), -len(asn) + len(asn.lstrip())), key=abs
            )

            return Dseq(sns.strip(), asn[::-1].strip(), ovhg=ovhg, linear=True)
        else:
            sl = slice(sl.start or 0, sl.stop or len(self), sl.step)
            if sl.start > len(self) or sl.stop > len(self):
                return Dseq("")
            if sl.start < sl.stop:
                return Dseq(
                    self.watson[sl], self.crick[::-1][sl][::-1], ovhg=0, linear=True
                )
            else:
                try:
                    stp = abs(sl.step)
                except TypeError:
                    stp = 1
                start = sl.start
                stop = sl.stop

                w = (
                    self.watson[(start or len(self)) :: stp]
                    + self.watson[: (stop or 0) : stp]
                )
                c = (
                    self.crick[len(self) - stop :: stp]
                    + self.crick[: len(self) - start : stp]
                )

                return Dseq(w, c, ovhg=0, linear=True)

    def __eq__(self, other):
        """Compare to another Dseq object OR an object that implements
        watson, crick and ovhg properties. This comparison is case
        insensitive.

        """
        try:
            same = (
                other.watson.lower() == self.watson.lower()
                and other.crick.lower() == self.crick.lower()
                and other.ovhg == self._ovhg
                and self.circular == other.circular
            )
            # Also test for alphabet ?
        except AttributeError:
            same = False
        return same

    def __repr__(self):
        """Returns a representation of the sequence, truncated if
        longer than 30 bp"""

        if len(self) > 30:

            if self._ovhg > 0:
                d = self.crick[-self._ovhg :][::-1]
                hej = len(d)
                if len(d) > 10:
                    d = "{}..{}".format(d[:4], d[-4:])
                a = len(d) * " "

            elif self._ovhg < 0:
                a = self.watson[: max(0, -self._ovhg)]
                hej = len(a)
                if len(a) > 10:
                    a = "{}..{}".format(a[:4], a[-4:])
                d = len(a) * " "
            else:
                a = ""
                d = ""
                hej = 0

            x = self._ovhg + len(self.watson) - len(self.crick)

            if x > 0:
                c = self.watson[len(self.crick) - self._ovhg :]
                y = len(c)
                if len(c) > 10:
                    c = "{}..{}".format(c[:4], c[-4:])
                f = len(c) * " "
            elif x < 0:
                f = self.crick[:-x][::-1]
                y = len(f)
                if len(f) > 10:
                    f = "{}..{}".format(f[:4], f[-4:])
                c = len(f) * " "
            else:
                c = ""
                f = ""
                y = 0

            L = len(self) - hej - y
            x1 = -min(0, self._ovhg)
            x2 = x1 + L
            x3 = -min(0, x)
            x4 = x3 + L

            b = self.watson[x1:x2]
            e = self.crick[x3:x4][::-1]

            if len(b) > 10:
                b = "{}..{}".format(b[:4], b[-4:])
                e = "{}..{}".format(e[:4], e[-4:])

            return _pretty_str(
                "{klass}({top}{size})\n" "{a}{b}{c}\n" "{d}{e}{f}"
            ).format(
                klass=self.__class__.__name__,
                top={True: "-", False: "o"}[self.linear],
                size=len(self),
                a=a,
                b=b,
                c=c,
                d=d,
                e=e,
                f=f,
            )

        else:
            return _pretty_str(
                "{}({}{})\n{}\n{}".format(
                    self.__class__.__name__,
                    {True: "-", False: "o"}[self.linear],
                    len(self),
                    self._ovhg * " " + self.watson,
                    -self._ovhg * " " + self.crick[::-1],
                )
            )

    def reverse_complement(self):
        """Returns a Dseq object where watson and crick have switched
        places. This represents the same double stranded sequence.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("catcgatc")
        >>> a
        Dseq(-8)
        catcgatc
        gtagctag
        >>> b=a.reverse_complement()
        >>> b
        Dseq(-8)
        gatcgatg
        ctagctac
        >>>

        """
        ovhg = len(self.watson) - len(self.crick) + self._ovhg
        return Dseq.quick(
            self.crick,
            self.watson,
            ovhg=ovhg,
            linear=self.linear,
            circular=self.circular,
        )

    rc = reverse_complement  # alias for reverse_complement

    def shifted(self, shift):
        """Returns a shifted version of a circular Dseq object."""
        if self.linear:
            raise TypeError("DNA is not circular.")
        shift = shift % len(self)
        if not shift:
            return self
        else:
            return (self[shift:] + self[:shift]).looped()

    def looped(self):
        """Returns a circularized Dseq object. This can only be done if the
        two ends are compatible, otherwise a TypeError is raised.

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("catcgatc")
        >>> a
        Dseq(-8)
        catcgatc
        gtagctag
        >>> a.looped()
        Dseq(o8)
        catcgatc
        gtagctag
        >>> a.T4("t")
        Dseq(-8)
        catcgat
         tagctag
        >>> a.T4("t").looped()
        Dseq(o7)
        catcgat
        gtagcta
        >>> a.T4("a")
        Dseq(-8)
        catcga
          agctag
        >>> a.T4("a").looped()
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
          File "/usr/local/lib/python2.7/dist-packages/pydna/dsdna.py", line 357, in looped
            if type5 == type3 and str(sticky5) == str(rc(sticky3)):
        TypeError: DNA cannot be circularized.
        5' and 3' sticky ends not compatible!
        >>>

        """
        if self.circular:
            return self
        type5, sticky5 = self.five_prime_end()
        type3, sticky3 = self.three_prime_end()
        if type5 == type3 and str(sticky5) == str(_rc(sticky3)):
            nseq = Dseq.quick(
                self.watson,
                self.crick[-self._ovhg :] + self.crick[: -self._ovhg],
                ovhg=0,
                linear=False,
                circular=True,
            )
            assert len(nseq.crick) == len(nseq.watson)
            return nseq
        else:
            raise TypeError(
                "DNA cannot be circularized.\n" "5' and 3' sticky ends not compatible!"
            )

    def tolinear(self):  # pragma: no cover
        """Returns a blunt, linear copy of a circular Dseq object. This can
        only be done if the Dseq object is circular, otherwise a
        TypeError is raised.

        This method is deprecated, use slicing instead. See example below.

        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("catcgatc", circular=True)
        >>> a
        Dseq(o8)
        catcgatc
        gtagctag
        >>> a[:]
        Dseq(-8)
        catcgatc
        gtagctag
        >>>

        """
        import warnings as _warnings
        from pydna import _PydnaDeprecationWarning

        _warnings.warn(
            "tolinear method is obsolete; "
            "please use obj[:] "
            "instead of obj.tolinear().",
            _PydnaDeprecationWarning,
        )
        if self.linear:
            raise TypeError("DNA is not circular.\n")
        selfcopy = _copy.copy(self)
        selfcopy._linear = True
        selfcopy._circular = False
        return selfcopy  # self.__class__(self.watson, linear=True)

    def five_prime_end(self):
        """Returns a tuple describing the structure of the 5' end of
        the DNA fragment

        Examples
        --------
        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.five_prime_end()
        ('blunt', '')
        >>> a=Dseq("aaa", "ttt", ovhg=1)
        >>> a
        Dseq(-4)
         aaa
        ttt
        >>> a.five_prime_end()
        ("3'", 't')
        >>> a=Dseq("aaa", "ttt", ovhg=-1)
        >>> a
        Dseq(-4)
        aaa
         ttt
        >>> a.five_prime_end()
        ("5'", 'a')
        >>>

        See also
        --------
        pydna.dseq.Dseq.three_prime_end

        """
        if self.watson and not self.crick:
            return "5'", self.watson.lower()
        if not self.watson and self.crick:
            return "3'", self.crick.lower()
        if self._ovhg < 0:
            sticky = self.watson[: -self._ovhg].lower()
            type_ = "5'"
        elif self._ovhg > 0:
            sticky = self.crick[-self._ovhg :].lower()
            type_ = "3'"
        else:
            sticky = ""
            type_ = "blunt"
        return type_, sticky

    def three_prime_end(self):
        """Returns a tuple describing the structure of the 5' end of
        the DNA fragment

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.three_prime_end()
        ('blunt', '')
        >>> a=Dseq("aaa", "ttt", ovhg=1)
        >>> a
        Dseq(-4)
         aaa
        ttt
        >>> a.three_prime_end()
        ("3'", 'a')
        >>> a=Dseq("aaa", "ttt", ovhg=-1)
        >>> a
        Dseq(-4)
        aaa
         ttt
        >>> a.three_prime_end()
        ("5'", 't')
        >>>

        See also
        --------
        pydna.dseq.Dseq.five_prime_end

        """

        ovhg = len(self.watson) - len(self.crick) + self._ovhg

        if ovhg < 0:
            sticky = self.crick[:-ovhg].lower()
            type_ = "5'"
        elif ovhg > 0:
            sticky = self.watson[-ovhg:].lower()
            type_ = "3'"
        else:
            sticky = ""
            type_ = "blunt"
        return type_, sticky

    def __add__(self, other):
        """Simulates ligation between two DNA fragments.

        Add other Dseq object at the end of the sequence.
        Type error is raised if any of the points below are fulfilled:

        * one or more objects are circular
        * if three prime sticky end of self is not the same type
          (5' or 3') as the sticky end of other
        * three prime sticky end of self complementary with five
          prime sticky end of other.

        Phosphorylation and dephosphorylation is not considered.

        DNA is allways presumed to have the necessary 5' phospate
        group necessary for ligation.

        """
        # test for circular DNA
        if self.circular:
            raise TypeError("circular DNA cannot be ligated!")
        try:
            if other.circular:
                raise TypeError("circular DNA cannot be ligated!")
        except AttributeError:
            pass

        self_type, self_tail = self.three_prime_end()
        other_type, other_tail = other.five_prime_end()

        if self_type == other_type and str(self_tail) == str(_rc(other_tail)):

            answer = Dseq.quick(
                self.watson + other.watson, other.crick + self.crick, self._ovhg
            )
        elif not self:
            answer = _copy.copy(other)
        elif not other:
            answer = _copy.copy(self)
        else:
            raise TypeError("sticky ends not compatible!")
        return answer

    def __mul__(self, number):
        if not isinstance(number, int):
            raise TypeError(
                "TypeError: can't multiply Dseq by non-int of type {}".format(
                    type(number)
                )
            )
        if number <= 0:
            return self.__class__("")
        new = _copy.copy(self)
        for i in range(number - 1):
            new += self
        return new

    def _fill_in_five_prime(self, nucleotides):
        stuffer = ""
        type, se = self.five_prime_end()
        if type == "5'":
            for n in _rc(se):
                if n in nucleotides:
                    stuffer += n
                else:
                    break
        return self.crick + stuffer, self._ovhg + len(stuffer)

    def _fill_in_three_prime(self, nucleotides):
        stuffer = ""
        type, se = self.three_prime_end()
        if type == "5'":
            for n in _rc(se):
                if n in nucleotides:
                    stuffer += n
                else:
                    break
        return self.watson + stuffer

    def fill_in(self, nucleotides=None):
        """Fill in of five prime protruding end with a DNA polymerase
        that has only DNA polymerase activity (such as exo-klenow [#]_)
        and any combination of A, G, C or T. Default are all four
        nucleotides together.

        Parameters
        ----------

        nucleotides : str

        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.fill_in()
        Dseq(-3)
        aaa
        ttt
        >>> b=Dseq("caaa", "cttt")
        >>> b
        Dseq(-5)
        caaa
         tttc
        >>> b.fill_in()
        Dseq(-5)
        caaag
        gtttc
        >>> b.fill_in("g")
        Dseq(-5)
        caaag
        gtttc
        >>> b.fill_in("tac")
        Dseq(-5)
        caaa
         tttc
        >>> c=Dseq("aaac", "tttg")
        >>> c
        Dseq(-5)
         aaac
        gttt
        >>> c.fill_in()
        Dseq(-5)
         aaac
        gttt
        >>>

        References
        ----------
        .. [#] http://en.wikipedia.org/wiki/Klenow_fragment#The_exo-_Klenow_fragment

        """
        if not nucleotides:
            nucleotides = "GATCRYWSMKHBVDN"
        nucleotides = set(nucleotides.lower() + nucleotides.upper())
        crick, ovhg = self._fill_in_five_prime(nucleotides)
        watson = self._fill_in_three_prime(nucleotides)
        return Dseq(watson, crick, ovhg)

    def mung(self):
        """
        Simulates treatment a nuclease with 5'-3' and 3'-5' single
        strand specific exonuclease activity (such as mung bean nuclease [#]_)

        ::

             ggatcc    ->     gatcc
              ctaggg          ctagg

              ggatcc   ->      ggatc
             tcctag            cctag

         >>> from pydna.dseq import Dseq
         >>> b=Dseq("caaa", "cttt")
         >>> b
         Dseq(-5)
         caaa
          tttc
         >>> b.mung()
         Dseq(-3)
         aaa
         ttt
         >>> c=Dseq("aaac", "tttg")
         >>> c
         Dseq(-5)
          aaac
         gttt
         >>> c.mung()
         Dseq(-3)
         aaa
         ttt



        References
        ----------
        .. [#] http://en.wikipedia.org/wiki/Mung_bean_nuclease


        """
        return Dseq(
            self.watson[
                max(0, -self.ovhg) : min(len(self.watson), len(self.crick) - self.ovhg)
            ]
        )

    def T4(self, nucleotides=None):
        """Fill in five prime protruding ends and chewing back
        three prime protruding ends by a DNA polymerase providing both
        5'-3' DNA polymerase activity and 3'-5' nuclease acitivty
        (such as T4 DNA polymerase). This can be done in presence of any
        combination of the four A, G, C or T. Removing one or more nucleotides
        can facilitate engineering of sticky ends. Default are all four nucleotides together.

        Parameters
        ----------
        nucleotides : str


        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("gatcgatc")
        >>> a
        Dseq(-8)
        gatcgatc
        ctagctag
        >>> a.T4()
        Dseq(-8)
        gatcgatc
        ctagctag
        >>> a.T4("t")
        Dseq(-8)
        gatcgat
         tagctag
        >>> a.T4("a")
        Dseq(-8)
        gatcga
          agctag
        >>> a.T4("g")
        Dseq(-8)
        gatcg
           gctag
        >>>

        """

        if not nucleotides:
            nucleotides = "GATCRYWSMKHBVDN"
        nucleotides = set(nucleotides.lower() + nucleotides.upper())
        type, se = self.five_prime_end()
        if type == "5'":
            crick, ovhg = self._fill_in_five_prime(nucleotides)
        else:
            if type == "3'":
                ovhg = 0
                crick = self.crick[: -len(se)]
            else:
                ovhg = 0
                crick = self.crick
        x = len(crick) - 1
        while x >= 0:
            if crick[x] in nucleotides:
                break
            x -= 1
        ovhg = x - len(crick) + 1 + ovhg
        crick = crick[: x + 1]
        if not crick:
            ovhg = 0
        watson = self.watson
        type, se = self.three_prime_end()
        if type == "5'":
            watson = self._fill_in_three_prime(nucleotides)
        else:
            if type == "3'":
                watson = self.watson[: -len(se)]
        x = len(watson) - 1
        while x >= 0:
            if watson[x] in nucleotides:
                break
            x -= 1
        watson = watson[: x + 1]
        return Dseq(watson, crick, ovhg)

    t4 = T4  # alias for the T4 method.

    def no_cutters(self, batch=CommOnly):
        """Returns the enzymes in a RestrictionBatch that do **not**
        cut the sequence."""
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if not sitelist}
        return _RestrictionBatch(ncut)

    def unique_cutters(self, batch=CommOnly):
        """Returns the enzymes in a RestrictionBatch that cut the sequence
        exactly once."""
        return self.n_cutters(n=1, batch=batch)

    once_cutters = unique_cutters  # alias for unique_cutters

    def twice_cutters(self, batch=CommOnly):
        """Returns the enzymes in a RestrictionBatch that cut the sequence
        exactly twice."""
        return self.n_cutters(n=2, batch=batch)

    def n_cutters(self, n=3, batch=CommOnly):
        """Returns the enzymes in a RestrictionBatch that cut the sequence
        n times."""
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if len(sitelist) == n}
        return _RestrictionBatch(ncut)

    def cutters(self, batch=CommOnly):
        """Returns the enzymes in a RestrictionBatch that cut the sequence
        at least once."""
        ana = batch.search(self)
        ncut = {enz: sitelist for (enz, sitelist) in ana.items() if sitelist}
        return _RestrictionBatch(ncut)

    def seguid(self):
        """Returns the SEGUID for the Dseq. The definition
        varies with the amount of stagger between the sequences."""
        rc_ovhg = len(self.watson) - len(self.crick) + self._ovhg
        if self._ovhg == rc_ovhg == 0:
            return _seg(min(self.watson, self.crick))
        if self._ovhg < rc_ovhg:
            w = self.watson
            c = self.crick
            o = self._ovhg
        elif self._ovhg > rc_ovhg:
            w = self.crick
            c = self.watson
            o = rc_ovhg
        elif self._ovhg == rc_ovhg:
            w, c = sorted((self.watson, self.crick))
            o = self._ovhg
        return _seg(_pretty_str(o) + w + "|" + c)

    def isblunt(self):
        """Returns True if Dseq is blunt and false if staggered or circular.

        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> a=Dseq("gat")
        >>> a
        Dseq(-3)
        gat
        cta
        >>> a.isblunt()
        True
        >>> a=Dseq("gat", "atcg")
        >>> a
        Dseq(-4)
         gat
        gcta
        >>> a.isblunt()
        False
        >>> a=Dseq("gat", "gatc")
        >>> a
        Dseq(-4)
        gat
        ctag
        >>> a.isblunt()
        False
        >>> a=Dseq("gat", circular=True)
        >>> a
        Dseq(o3)
        gat
        cta
        >>> a.isblunt()
        False"""

        return ( self._ovhg==0 and
                len(self.watson) == len(self.crick)
                and self._linear)

    def cut(self, *enzymes):
        """Returns a list of linear Dseq fragments produced in the digestion.
        If there are no cuts, an empty list is returned.

        Parameters
        ----------

        enzymes : enzyme object or iterable of such objects
            A Bio.Restriction.XXX restriction objects or iterable.

        Returns
        -------
        frags : list
            list of Dseq objects formed by the digestion


        Examples
        --------

        >>> from pydna.dseq import Dseq
        >>> seq=Dseq("ggatccnnngaattc")
        >>> seq
        Dseq(-15)
        ggatccnnngaattc
        cctaggnnncttaag
        >>> from Bio.Restriction import BamHI,EcoRI
        >>> type(seq.cut(BamHI))
        <class 'tuple'>
        >>> for frag in seq.cut(BamHI): print(repr(frag))
        Dseq(-5)
        g
        cctag
        Dseq(-14)
        gatccnnngaattc
            gnnncttaag
        >>> seq.cut(EcoRI, BamHI) ==  seq.cut(BamHI, EcoRI)
        True
        >>> a,b,c = seq.cut(EcoRI, BamHI)
        >>> a+b+c
        Dseq(-15)
        ggatccnnngaattc
        cctaggnnncttaag
        >>>

        """

        pad = "n" * 50

        if self.linear:
            dsseq = self.mung()
        else:
            dsseq = Dseq.from_string(self._data, linear=True, circular=False)

        if len(enzymes) == 1 and hasattr(
            enzymes[0], "intersection"
        ):  # RestrictionBatch
            enzymecuts = []
            for e in enzymes[0]:
                # cuts = e.search(dsseq+dsseq[:e.size-1] if self.circular else dsseq)
                cuts = e.search(
                    _Seq(pad + dsseq.watson + dsseq.watson[: e.size - 1] + pad)
                    if self.circular
                    else dsseq
                )
                enzymecuts.append((cuts, e))
            enzymecuts.sort()
            enzymes = [e for (c, e) in enzymecuts if c]
        else:
            enzymes = [
                e
                for e in list(dict.fromkeys(_flatten(enzymes)))
                if e.search(
                    _Seq(pad + dsseq.watson + dsseq.watson[: e.size - 1] + pad)
                    if self.circular
                    else dsseq
                )
            ]  # flatten

        if not enzymes:
            return ()

        if self.linear:
            frags = [self]
        else:
            l = len(self)
            for e in enzymes:
                wpos = [
                    x - len(pad) - 1
                    for x in e.search(
                        _Seq(pad + self.watson + self.watson[: e.size - 1]) + pad
                    )
                ][::-1]
                cpos = [
                    x - len(pad) - 1
                    for x in e.search(
                        _Seq(pad + self.crick + self.crick[: e.size - 1]) + pad
                    )
                ][::-1]

                for w, c in _itertools.product(wpos, cpos):
                    if w % len(self) == (self.length - c + e.ovhg) % len(self):
                        frags = [
                            Dseq(
                                self.watson[w % l :] + self.watson[: w % l],
                                self.crick[c % l :] + self.crick[: c % l],
                                ovhg=e.ovhg,
                                pos=w,
                            )
                        ]
                        break
                else:
                    continue
                break

        newfrags = []

        for enz in enzymes:
            for frag in frags:

                ws = [x - 1 for x in enz.search(_Seq(frag.watson) + "N")]
                cs = [x - 1 for x in enz.search(_Seq(frag.crick) + "N")]

                sitepairs = [
                    (sw, sc)
                    for sw, sc in _itertools.product(ws, cs[::-1])
                    if (
                        sw + max(0, frag.ovhg) - max(0, enz.ovhg)
                        == len(frag.crick) - sc - min(0, frag.ovhg) + min(0, enz.ovhg)
                    )
                ]

                sitepairs.append((self.length, 0))

                w2, c1 = sitepairs[0]

                newfrags.append(
                    Dseq(
                        frag.watson[:w2], frag.crick[c1:], ovhg=frag.ovhg, pos=frag.pos
                    )
                )

                for (w1, c2), (w2, c1) in zip(sitepairs[:-1], sitepairs[1:]):
                    newfrags.append(
                        Dseq(
                            frag.watson[w1:w2],
                            frag.crick[c1:c2],
                            ovhg=enz.ovhg,
                            pos=frag.pos + w1 - max(0, enz.ovhg),
                        )
                    )

            frags = newfrags
            newfrags = []

        return tuple(frags)


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
