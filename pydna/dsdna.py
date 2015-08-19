#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2013 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''Provides two classes, Dseq and Dseqrecord, for handling double stranded
DNA sequences. Dseq and Dseqrecord are subclasses of Biopythons
Seq and SeqRecord classes, respectively. These classes support the
notion of circular and linear DNA.

'''
import cPickle
import shelve

import copy
import datetime
import itertools
import operator
import os
import re
import StringIO
import sys
import textwrap
import math
import glob
import colorsys

from warnings import warn

from prettytable import PrettyTable

from Bio                    import SeqIO
from Bio.Alphabet.IUPAC     import IUPACAmbiguousDNA
from Bio.Seq                import Seq
from Bio.Seq                import translate
from Bio.Seq                import _maketrans
from Bio.Data.IUPACData     import ambiguous_dna_complement as amb_compl
amb_compl.update({"U":"A"})
_complement_table = _maketrans(amb_compl)
from Bio.SeqRecord          import SeqRecord
from Bio.SeqFeature         import SeqFeature
from Bio.SeqFeature         import FeatureLocation
from Bio.SeqFeature         import CompoundLocation
from Bio.SeqUtils           import GC
from Bio.GenBank            import RecordParser
from Bio.Data.CodonTable    import TranslationError

from _sequencetrace         import SequenceTraceFactory

from pydna.findsubstrings_suffix_arrays_python import common_sub_strings
from pydna.utils  import seguid  as seg
from pydna.utils  import cseguid as cseg
from pydna.pretty import pretty_str, pretty_string #, pretty_unicode

try:
    import IPython
except ImportError:
    def display(item): return item
else:
    from IPython.display import Markdown as display


def rc(sequence):
    '''returns the reverse complement of sequence (string)
    accepts mixed DNA/RNA
    '''
    return sequence.translate(_complement_table)[::-1]

class Dseq(Seq):
    '''Dseq is a class designed to hold information for a double stranded
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

    alphabet : Bio.Alphabet, optional
        Bio.Alphabet.IUPAC.IUPACAmbiguousDNA from the Biopython package is set as default.

    Examples
    --------

    Dseq is a subclass of the Biopython Seq object. It stores two
    strings representing the watson (sense) and crick(antisense) strands.
    two properties called linear and circular, and a numeric value ovhg
    (overhang) describing the stagger for the watson and crick strand
    in the 5' end of the fragment.

    The most common usage is probably to create a Dseq object as a
    part of a Dseqrecord object (see :class:`pydna.dsdna.Dseqrecord`).

    There are three ways of creating a Dseq object directly:

    Only one argument (string):

    >>> import pydna
    >>> pydna.Dseq("aaa")
    Dseq(-3)
    aaa
    ttt

    The given string will be interpreted as the watson strand of a
    blunt, linear double stranded sequence object. The crick strand
    is created automatically from the watson strand.

    Two arguments (string, string):

    >>> import pydna
    >>> pydna.Dseq("gggaaat","ttt")
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

        ovhg=-2

        XXXXX
          XXXXX

        ovhg=-1

        XXXXX
         XXXXX


        ovhg=0

        XXXXX
        XXXXX

        ovhg=1

         XXXXX
        XXXXX

        ovhg=2

          XXXXX
        XXXXX

    Example of creating Dseq objects with different amounts of stagger:

    >>> pydna.Dseq(watson="agt",crick="actta",ovhg=-2)
    Dseq(-7)
    agt
      attca
    >>> pydna.Dseq(watson="agt",crick="actta",ovhg=-1)
    Dseq(-6)
    agt
     attca
    >>> pydna.Dseq(watson="agt",crick="actta",ovhg=0)
    Dseq(-5)
    agt
    attca
    >>> pydna.Dseq(watson="agt",crick="actta",ovhg=1)
    Dseq(-5)
     agt
    attca
    >>> pydna.Dseq(watson="agt",crick="actta",ovhg=2)
    Dseq(-5)
      agt
    attca

    If the ovhg parameter is psecified a crick strand also needs to be supplied,
    otherwise an exception is raised.

    >>> pydna.Dseq(watson="agt",ovhg=2)
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
      File "/usr/local/lib/python2.7/dist-packages/pydna_/dsdna.py", line 169, in __init__
        else:
    Exception: ovhg defined without crick strand!


    The default alphabet is set to Biopython IUPACAmbiguousDNA

    The shape of the fragment is set by either:

    1. linear   = False, True
    2. circular = True, False

    Note that both ends of the DNA fragment has to be compatible to set
    circular = True (or linear = False).


    >>> pydna.Dseq("aaa","ttt")
    Dseq(-3)
    aaa
    ttt
    >>> pydna.Dseq("aaa","ttt",ovhg=0)
    Dseq(-3)
    aaa
    ttt
    >>> pydna.Dseq("aaa", "ttt", linear = False ,ovhg=0)
    Dseq(o3)
    aaa
    ttt
    >>> pydna.Dseq("aaa", "ttt", circular = True , ovhg=0)
    Dseq(o3)
    aaa
    ttt

    Coercing to string

    >>> a=pydna.Dseq("tttcccc","aaacccc")
    >>> a
    Dseq(-11)
        tttcccc
    ccccaaa

    >>> str(a)
    'ggggtttcccc'

    The double stranded part is accessible through the dsdata property:

    >>> a.dsdata
    'ttt'

    A Dseq object can be longer that either the watson or crick strands.

    ::

        <-- length -->
        GATCCTTT
             AAAGCCTAG



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
    >>> import pydna
    >>> d=pydna.Dseq(s, linear=True)
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
    >>> d=pydna.Dseq(s, circular=True)
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
    >>> d=pydna.Dseq(s, circular=True)
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
    pydna.dsdna.Dseqrecord

    '''

    def __init__(self,
                  watson,
                  crick         = None,
                  ovhg          = None,
                  linear        = None,
                  circular      = None,
                  alphabet      = IUPACAmbiguousDNA() ):

        watson = "".join(watson.split())

        if ovhg is None:
            if crick is None:
                self.crick = rc(watson)
                self._ovhg = 0
            else:
                crick = "".join(crick.split())

                olaps = common_sub_strings(str(watson).lower(),
                                           str(rc(crick).lower()),
                                           int(math.log(len(watson))/math.log(4)))
                try:
                    F,T,L = olaps[0]
                except IndexError:
                    raise Exception("Could not anneal the two strands! "
                                    "ovhg should be provided")
                ovhgs = [ol[1]-ol[0] for ol in olaps if ol[2]==L]
                if len(ovhgs)>1:
                    for o in ovhgs:
                        print o
                    raise Exception("More than one way of annealing the strands "
                                    "ovhg should be provided")

                self._ovhg = T-F
                self.crick = crick
        elif crick is None:
            raise Exception("ovhg defined without crick strand!")
        else:
            self._ovhg = ovhg
            self.crick = "".join(crick.split())

        self.watson = watson

        sns = ((self._ovhg*" ")  + str(self.watson))
        asn = ((-self._ovhg*" ") + str(rc(self.crick)))

        self.todata = "".join([a.strip() or b.strip() for a,b in itertools.izip_longest(sns,asn, fillvalue=" ")])
        self.dsdata = "".join([a for a, b in itertools.izip_longest(sns,asn, fillvalue=" ") if a.lower()==b.lower()])

        if circular == None and linear in (True, False,):
            self._linear   = linear
            self._circular = not linear
        elif linear == None and circular in (True, False,):
            self._circular = circular
            self._linear   = not circular
        elif circular == linear == None:
            self._circular = False
            self._linear   = True
        elif linear in (True, False,) and circular in (True, False,) and circular != linear:
            self._circular = circular
            self._linear   = not circular
        else:
            raise Exception("circular and linear argument set to {} and {}, respectively\n".format(circular,linear)+
                            "circular and linear are each others opposites.")

        assert self._circular != self._linear

        if (self.circular and
            self.five_prime_end()[0]  != "blunt" and
            self.three_prime_end()[0] != "blunt"):
            raise Exception("DNA is circular, but has staggered ends!\n")

        Seq.__init__(self, self.todata, alphabet)

    def find(self, sub, start=0, end=sys.maxint):
        """Find method, like that of a python string.

        This behaves like the python string method of the same name.

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
        >>> import pydna
        >>> seq = Dseq("atcgactgacgtgtt")
        >>> seq
        Dseq(-15)
        atcgactgacgtgtt
        tagctgactgcacaa
        >>> seq.find("gac")
        3
        >>> seq = pydna.Dseq(watson="agt",crick="actta",ovhg=-2)
        >>> seq
        Dseq(-7)
        agt
          attca
        >>> seq.find("taa")
        2
        """

        if self.linear:
            return Seq.find(self, sub, start, end)

        sub_str = self._get_seq_str_and_check_alphabet(sub)

        return (str(self)+str(self)).find(sub_str, start, end)



    def __getitem__(self, sl):
        '''Returns a subsequence.
        '''

        if self.linear:
            sns = (self._ovhg*" " + self.watson)[sl]
            asn = (-self._ovhg*" " + self.crick[::-1])[sl]
            ovhg = max((len(sns) - len(sns.lstrip()),
                        -len(asn) + len(asn.lstrip())),
                        key=abs)
            return Dseq(sns.strip(), asn[::-1].strip(), ovhg=ovhg, linear=True)
        else:
            sl = slice(sl.start or 0,
                       sl.stop  or len(self),
                       sl.step)

            if sl.start<sl.stop:
                return Dseq(self.watson[sl],self.crick[::-1][sl][::-1], ovhg=0, linear=True)
            else:
                try:
                    stp = abs(sl.step)
                except TypeError:
                    stp = 1
                start=sl.start
                stop=sl.stop
                if not start:
                    start=0
                if not stop:
                    stop=len(self)

                w = self.watson[(start or len(self))::stp] + self.watson[:(stop or 0):stp]
                c = self.crick[len(self)-stop::stp] + self.crick[:len(self)-start:stp]

                return Dseq(w, c, ovhg=0, linear=True)

    def __eq__( self, other ):
        '''Compare to another Dseq object OR an object that implements
        watson, crick and ovhg properties. This comparison is case
        insensitive.

        '''
        try:
            same = (other.watson.lower() == self.watson.lower() and
                    other.crick.lower()  == self.crick.lower()  and
                    other.ovhg == self._ovhg)
        except AttributeError:
            same = False
        return same


    def fig(self):
        '''Returns a representation of the sequence, truncated if
       longer than 30 bp:

       Examples
       --------

       >>> import pydna
       >>> a=pydna.Dseq("atcgcttactagcgtactgatcatctgact")
       >>> a
       Dseq(-30)
       atcgcttactagcgtactgatcatctgact
       tagcgaatgatcgcatgactagtagactga
       >>> a+=Dseq("A")
       >>> a
       Dseq(-31)
       atcg..actA
       tagc..tgaT


       '''
        return self.__repr__()

    def __repr__(self):
        '''Returns a representation of the sequence, truncated if
        longer than 30 bp'''

        if len(self) > 30:

            if self.ovhg>0:
                d = self.crick[-self.ovhg:][::-1]
                hej = len(d)
                if len(d)>10:
                    d = "{}..{}".format(d[:4], d[-4:])
                a = len(d)*" "

            elif self.ovhg<0:
                a = self.watson[:max(0,-self.ovhg)]
                hej = len(a)
                if len(a)>10:
                    a = "{}..{}".format(a[:4], a[-4:])
                d = len(a)*" "
            else:
                a = ""
                d = ""
                hej=0

            x = self.ovhg+len(self.watson)-len(self.crick)

            if x>0:
                c=self.watson[len(self.crick)-self.ovhg:]
                y=len(c)
                if len(c)>10:
                    c = "{}..{}".format(c[:4], c[-4:])
                f=len(c)*" "
            elif x<0:
                f=self.crick[:-x][::-1]
                y=len(f)
                if len(f)>10:
                    f = "{}..{}".format(f[:4], f[-4:])
                c=len(f)*" "
            else:
                c = ""
                f = ""
                y=0

            L = len(self)-hej-y
            x1 = -min(0, self.ovhg)
            x2 = x1+L
            x3 = -min(0, x)
            x4 = x3+L

            b=self.watson[x1:x2]
            e=self.crick[x3:x4][::-1]

            if len(b)>10:
                b = "{}..{}".format(b[:4], b[-4:])
                e = "{}..{}".format(e[:4], e[-4:])

            #import sys;sys.exit()

            return ("{klass}({top}{size})\n"
                    "{a}{b}{c}\n"
                    "{d}{e}{f}").format(klass = self.__class__.__name__,
                                          top = {True:"-", False:"o"}[self.linear],
                                          size = len(self),
                                          a=a,
                                          b=b,
                                          c=c,
                                          d=d,
                                          e=e,
                                          f=f,)


        else:
            return "{}({}{})\n{}\n{}".format(self.__class__.__name__,
                                                {True:"-", False:"o"}[self.linear],
                                                len(self),
                                                self._ovhg*" " + self.watson,
                                               -self._ovhg*" "+ self.crick[::-1])

    def rc(self):
        '''Alias of the reverse_complement method'''
        return self.reverse_complement()

    def reverse_complement(self):
        '''Returns a Dseq object where watson and crick have switched
        places.

        Examples
        --------
        >>> import pydna
        >>> a=pydna.Dseq("catcgatc")
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

       '''
        ovhg = len(self.watson) - len(self.crick) + self._ovhg
        return Dseq(self.crick, self.watson, ovhg=ovhg, circular = self.circular)

    def looped(self):
        '''Returns a circularized Dseq object. This can only be done if the
        two ends are compatible, otherwise a TypeError is raised.

        Examples
        --------
        >>> import pydna
        >>> a=pydna.Dseq("catcgatc")
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

       '''
        if self.circular:
            return self
        type5, sticky5 = self.five_prime_end()
        type3, sticky3 = self.three_prime_end()
        if type5 == type3 and str(sticky5) == str(rc(sticky3)):
            nseq = Dseq(self.watson, self.crick[-self._ovhg:] + self.crick[:-self._ovhg], 0, circular=True)
            assert len(nseq.crick) == len(nseq.watson)
            return nseq
        else:
            raise TypeError("DNA cannot be circularized.\n"
                             "5' and 3' sticky ends not compatible!")

    def tolinear(self):
        '''Returns a blunt, linear copy of a circular Dseq object. This can
       only be done if the Dseq object is circular, otherwise a
       TypeError is raised.

       Examples
       --------

       >>> import pydna
       >>> a=pydna.Dseq("catcgatc", circular=True)
       >>> a
       Dseq(o8)
       catcgatc
       gtagctag
       >>> a.tolinear()
       Dseq(-8)
       catcgatc
       gtagctag
       >>>

       '''

        if self.linear:
            raise TypeError("DNA is not circular.\n")
        return self.__class__(self.watson, self.crick, ovhg=0, linear=True)

    def five_prime_end(self):
        '''Returns a tuple describing the structure of the 5' end of
        the DNA fragment

        Examples
        --------
        >>> import pydna
        >>> a=pydna.Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.five_prime_end()
        ('blunt', '')
        >>> a=pydna.Dseq("aaa", "ttt", ovhg=1)
        >>> a
        Dseq(-4)
         aaa
        ttt
        >>> a.five_prime_end()
        ("3'", 't')
        >>> a=pydna.Dseq("aaa", "ttt", ovhg=-1)
        >>> a
        Dseq(-4)
        aaa
         ttt
        >>> a.five_prime_end()
        ("5'", 'a')
        >>>

        See also
        --------
        pydna.dsdna.Dseq.three_prime_end

        '''
        if self.watson and not self.crick:
            return "5'",self.watson.lower()
        if not self.watson and self.crick:
            return "3'",self.crick.lower()
        if self._ovhg < 0:
            sticky = self.watson[:-self._ovhg].lower()
            type_ = "5'"
        elif self._ovhg > 0:
            sticky = self.crick[-self._ovhg:].lower()
            type_ = "3'"
        else:
            sticky = ""
            type_ = "blunt"
        return type_, sticky

    def three_prime_end(self):
        '''Returns a tuple describing the structure of the 5' end of
        the DNA fragment

        >>> import pydna
        >>> a=pydna.Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.three_prime_end()
        ('blunt', '')
        >>> a=pydna.Dseq("aaa", "ttt", ovhg=1)
        >>> a
        Dseq(-4)
         aaa
        ttt
        >>> a.three_prime_end()
        ("3'", 'a')
        >>> a=pydna.Dseq("aaa", "ttt", ovhg=-1)
        >>> a
        Dseq(-4)
        aaa
         ttt
        >>> a.three_prime_end()
        ("5'", 't')
        >>>

        See also
        --------
        pydna.dsdna.Dseq.five_prime_end

        '''

        ovhg = len(self.watson)-len(self.crick)+self._ovhg

        if ovhg < 0:
            sticky = self.crick[:-ovhg].lower()
            type_ = "5'"
        elif ovhg > 0:
            sticky = self.watson[-ovhg:].lower()
            type_ = "3'"
        else:
            sticky = ''
            type_ = "blunt"
        return type_, sticky

    def __add__(self, other):
        '''Simulates ligation between two DNA fragments.

        Add other Dseq object at the end of the sequence.
        Type error if all of the points below are fulfilled:

        * either objects are circular
        * if three prime sticky end of self is not the same type
          (5' or 3') as the sticky end of other
        * three prime sticky end of self complementary with five
          prime sticky end of other.

        Phosphorylation and dephosphorylation is not considered.
        DNA is allways presumed to have the necessary 5' phospate
        group necessary for ligation.

       '''
        # test for circular DNA
        if self.circular:
            raise TypeError("circular DNA cannot be ligated!")
        try:
            if other.circular:
                raise TypeError("circular DNA cannot be ligated!")
        except AttributeError:
            pass

        self_type,  self_tail  = self.three_prime_end()
        other_type, other_tail = other.five_prime_end()

        if (self_type == other_type and
            str(self_tail) == str(rc(other_tail))):
            answer = Dseq(self.watson + other.watson,
                          other.crick + self.crick,
                          self._ovhg,)
        elif not self:
            answer = copy.copy(other)
        elif not other:
            answer = copy.copy(self)
        else:
            raise TypeError("sticky ends not compatible!")
        return answer

    def __mul__(self, number):
        if not isinstance(number, int):
            raise TypeError("TypeError: can't multiply Dseq by non-int of type {}".format(type(number)))
        if number<=0:
            return self.__class__("")
        new = copy.copy(self)
        for i in range(number-1):
            new += self
        return new

    def _fill_in_five_prime(self, nucleotides):
        stuffer = ''
        type, se = self.five_prime_end()
        if type == "5'":
            for n in rc(se):
                if n in nucleotides:
                    stuffer+=n
                else:
                    break
        return self.crick+stuffer, self._ovhg+len(stuffer)

    def _fill_in_three_prime(self, nucleotides):
        stuffer = ''
        type, se = self.three_prime_end()
        if type == "5'":
            for n in rc(se):
                if n in nucleotides:
                    stuffer+=n
                else:
                    break
        return self.watson+stuffer

    def fill_in(self, nucleotides=None):
        '''Fill in of five prime protruding end with a DNA polymerase
        that has only DNA polymerase activity (such as exo-klenow [#]_)
        and any combination of A, G, C or T. Default are all four
        nucleotides together.

        Parameters
        ----------

        nucleotides : str

        Examples
        --------

        >>> import pydna
        >>> a=pydna.Dseq("aaa", "ttt")
        >>> a
        Dseq(-3)
        aaa
        ttt
        >>> a.fill_in()
        Dseq(-3)
        aaa
        ttt
        >>> b=pydna.Dseq("caaa", "cttt")
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
        >>> c=pydna.Dseq("aaac", "tttg")
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

        '''
        if not nucleotides:
            nucleotides = self.alphabet.letters
        nucleotides = set(nucleotides.lower()+nucleotides.upper())
        crick, ovhg = self._fill_in_five_prime(nucleotides)
        watson = self._fill_in_three_prime(nucleotides)
        return Dseq(watson, crick, ovhg)

    def mung(self):
        '''
       Simulates treatment a nuclease with 5'-3' and 3'-5' single
       strand specific exonuclease activity (such as mung bean nuclease [#]_)

       ::

            ggatcc    ->     gatcc
             ctaggg          ctagg

             ggatcc   ->      ggatc
            tcctag            cctag

        >>> import pydna
        >>> b=pydna.Dseq("caaa", "cttt")
        >>> b
        Dseq(-5)
        caaa
         tttc
        >>> b.mung()
        Dseq(-3)
        aaa
        ttt
        >>> c=pydna.Dseq("aaac", "tttg")
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


        '''
        return Dseq(self.dsdata)

    def t4(self,*args,**kwargs):
        '''Alias for the :func:`T4` method '''
        return self.T4(*args,**kwargs)

    def T4(self, nucleotides=None):
        '''Fill in of five prime protruding ends and chewing back of
       three prime protruding ends by a DNA polymerase providing both
       5'-3' DNA polymerase activity and 3'-5' nuclease acitivty
       (such as T4 DNA polymerase). This can be done in presence of any
       combination of the four A, G, C or T. Default are all four
       nucleotides together.

       Alias for the :func:`t4` method

       Parameters
       ----------

       nucleotides : str


       Examples
       --------

       >>> import pydna
       >>> a=pydna.Dseq("gatcgatc")
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

       '''

        if not nucleotides: nucleotides = self.alphabet.letters
        nucleotides = set(nucleotides.lower() + nucleotides.upper())
        type, se = self.five_prime_end()
        crick = self.crick
        if type == "5'":
            crick, ovhg = self._fill_in_five_prime(nucleotides)
        else:
            if type == "3'":
                ovhg = 0
                crick = self.crick[:-len(se)]
        x = len(crick)-1
        while x>=0:
            if crick[x] in nucleotides:
                break
            x-=1
        ovhg = x-len(crick)+1
        crick = crick[:x+1]
        if not crick: ovhg=0
        watson = self.watson
        type, se = self.three_prime_end()
        if type == "5'":
            watson = self._fill_in_three_prime(nucleotides)
        else:
            if type == "3'":
                watson = self.watson[:-len(se)]
        x = len(watson)-1
        while x>=0:
            if watson[x] in nucleotides:
                break
            x-=1
        watson=watson[:x+1]
        return Dseq(watson, crick, ovhg)

    def _cut(self, *enzymes):

        output = []
        stack = []
        stack.extend(reversed(enzymes))
        while stack:
            top = stack.pop()
            if hasattr(top, "__iter__"):
                stack.extend(reversed(top))
            else:
                output.append(top)
        enzymes = output

        if not hasattr(enzymes, '__iter__'):
            enzymes = (enzymes,)

        if self.circular:
            frags=[self.tolinear()*3,]
        else:
            frags=[self,]

        newfrags=[]

        enzymes = [e for (p,e) in sorted([(enzyme.search(Seq(frags[0].dsdata))[::-1], enzyme) for enzyme in enzymes], reverse=True) if p]

        if not enzymes:
            return [self,]

        for enzyme in enzymes:
            for frag in frags:

                if enzyme.search(Seq(frag.dsdata)):

                    watson_fragments = [str(s) for s in enzyme.catalyze(Seq(frag.watson+"N"))]
                    crick_fragments  = [str(s) for s in enzyme.catalyze(Seq(frag.crick+"N" ))[::-1]]

                    watson_fragments[-1] = watson_fragments[-1][:-1]
                    crick_fragments[0]   = crick_fragments[0][:-1]

                    s = zip(watson_fragments, crick_fragments)

                    if frag.linear:
                        newfrags.append(Dseq(*s.pop(0),
                                             ovhg = frag.ovhg,
                                             linear = True))
                        for seqs in s:
                            newfrags.append(Dseq(*seqs,
                                                 ovhg = enzyme.ovhg,
                                                 linear = True))
                    else:
                        for seqs in s:
                            newfrags.append(Dseq(*seqs,
                                                 ovhg=enzyme.ovhg,
                                                 linear=True))
                else:
                    newfrags.append(frag)
            frags=newfrags
            newfrags=[]

        if self.circular:
            swl = len(self.watson)
            frags = frags[1:-1]
            newfrags = [frags.pop(0),]
            while sum(len(f.watson) for f in newfrags) < swl:
                newfrags.append(frags.pop(0))
            frags = newfrags

        return frags

    def cut(self, *enzymes):
        '''Returns a list of linear Dseq fragments produced in the digestion.
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

        >>> from pydna import Dseq
        >>> seq=Dseq("ggatccnnngaattc")
        >>> seq
        Dseq(-15)
        ggatccnnngaattc
        cctaggnnncttaag
        >>> from Bio.Restriction import BamHI,EcoRI
        >>> type(seq.cut(BamHI))
        <type 'list'>
        >>> for frag in seq.cut(BamHI): print(frag.fig())
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

        '''

        output = []
        stack = []
        stack.extend(reversed(enzymes))
        while stack:
            top = stack.pop()
            if hasattr(top, "__iter__"):
                stack.extend(reversed(top))
            else:
                output.append(top)
        enzymes = output
        if not hasattr(enzymes, '__iter__'):
            enzymes = (enzymes,)

        if self.circular:
            frags=[self.tolinear()*3,]
        else:
            frags=[self,]

        newfrags=[]

        enzymes = [e for (p,e) in sorted([(enzyme.search(Seq(frags[0].dsdata))[::-1], enzyme) for enzyme in enzymes], reverse=True) if p]


        if not enzymes:
            return []

        for enz in enzymes:
            for frag in frags:

                ws = [x-1 for x in enz.search(Seq(frag.watson)+"N")] #, linear = frag.linear
                cs = [x-1 for x in enz.search(Seq(frag.crick) +"N")] #, linear = frag.linear

                sitepairs = [(sw, sc) for sw, sc in zip(ws,cs[::-1])
                             if (sw + max(0, frag.ovhg) -
                             max(0, enz.ovhg)
                             ==
                             len(frag.crick)-sc -
                             min(0, frag.ovhg) +
                             min(0, enz.ovhg))]

                sitepairs = sitepairs + [(len(frag.watson), 0)]

                w2, c1 = sitepairs[0]

                nwat = frag.watson[:w2]
                ncrk = frag.crick[c1:]

                newfrags.append(Dseq(nwat, ncrk, ovhg=frag.ovhg))

                for (w1, c2), (w2, c1)  in zip(sitepairs[:-1], sitepairs[1:]):
                    nwat = frag.watson[w1:w2]
                    ncrk = frag.crick[c1:c2]
                    newfrag = Dseq(nwat, ncrk, ovhg = enz.ovhg)
                    newfrags.append(newfrag)

                if not newfrags:
                    newfrags.append(frag)

            frags=newfrags
            newfrags=[]

        if self.circular:
            swl = len(self.watson)
            frags = frags[1:-1]
            newfrags = [frags.pop(0),]
            while sum(len(f.watson) for f in newfrags) < swl:
                newfrags.append(frags.pop(0))
            frags = newfrags[-1:] + newfrags[:-1]
        return frags

    def seguid(self):
        rc_ovhg = len(self.watson) - len(self.crick) + self._ovhg
        if self.ovhg<rc_ovhg:
            w = self.watson
            c = self.crick
            o  =self.ovhg
        elif self.ovhg>rc_ovhg:
            w = self.crick
            c = self.watson
            o = rc_ovhg
        elif self.ovhg==rc_ovhg:
            w, c = sorted((self.watson, self.crick))
            o = self.ovhg
        return seg( str(o) + w + "|" + c)

    @property
    def ovhg(self):
        '''The ovhg property'''
        return self._ovhg

    @property
    def linear(self):
        '''The linear property'''
        return self._linear

    @property
    def circular(self):
        '''The circular property'''
        return self._circular


class Dseqrecord(SeqRecord):
    '''Dseqrecord is a double stranded version of the Biopython SeqRecord [#]_ class.
    The Dseqrecord object holds a Dseq object describing the sequence.
    Additionally, Dseqrecord hold meta information about the sequence in the
    from of a list of SeqFeatures, in the same way as the SeqRecord does.
    The Dseqrecord can be initialized with a string, Seq, Dseq, SeqRecord
    or another Dseqrecord. The sequence information will be stored in a
    Dseq object in all cases. Dseqrecord objects can be read or parsed
    from sequences in Fasta, Embl or Genbank format.

    There is a short representation associated with the Dseqrecord.
    ``Dseqrecord(-3)`` represents a linear sequence of length 2
    while ``Dseqrecord(o7)``
    represents a circular sequence of length 7.

    Dseqrecord and Dseq share the same concept of length

    ::

        <-- length -->
        GATCCTTT
             AAAGCCTAG




    Parameters
    ----------
    record  : string, Seq, SeqRecord, Dseq or other Dseqrecord object
        This data will be used to form the seq property

    circular : bool, optional
        True or False reflecting the shape of the DNA molecule

    linear : bool, optional
        True or False reflecting the shape of the DNA molecule


    Examples
    --------

    >>> from pydna import Dseqrecord
    >>> a=Dseqrecord("aaa")
    >>> a
    Dseqrecord(-3)
    >>> a.seq
    Dseq(-3)
    aaa
    ttt
    >>> from Bio.Seq import Seq
    >>> b=Dseqrecord(Seq("aaa"))
    >>> b
    Dseqrecord(-3)
    >>> b.seq
    Dseq(-3)
    aaa
    ttt
    >>> from Bio.SeqRecord import SeqRecord
    >>> c=Dseqrecord(SeqRecord(Seq("aaa")))
    >>> c
    Dseqrecord(-3)
    >>> c.seq
    Dseq(-3)
    aaa
    ttt
    >>> a.seq.alphabet
    IUPACAmbiguousDNA()
    >>> b.seq.alphabet
    IUPACAmbiguousDNA()
    >>> c.seq.alphabet
    IUPACAmbiguousDNA()
    >>>

    References
    ----------

    .. [#] http://biopython.org/wiki/SeqRecord

    '''

    def __init__(self, record,
                       circular               = None,
                       linear                 = None,
                       n                      = 10E-12, # pmols
                       *args,
                       **kwargs):
        self.n = n
        if circular == None and linear in (True, False,):
            circular = not linear
        elif linear == None and circular in (True, False,):
            linear   = not circular

        try:
            record.letter_annotations = {}
        except AttributeError:
            pass

        if isinstance(record, basestring):  # record is a string
            SeqRecord.__init__(self,
                               Dseq(record,
                                    rc(record),
                                    ovhg=0 ,
                                    linear=linear,
                                    circular=circular),
                               *args,
                               **kwargs)
        elif hasattr(record, "features"): # record is SeqRecord or Dseqrecord?
            for key, value in record.__dict__.items():
                setattr(self, key, value )
            if hasattr(record.seq, "watson"): # record.seq is a Dseq, so record is Dseqrecord
                new_seq = copy.copy(record.seq)
                if new_seq.circular and linear:
                    new_seq = new_seq.tolinear()
                if new_seq.linear and circular:
                    new_seq = new_seq.looped()
                self.seq=new_seq
            else:                             # record is Bio.SeqRecord
                self.seq=Dseq(str(self.seq),
                              rc(str(self.seq)),
                              ovhg=0 ,
                              linear=linear,
                              circular=circular)
        elif hasattr(record, "watson"):      # record is Dseq ?
            if record.circular and linear:
                record = record.tolinear()
            if record.linear and circular:
                record = record.looped()
            SeqRecord.__init__(self, record, *args, **kwargs)
        elif isinstance(record, Seq):         # record is Bio.Seq ?
            SeqRecord.__init__(self, Dseq(str(record),
                                          str(record.reverse_complement()),
                                          ovhg=0 ,
                                          linear=linear,
                                          circular=circular),
                                          *args,
                                          **kwargs)
        else:
            raise TypeError(("record argument needs to be a string,"
                              "Seq, SeqRecord, Dseq or Dseqrecord object,"
                              " got {}").format(type(record)))

        if len(self.name)>16:
            short_name = self.name[:16]
            warn("name property {} truncated to 16 chars {}".format(self.name, short_name))
            self.name = short_name

        if self.name == "<unknown name>":
            self.name = "na"

        if self.id == "<unknown id>":
            self.id = "-"

        if self.description =="<unknown description>":
            self.description = "@"

        if not 'date' in self.annotations:
            self.annotations.update({"date": datetime.date.today().strftime("%d-%b-%Y").upper()})

        self.map_target = None

    @property
    def linear(self):
        '''The linear property'''
        return self.seq.linear

    @property
    def circular(self):
        '''The circular property'''
        return self.seq.circular



    @property
    def locus(self):
        ''' alias for name property '''
        return self.name

    @locus.setter
    def locus(self, value):
        ''' alias for name property '''
        if len(value)>16:
            raise Exception()
        self.name = value
        return



    @property
    def accession(self):
        ''' alias for id property '''
        return self.id

    @accession.setter
    def accession(self, value):
        ''' alias for id property '''
        self.id = value
        return

    @property
    def definition(self):
        ''' alias for description property '''
        return self.description

    @definition.setter
    def definition(self, value):
        ''' alias for id property '''
        self.description = value
        return

    def seguid(self):
        '''Returns the url safe SEGUID [#]_ for the sequence.
           This checksum is the same as seguid but with base64.urlsafe
           encoding [#]_ instead of the normal base 64. This means that
           the characters + and / are replaced with - and _ so that
           the checksum can be a pert of and URL or a filename.

           Examples
           --------
           >>> import pydna
           >>> a=pydna.Dseqrecord("aaaaaaa")
           >>> a.seguid() # original seguid is +bKGnebMkia5kNg/gF7IORXMnIU
           '-bKGnebMkia5kNg_gF7IORXMnIU'

           References
           ----------

       .. [#] http://wiki.christophchamp.com/index.php/SEGUID

       '''
        return seg(self.seq)

    def isorf(self, table=1):
        '''Detects if sequence is an open reading frame (orf) in the 5'-3' direction.
        Translation tables are numbers according to the NCBI numbering [#]_.

        Parameters
        ----------
        table  : int
            Sets the translation table, default is 1 (standard code)

        Returns
        -------
        bool
            True if sequence is an orf, False otherwise.


        Examples
        --------

        >>> from pydna import Dseqrecord
        >>> a=Dseqrecord("atgtaa")
        >>> a.isorf()
        True
        >>> b=Dseqrecord("atgaaa")
        >>> b.isorf()
        False
        >>> c=Dseqrecord("atttaa")
        >>> c.isorf()
        False

        References
        ----------

        .. [#] http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c

        '''

        try:
            self.seq.translate(table=table, cds=True)
        except TranslationError:
            return False
        else:
            return True

    def add_feature(self, x=None, y=None, seq=None, label=None, type="misc", **kwargs):
        '''Adds a feature of type misc to the feature list of the sequence.

        Parameters
        ----------
        x  : int
            Indicates start of the feature
        y  : int
            Indicates end of the feature

        Examples
        --------

        >>> from pydna import Dseqrecord
        >>> a=Dseqrecord("atgtaa")
        >>> a.features
        []
        >>> a.add_feature(2,4)
        >>> a.features
        [SeqFeature(FeatureLocation(ExactPosition(2), ExactPosition(4)), type='misc')]
        '''
        qualifiers = {"label": label}
        if seq:
            seq = Dseqrecord(seq)
            x = self.seq.lower().find(seq.seq.lower())
            if x==-1:
                return
            y = x + len(seq)
        self.features.append( SeqFeature(FeatureLocation(x, y), type=type, qualifiers = qualifiers, **kwargs))

        '''

        In [11]: a.seq.translate()
        Out[11]: Seq('K', ExtendedIUPACProtein())

        In [12]:
        '''


    def extract_feature(self, n):
        '''Extracts a feature and creates a new Dseqrecord object.

        Parameters
        ----------
        n  : int
            Indicates the feature to extract

        Examples
        --------

        >>> from pydna import Dseqrecord
        >>> a=Dseqrecord("atgtaa")
        >>> a.add_feature(2,4)
        >>> b=a.extract_feature(0)
        >>> b
        Dseqrecord(-2)
        >>> b.seq
        Dseq(-2)
        gt
        ca

        '''
        return self.features[n].extract(self)

    def spread_ape_colors(self):
    	''' This method assigns random colors compatible with the ApE editor
    	to features.
    	'''

        def get_N_HexCol(N):
            HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
            hex_out = []
            for rgb in HSV_tuples:
                rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
                hex_out.append("".join(map(lambda x: chr(x).encode('hex'),rgb)))
            return hex_out

        for i, color in enumerate(get_N_HexCol(len(self.features))):
            self.features[i].qualifiers['ApEinfo_fwdcolor'] = "#"+color
            self.features[i].qualifiers['ApEinfo_revcolor'] = "#"+color

    def olaps(self, other, *args, **kwargs):
        other = Dseqrecord(other)
        olaps = common_sub_strings(str(self.seq).lower(), str(other.seq).lower(), **kwargs)
        return [ self[olap[0]:olap[0]+olap[2]] for olap in olaps ]

    def list_features(self):
        '''Prints an ASCII table with all features.

        Examples
        --------

        >>> from pydna import Dseqrecord
        >>> a=Dseqrecord("atgtaa")
        >>> a.add_feature(2,4)
        >>> print(a.list_features())
        +----------+-----------+-------+-----+--------+--------------+------+------+
        | Feature# | Direction | Start | End | Length | id           | type | orf? |
        +----------+-----------+-------+-----+--------+--------------+------+------+
        | 0        |    None   |   2   |  4  |      2 | <unknown id> | misc |  no  |
        +----------+-----------+-------+-----+--------+--------------+------+------+
        >>>
        '''

        x = PrettyTable(["Feature#", "Direction", "Start", "End", "Length", "id", "type", "orf?"])
        x.align["Feature#"] = "l" # Left align
        x.align["Length"] = "r"
        x.align["id"] = "l"
        x.align["type"] = "l"
        x.padding_width = 1 # One space between column edges and contents
        for i, sf in enumerate(self.features):
            x.add_row([ i,
                        {1:"-->", -1:"<--", 0:"---", None:"None"}[sf.strand],
                        sf.location.start,
                        sf.location.end,
                        len(sf), sf.id, sf.type,
                        {True:"yes",False:"no"}[self.extract_feature(i).isorf() or self.extract_feature(i).rc().isorf()]])
        return pretty_str(x)

    def gc(self):
    	'''Returns GC content '''
        return pretty_string(round(GC(str(self.seq)), 1))

    def cseguid(self):
        '''Returns the url safe cSEGUID for the sequence.

        Only defined for circular sequences.

        The cSEGUID checksum uniqely identifies a circular
        sequence regardless of where the origin is set.
        The two Dseqrecord objects below are circular
        permutations.

        Examples
        --------

        >>> import pydna
        >>> a=pydna.Dseqrecord("agtatcgtacatg", circular=True)
        >>> a.cseguid() # cseguid is CTJbs6Fat8kLQxHj+/SC0kGEiYs
        'CTJbs6Fat8kLQxHj-_SC0kGEiYs'

        >>> a=pydna.Dseqrecord("gagtatcgtacat", circular=True)
        >>> a.cseguid()
        'CTJbs6Fat8kLQxHj-_SC0kGEiYs'

       '''
        if self.linear:
            raise Exception("cseguid is only defined for circular sequences.")
        return cseg(self.seq)

    def lseguid(self):
        '''Returns the url safe lSEGUID for the sequence.

        Only defined for linear double stranded sequences.

        The lSEGUID checksum uniqely identifies a linear
        sequence independent of the direction.
        The two Dseqrecord objects below are each others
        reverse complements, so they do in fact refer to
        the same molecule.

        Examples
        --------

        >>> import pydna
        >>> a=pydna.Dseqrecord("agtatcgtacatg")
        >>> a.lseguid() # cseguid is CTJbs6Fat8kLQxHj+/SC0kGEiYs
        'QFuP3noYs92MGFJ2YGymCrxXFU4'

        >>> b=pydna.Dseqrecord("catgtacgatact")
        >>> a.lseguid()
        'QFuP3noYs92MGFJ2YGymCrxXFU4'

       '''
        if self.circular:
            raise Exception("lseguid is only defined for linear sequences.")
        return self.seq.seguid()

    def stamp(self, chksum = (("SEGUID", seg),("cSEGUID", cseg))):
        '''Adds a checksum to the description property. This will
        show in the genbank format. Default is seguid for linear sequences
        and cseguid for circular.

        The following string:

        ``<type><seguid>``

        For example:

        ``SEGUID_<seguid>``

        for linear sequences or:

        ``cSEGUID_<seguid>``

        for circular sequences will be appended to the description property
        of the Dseqrecord object (string).

        https://xkcd.com/1179/

        The stamp can be verified with :func:`verify_stamp`

        Examples
        --------

        >>> import pydna
        >>> a=pydna.Dseqrecord("aaa")
        >>> a.stamp()
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE...'
        >>> a.description
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE...'
        >>> a.verify_stamp()
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE'

        See also
        --------
        pydna.dsdna.Dseqrecord.verify_stamp
        '''

        name, alg = {True:chksum[0], False:chksum[1]}[self.linear]

        now = datetime.datetime.utcnow().isoformat("T")

        pattern = "({name})_\s*\S{{27}}_".format(name=name)

        stamp = re.search(pattern, self.description)

        if not stamp:
            stamp = "{}_{}_{}".format(name,
                                      alg(str(self.seq)),
                                      now)
            if not self.description or self.description=="@":
                self.description = stamp
            elif not re.search(pattern, self.description):
                self.description += " "+stamp
        else:
            raise Exception("sequence already stamped {}")

        return pretty_str(stamp)

    def verify_stamp(self, chksum = (("SEGUID", seg),("cSEGUID", cseg))):
        '''Verifies the SEGUID stamp in the description property is
       valid. True if stamp match the sequid calculated from the sequence.
       Exception raised if no stamp can be found.

        >>> import pydna
        >>> b=pydna.read(">a\\naaa")
        >>> b.annotations['date'] = '02-FEB-2013'
        >>> b.seguid()
        'YG7G6b2Kj_KtFOX63j8mRHHoIlE'
        >>> print(b.format("gb"))
        LOCUS       a                          3 bp    DNA     linear   UNK 02-FEB-2013
        DEFINITION  a
        ACCESSION   a
        VERSION     a
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 aaa
        //
        >>> b.stamp()
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE_...'
        >>> b
        Dseqrecord(-3)
        >>> print(b.format("gb"))
        LOCUS       a                          3 bp    DNA     linear   UNK 02-FEB-2013
        DEFINITION  a SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE_...
        ACCESSION   a
        VERSION     a
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 aaa
        //
        >>> b.verify_stamp()
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE'
        >>>

       See also
       --------
       pydna.dsdna.Dseqrecord.stamp

       '''
        name, alg = {True:chksum[0], False:chksum[1]}[self.linear]
        pattern = "{name}_{chksum}".format(name=name, chksum=alg(self.seq))

        if not pattern in self.description:
            raise Exception("No stamp present in the description property.")
        else:
            return pretty_str(pattern)



    def looped(self):
        '''
        Returns a circular version of the Dseqrecord object. The
        underlying Dseq object has to have compatible ends.


        Examples
        --------
        >>> import pydna
        >>> a=pydna.Dseqrecord("aaa")
        >>> a
        Dseqrecord(-3)
        >>> b=a.looped()
        >>> b
        Dseqrecord(o3)
        >>>

        See also
        --------
        pydna.dsdna.Dseq.looped
        '''
        new = copy.copy(self)
        for key, value in self.__dict__.items():
            setattr(new, key, value )
        new._seq = self.seq.looped()
        for fn, fo in zip(new.features, self.features):
            fn.qualifiers = fo.qualifiers

        return new

    def tolinear(self):
        '''
        Returns a linear, blunt copy of a circular Dseqrecord object. The
        underlying Dseq object has to be circular.

        Examples
        --------
        >>> import pydna
        >>> a=pydna.Dseqrecord("aaa", circular = True)
        >>> a
        Dseqrecord(o3)
        >>> b=a.tolinear()
        >>> b
        Dseqrecord(-3)
        >>>

        '''

        new = copy.copy(self)
        for key, value in self.__dict__.items():
            setattr(new, key, value )
        new._seq = self.seq.tolinear()
        for fn, fo in zip(new.features, self.features):
            fn.qualifiers = fo.qualifiers

        return new

    def format(self, f="gb"):
        '''Returns the sequence as a string using a format supported by Biopython
        SeqIO [#]_. Default is "gb" which is short for Genbank.

        Examples
        --------
        >>> import pydna
        >>> x=pydna.Dseqrecord("aaa")
        >>> x.annotations['date'] = '02-FEB-2013'
        >>> x
        Dseqrecord(-3)
        >>> print(x.format("gb"))
        LOCUS       na                         3 bp    DNA     linear   UNK 02-FEB-2013
        DEFINITION  @
        ACCESSION   -
        VERSION     -
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 aaa
        //


        References
        ----------

        .. [#] http://biopython.org/wiki/SeqIO


        '''

        s = SeqRecord.format(self, f).strip()

        if f in ("genbank","gb"):
            if self.circular:
                return pretty_string(s[:55]+"circular"+s[63:])
            else:
                return pretty_string(s[:55]+"linear"+s[61:])
        else:
            return pretty_string(s)

    def write(self, filename=None, f="gb"):
        '''Writes the Dseqrecord to a file using the format f, which must
        be a format supported by Biopython SeqIO for writing [#]_. Default
        is "gb" which is short for Genbank. Note that Biopython SeqIO reads
        more formats than it writes.

        Filename is the path to the file where the sequece is to be
        written. The filename is optional, if it is not given, the
        description property (string) is used together with the format.

        If obj is the Dseqrecord object, the default file name will be:

        ``<obj.description>.<f>``

        Where <f> is "gb" by default. If the filename already exists and
        AND the sequence it contains is different, a new file name will be
        used so that the old file is not lost:

        ``<obj.description>_NEW.<f>``

        References
        ----------
        .. [#] http://biopython.org/wiki/SeqIO

       '''
        if not filename:
            filename = "{name}.{type}".format(name=self.description, type=f)
            # invent a name if none given
        if isinstance(filename, basestring):
            name, ext = os.path.splitext(filename)
            result = "###[{name}]({filename})".format(name=name, filename=filename)
            if os.path.isfile(filename):
                len_new    = len(self)
                seguid_new = self.seguid()
                old_file   = read(filename)
                len_old    = len(old_file)
                seguid_old = old_file.seguid()
                if seguid_new != seguid_old or self.circular != old_file.circular:
                    # If new sequence is different, the old file is saved with "OLD" suffix
                    old_filename = "{}_OLD{}".format(name, ext)
                    os.rename(filename, old_filename)
                    result = ('#Sequence changed!\n'
                              '[{old_filename}]({old_filename}) {len_old} bp seguid {seguid_old}\n\n'
                              '[{filename}]({filename}) {len_new} bp seguid {seguid_new}\n').format(old_filename=old_filename,
                                                                                      len_old=len_old,
                                                                                      seguid_old=seguid_old,
                                                                                      filename=filename,
                                                                                      len_new=len_new,
                                                                                      seguid_new=seguid_new,
                                                                                      w=max(len(filename),len(old_filename)),
                                                                                      wn=len(str(max(len_new, len_old))))
            with open(filename, "w") as fp: fp.write(self.format(f))

        else:
            raise Exception("filename has to be a string, got", type(filename))
        return display(result)


    def __str__(self):
        return ( "Dseqrecord\n"
                 "circular: {}\n"
                 "size: {}\n").format(self.circular, len(self))+SeqRecord.__str__(self)

    def __contains__(self, other):
        if other.lower() in str(self.seq).lower():
            return True
        else:
            s = self.seq.watson.replace(" ","")
            ln  = len(s)
            spc = 3-ln%3 if ln%3 else 0
            s = "n" * spc + s + "nnn"
            for frame in range(3):
                if other.lower() in translate(s[frame:frame+spc+ln]).lower():
                    return True

            #s = self.seq.crick.replace(" ","")
            #ln=len(s)
            #spc = 3-ln%3 if ln%3 else 0
            #s = "n" * spc + s + "nnn"
            #for frame in range(3):
            #    if other.lower() in translate(s[frame:frame+spc+ln]).lower():
            #        return True
        return False

    def find_aa(self, other):
        return self.find_aminoacids(other)

    def find_aminoacids(self, other):
        '''
        >>> import pydna
        >>> s=pydna.Dseqrecord("atgtacgatcgtatgctggttatattttag")
        >>> s.seq.translate()
        Seq('MYDRMLVIF*', HasStopCodon(ExtendedIUPACProtein(), '*'))
        >>> "RML" in s
        True
        >>> "MMM" in s
        False
        >>> s.seq.rc().translate()
        Seq('LKYNQHTIVH', ExtendedIUPACProtein())
        >>> "QHT" in s.rc()
        True
        >>> "QHT" in s
        False
        >>> slc = s.find_aa("RML")
        >>> slc
        slice(9, 18, None)
        >>> s[slc]
        Dseqrecord(-9)
        >>> code = s[slc].seq
        >>> code
        Dseq(-9)
        cgtatgctg
        gcatacgac
        >>> code.translate()
        Seq('RML', ExtendedIUPACProtein())
        '''
        other = str(other).lower()
        assert self.seq.watson == "".join(self.seq.watson.split())
        s = self.seq.watson
        ln  = len(s)
        spc = 3-ln%3 if ln%3 else 0
        s = s + "n"*spc + "nnn"
        start = None
        for frame in range(3):
            try:
                start = translate(s[frame:frame+ln+spc]).lower().index(other)
                break
            except ValueError:
                pass
        oh = self.seq.ovhg if self.seq.ovhg>0 else 0
        if start == None:
            return None
        else:
            return slice(frame+start*3+oh, frame+(start+len(other))*3+oh)

    def map_trace_files(self, pth):
        import glob
        traces = []
        stf = SequenceTraceFactory()
        for name in glob.glob(pth):
            traces.append( stf.loadTraceFile( name ))
        if not traces:
            raise(Exception("no trace files found!"))
        if hasattr( self.map_target, "step" ):
            area = self.map_target
        elif hasattr( self.map_target, "extract" ):
            area = slice(self.map_target.location.start, self.map_target.location.end)
        else:
            area = None

        if area:
            self.matching_reads     = []
            self.not_matching_reads = []
            target    = str(self[area].seq).lower()
            target_rc = str(self[area].seq.rc()).lower()
            for trace in traces:
                if target in str(trace.basecalls).lower() or target_rc in str(trace.basecalls).lower():
                   self.matching_reads.append(trace)
                else:
                    self.not_matching_reads.append(trace)
            reads = self.matching_reads
        else:
            self.matching_reads     = None
            self.not_matching_reads = None
            reads = traces

        for read in reads:
            matches = common_sub_strings(str(self.seq).lower(), read.basecalls.lower(), 25)
            if len(matches)>1:
                newmatches = [matches[0],]
                for i, x in enumerate(matches[1:]):
                    g,f,h = matches[i]
                    if g+h < x[0] and f+h < x[1]:
                        newmatches.append(x)
            elif len(matches)==1:
                newmatches = matches
            else:
                continue

            if len(newmatches)>1:
                ms = []
                for m in newmatches:
                    ms.append(FeatureLocation(m[0], m[0]+m[2]))
                loc = CompoundLocation(ms)
            else:
                a,b,c = newmatches[0]
                loc = FeatureLocation(a,a+c)

            self.features.append( SeqFeature(loc,
                                             qualifiers = {"label": read.getFileName()},
                                             type="trace") )

        return [x.getFileName() for x in reads]

    def __repr__(self):
        return "Dseqrecord({}{})".format({True:"-", False:"o"}[self.linear],len(self))

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text('Dseqrecord(...)')
        else:
            p.text("Dseqrecord({}{})".format({True:"-", False:"o"}[self.linear],len(self)))

    def __add__(self, other):
        if hasattr(other, "seq") and hasattr(other.seq, "watson"):
            offset = other.seq.ovhg
            other.features = [f._shift(offset) for f in other.features]
            #answer = self.__class__(SeqRecord.__add__(self, other))
            answer = Dseqrecord(SeqRecord.__add__(self, other))
            answer.n = min(self.n, other.n)
        else:
            #answer = self.__class__(SeqRecord.__add__(self, Dseqrecord(other)))
            answer = Dseqrecord(SeqRecord.__add__(self, Dseqrecord(other)))
            answer.n = self.n
        return answer

    def __mul__(self, number):
        if not isinstance(number, int):
            raise TypeError("TypeError: can't multiply Dseqrecord by non-int of type {}".format(type(number)))
        if self.circular:
            raise TypeError("TypeError: can't multiply circular Dseqrecord")
        if number>0:
            new = copy.copy(self)
            for i in range(1, number):
                new += self
            return new
        else:
            return self.__class__("")

    def __getitem__(self, sl):
        answer = copy.copy(self)
        answer.seq = answer.seq.__getitem__(sl)
        answer.seq.alphabet = self.seq.alphabet
        if self.linear or sl.start<sl.stop:
            answer.features = SeqRecord.__getitem__(self, sl).features
        else:
            try:
                answer.features = self.shifted(sl.stop).features
            except Exception:
                answer.features = self.features
            answer.features = [f for f in answer.features if f.location.parts == sorted(f.location.parts)]
        return answer

    def __eq__( self, other ):
        try:
            if self.seq == other.seq and str(self.__dict__) == str(other.__dict__):
                return True
        except AttributeError:
            pass
        return False

    def __ne__( self, other ):
        return not self.__eq__(other)

    def linearize(self, *enzymes):
        '''This method is similar to :func:`cut` but throws an exception if there
        is not excactly on cut i.e. none or more than one digestion products.

        '''

        if self.seq._linear:
            raise Exception("Can only linearize circular molecules!")
        fragments = self.cut(*enzymes)
        if len(fragments)>1:
            raise Exception("More than one fragment is formed!")
        if not fragments:
            raise Exception("The enzyme(s) do not cut!")
        return fragments.pop()

    def cut(self, *enzymes):
        '''Digest the Dseqrecord object with one or more restriction enzymes.
        returns a list of linear Dseqrecords. If there are no cuts, an empty
        list is returned.

        See also :func:`Dseq.cut`

        Parameters
        ----------

        enzymes : enzyme object or iterable of such objects
            A Bio.Restriction.XXX restriction object or iterable of such.

        Returns
        -------
        Dseqrecord_frags : list
            list of Dseqrecord objects formed by the digestion

        Examples
        --------
        >>> import pydna
        >>> a=pydna.Dseqrecord("ggatcc")
        >>> from Bio.Restriction import BamHI
        >>> a.cut(BamHI)
        [Dseqrecord(-5), Dseqrecord(-5)]
        >>> frag1, frag2 = a.cut(BamHI)
        >>> frag1.seq
        Dseq(-5)
        g
        cctag
        >>> frag2.seq
        Dseq(-5)
        gatcc
            g


        '''

        output, stack = [], []
        stack.extend(reversed(enzymes))
        while stack:
            top = stack.pop()
            if hasattr(top, "__iter__"):
                stack.extend(reversed(top))
            else:
                output.append(top)
        enzymes = output
        if not hasattr(enzymes, '__iter__'):
            enzymes = (enzymes,)

        frags = self.seq.cut(enzymes)

        if not frags:
            return []

        if self.linear:
            last_pos=0
            #template = self.__class__(self, linear=True)
            #template = copy.copy(self)
            template = self
        else:
            last_pos = [p.pop(0)-max(0,enzyme.ovhg)-1 for (p,e) in
                         sorted([(enzyme.search(Seq(self.seq.dsdata),
                                                linear = self.linear)[::-1],
                                   enzyme) for enzyme in enzymes]) if p]
            if not last_pos:
                return [self]
            if 0 in last_pos:
                last_pos=0
            else:
                last_pos = last_pos.pop()
            template = self._multiply_circular(3)

        Dseqrecord_frags = []
        start = last_pos

        for f in frags:

            end = start + len(str(f))
            Dseqrecord_frag = Dseqrecord(f, linear=True, n=self.n)

            Dseqrecord_frag.features = template[start:end].features
            Dseqrecord_frag.annotations         = copy.copy(self[start:end].annotations)
            Dseqrecord_frag.name                = copy.copy(self.name)
            Dseqrecord_frag.dbxrefs             = copy.copy(self[start:end].dbxrefs)
            Dseqrecord_frag.id                  = copy.copy(self.id)
            Dseqrecord_frag.letter_annotations  = copy.copy(self[start:end].letter_annotations)

            Dseqrecord_frag.description = self.description+"_"+"_".join(str(e) for e in enzymes)

            Dseqrecord_frags.append(Dseqrecord_frag)
            start = end
            start-= len(f.three_prime_end()[1])

        return Dseqrecord_frags

    def number_of_cuts(self, *enzymes):
        output = []
        stack = []
        stack.extend(reversed(enzymes))
        while stack:
            top = stack.pop()
            if hasattr(top, "__iter__"):
                stack.extend(reversed(top))
            else:
                output.append(top)
        enzymes = output
        if not hasattr(enzymes, '__iter__'):
            enzymes = (enzymes,)
        return sum([len(enzyme.search(self.seq)) for enzyme in enzymes])

    def reverse_complement(self):
        '''Returns the reverse complement.

        Examples
        --------
        >>> import pydna
        >>> a=pydna.Dseqrecord("ggaatt")
        >>> a
        Dseqrecord(-6)
        >>> a.seq
        Dseq(-6)
        ggaatt
        ccttaa
        >>> a.reverse_complement().seq
        Dseq(-6)
        aattcc
        ttaagg
        >>>

        See also
        --------
        pydna.dsdna.Dseq.reverse_complement

        '''

        return self.rc()

    def rc(self):
        '''alias of the reverse_complement method'''
        answer = Dseqrecord(super(Dseqrecord, self).reverse_complement())
        assert answer.circular == self.circular
        answer.name       = "{}_rc".format(self.name[:13])
        answer.description= self.description+"_rc"
        answer.id         = self.id+"_rc"
        return answer
        #return Dseqrecord(self.seq.rc())


    def _multiply_circular(self, multiplier):
        '''returns a linearised version of a circular sequence multiplied by
       multiplier '''

        if self.linear:
            raise TypeError("sequence has to be circular!")
        if not isinstance(multiplier, int):
            raise TypeError("TypeError: can't multiply Dseq by non-int of type {}".format(type(multiplier)))
        if multiplier<=0:
            return self.__class__("")

        new_features = []

        for feature in self.features:
            new_feature = copy.deepcopy(feature)
            if len(new_feature.location.parts)>1:    # CompoundFeature
                j=0
                while (j+1)<=len(new_feature.location.parts):
                    if new_feature.location.parts[j].end == len(self) and new_feature.location.parts[j+1].start==0:
                        new_feature.location.parts[j] = FeatureLocation(new_feature.location.parts[j].start,
                                                                        new_feature.location.parts[j].end+len(new_feature.location.parts[j+1]))
                        del new_feature.location.parts[j+1]
                    j+=1
                slask = [new_feature.location.parts.pop(0)]
                for fl in new_feature.location.parts:
                    if fl.start < slask[-1].start:
                        slask.append(fl+len(self))
                    else:
                        slask.append(fl)
                if len(slask)>1:
                    new_feature.location.parts=slask
                else:
                    new_feature.location=slask[0]
            new_features.append(new_feature)

        sequence = self.tolinear()
        sequence.features = new_features
        sequence = sequence * multiplier

        sequence.features = [f for f in sequence.features if f.location.end <= len(sequence)]
        sequence.features.sort(key = operator.attrgetter('location.start'))
        return sequence

    def shifted(self, shift):
        '''Returns a circular Dseqrecord with a new origin <shift>.
         This only works on circular Dseqrecords. If we consider the following
         circular sequence:


         | ``GAAAT   <-- watson strand``
         | ``CTTTA   <-- crick strand``

         The T and the G on the watson strand are linked together as well
         as the A and the C of the of the crick strand.

         if ``shift`` is 1, this indicates a new origin at position 1:

         |    new origin at the | symbol:
         |
         | ``G|AAAT``
         | ``C|TTTA``

         new sequence:

         | ``AAATG``
         | ``TTTAC``

         Shift is always positive and 0<shift<length, so in the example
         below, permissible values of shift are 1,2 and 3

         Examples
         --------

         >>> import pydna
         >>> a=pydna.Dseqrecord("aaat",circular=True)
         >>> a
         Dseqrecord(o4)
         >>> a.seq
         Dseq(o4)
         aaat
         ttta
         >>> b=a.shifted(1)
         >>> b
         Dseqrecord(o4)
         >>> b.seq
         Dseq(o4)
         aata
         ttat

        '''
        if self.linear:
            raise Exception("Sequence is linear.\n"
                             "The origin can only be\n"
                             "shifted for a circular sequence!\n")

        length=len(self)

        if not 0<shift<length:
            raise Exception("shift is {}, has to be 0<shift<{}".format(shift, length))

        new = self._multiply_circular(3)[shift:]

        features_to_fold = [f for f in new.features if f.location.start<length<f.location.end]

        folded_features = []

        for feature in features_to_fold:

            if len(feature.location.parts)>1: # CompoundFeature
                nps=[]
                for part in feature.location.parts:

                    if part.start<part.end<=length:
                        nps.append(part)

                    elif part.start<length<part.end:
                        nps.append(FeatureLocation(part.start,length))
                        nps.append(FeatureLocation(0, part.end-length))

                    elif length<=part.start<part.end:
                        nps.append(FeatureLocation(part.start-length, part.end-length))

                folded_features.append(SeqFeature(CompoundLocation(nps),
                                                  qualifiers = feature.qualifiers,
                                                  type=feature.type))

            else:
                folded_features.append(SeqFeature(CompoundLocation([FeatureLocation(feature.location.start, length),
                                                                    FeatureLocation(0, feature.location.end-length)]),
                                                  qualifiers = feature.qualifiers,
                                                  type=feature.type))

        new = new[:length].looped()
        new.features.extend(folded_features)
        new.features.sort(key = operator.attrgetter('location.start'))
        new.description = self.description #!
        return new

    def synced(self, ref, limit = 25):
        '''This function returns a new circular sequence (Dseqrecord object), which has ben rotated
        in such a way that there is maximum overlap between the sequence and
        ref, which may be a string, Biopython Seq, SeqRecord object or
        another Dseqrecord object.

        The reason for using this could be to rotate a recombinant plasmid so
        that it starts at the same position after cloning. See the example below:


        Examples
        --------

        >>> import pydna
        >>> a=pydna.Dseqrecord("gaat",circular=True)
        >>> a.seq
        Dseq(o4)
        gaat
        ctta
        >>> d = a[2:] + a[:2]
        >>> d.seq
        Dseq(-4)
        atga
        tact
        >>> insert=pydna.Dseqrecord("CCC")
        >>> recombinant = (d+insert).looped()
        >>> recombinant.seq
        Dseq(o7)
        atgaCCC
        tactGGG
        >>> recombinant.synced(a).seq
        Dseq(o7)
        gaCCCat
        ctGGGta

        '''


        if self.linear:
            raise Exception("Only circular DNA can be synced!")
        try:
            rs = ref.seguid()
        except AttributeError:
            rs = seg(ref)

        refresh = False
        cached  = None

        csh = os.environ["pydna_cache"]

        key = str(self.seguid())+"|"+rs+"|"+str(limit)

        if csh in ("compare", "cached"):
            cache = shelve.open(os.path.join(os.environ["pydna_data_dir"],"synced"), protocol=cPickle.HIGHEST_PROTOCOL, writeback=False)
            try:
                cached = cache[str(key)]
            except:
                if os.environ["pydna_cache"] == "compare":
                    raise Exception("no result for this key!")
                else:
                    refresh = True

        if refresh or os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):

            newseq = copy.copy(self)

            s    = str(self.seq.watson).lower()
            s_rc = str(self.seq.crick).lower()

            if hasattr(ref, "seq"):
                r=ref.seq
                if hasattr(ref, "watson"):
                    r = str(r.watson).lower()
                else:
                    r = str(r).lower()
            else:
                r = str(ref.lower())

            try:
                circular_ref = ref.circular
            except AttributeError:
                circular_ref = False

            lim = min(limit, limit*(len(s)/limit)+1)

            c = common_sub_strings(s+s,       r, limit = lim)
            d = common_sub_strings(s_rc+s_rc, r, limit = lim)

            c =  [(x[0],x[2]) for x in c if x[1]==0]
            d =  [(x[0],x[2]) for x in d if x[1]==0]

            if not c and not d:
                raise Exception("There is no overlap between sequences!")

            if c:
                start, length = c.pop(0)
            else:
                start, length = 0,0

            if d:
                start_rc, length_rc = d.pop(0)
            else:
                start_rc, length_rc = 0,0

            if length_rc>length:
                start = start_rc
                newseq = newseq.rc()

            if start == 0:
                result = newseq
            else:
                result = newseq.shifted(start)

        if os.environ["pydna_cache"] == "compare":
            if result!=cached:
                module_logger.warning('dsdna error')

        if refresh or os.environ["pydna_cache"] == "refresh":
            cache[key] = result
        elif cached and os.environ["pydna_cache"] not in ("nocache", "refresh"):
            result = cached
            cache.close()

        return result




def read(data, ds = True):
    '''This function is similar the :func:`parse` function but expects one and only
    one sequence or and exception is thrown.

    Parameters
    ----------
    data : string
        see below
    ds : bool
        Double stranded or single stranded DNA, Return
        "Dseqrecord" or "SeqRecord" objects.

    Returns
    -------
    Dseqrecord
        contains the first Dseqrecord or SeqRecord object parsed.

    Notes
    -----

    The data parameter is similar to the data parameter for :func:`parse`.

    See Also
    --------
    parse

    '''

    results = parse(data, ds)
    try:
        results = results.pop()
    except IndexError:
        raise ValueError("No sequences found in data:\n({})".format(data[:79]))
    return results




def parse2(data, ds = True):
    '''experimental'''

    pattern =  r"(?:>.+\n^(?:^[^>]+?)(?=\n\n|>|LOCUS|ID))|(?:(?:LOCUS|ID)(?:(?:.|\n)+?)^//)"

    def extract_seqs(raw):
        raw = textwrap.dedent(raw).strip()
        raw = raw.replace( '\r\n', '\n')
        raw = raw.replace( '\r',   '\n')
        return re.findall(pattern, textwrap.dedent(raw+ "\n\n"),re.MULTILINE)

    files=[]
    rawseqs=[]

    if not hasattr(data, '__iter__'):
        data = (data,)
    for item in data:
        for pth in glob.glob(item):
            if os.path.isfile(pth):
                files.append(os.path.abspath(pth))
            else:
                for dirpath,_,filenames in os.walk(pth):
                    for f in filenames:
                        files.append( os.path.abspath(os.path.join(dirpath, f)))
            for file_ in files:
                with open(file_,'r') as f:
                    rawseqs.extend(extract_seqs(f.read()))
            files=[]
        else:
            rawseqs.extend(extract_seqs(item))
    sequences = []
    while rawseqs:
        circular = False
        rawseq = rawseqs.pop(0)
        handle = StringIO.StringIO(rawseq)
        try:
            parsed = SeqIO.read(handle, "embl", alphabet=IUPACAmbiguousDNA())
            #original_format = "embl"
            if "circular" in rawseq.splitlines()[0]:
                circular = True
        except ValueError:
            handle.seek(0)
            try:
                parsed = SeqIO.read(handle, "genbank", alphabet=IUPACAmbiguousDNA())
                #original_format = "genbank"
                handle.seek(0)
                parser = RecordParser()
                residue_type = parser.parse(handle).residue_type
                if "circular" in residue_type:
                    circular = True
            except ValueError:
                handle.seek(0)
                try:
                    parsed = SeqIO.read(handle, "fasta", alphabet=IUPACAmbiguousDNA())
                    if "circular" in rawseq.splitlines()[0]:
                        circular = True
                except ValueError:
                    continue

        if ds:
            sequences.append(Dseqrecord(parsed, circular = circular))
        else:
            sequences.append(parsed)
        handle.close()
    return sequences

def parse(data, ds = True):
    '''This function returns *all* DNA sequences found in data. If no
    sequences are found, an empty list is returned. This is a greedy
    function, use carefully.

    Parameters
    ----------
    data : string or iterable
        The data parameter is a string containing:

        1. an absolute path to a local file.
           The file will be read in text
           mode and parsed for EMBL, FASTA
           and Genbank sequences.

        2. an absolute path to a local directory.
           All files in the directory will be
           read and parsed as in 1.

        3. a string containing one or more
           sequences in EMBL, GENBANK,
           or FASTA format. Mixed formats
           are allowed.

        4. data can be a list or other iterable of 1 - 3

    ds : bool
        If True double stranded :class:`Dseqrecord` objects are returned.
        If False single stranded :class:`Bio.SeqRecord` [#]_ objects are returned.

    Returns
    -------
    list
        contains Dseqrecord or SeqRecord objects

    References
    ----------

    .. [#] http://biopython.org/wiki/SeqRecord

    See Also
    --------
    read

    '''
    raw= ""
    if not hasattr(data, '__iter__'):
        data = (data,)

    for item in data:
        #fn = os.path.join(dr, item )
        try:
            with open(item, 'rU') as f: raw+= f.read()
        except IOError:
            raw+=textwrap.dedent(item).strip()

    pattern =  r"(?:>.+\n^(?:^[^>]+?)(?=\n\n|>|LOCUS|ID))|(?:(?:LOCUS|ID)(?:(?:.|\n)+?)^//)"
    #raw = raw.replace( '\r\n', '\n')
    #raw = raw.replace( '\r',   '\n')
    rawseqs = re.findall(pattern, textwrap.dedent(raw + "\n\n"), flags=re.MULTILINE)
    sequences=[]

    while rawseqs:
        circular = False
        rawseq = rawseqs.pop(0)
        handle = StringIO.StringIO(rawseq)
        try:
            parsed = SeqIO.read(handle, "embl", alphabet=IUPACAmbiguousDNA())
            if "circular" in rawseq.splitlines()[0]:
                circular = True
        except ValueError:
            handle.seek(0)
            try:
                parsed = SeqIO.read(handle, "genbank", alphabet=IUPACAmbiguousDNA())
                handle.seek(0)
                parser = RecordParser()
                residue_type = parser.parse(handle).residue_type
                if "circular" in residue_type:
                    circular = True
            except ValueError:
                handle.seek(0)
                try:
                    parsed = SeqIO.read(handle, "fasta", alphabet=IUPACAmbiguousDNA())
                    if "circular" in rawseq.splitlines()[0]:
                        circular = True
                except ValueError:
                    continue

        if ds:
            sequences.append( Dseqrecord(parsed, circular = circular) )
        else:
            sequences.append( parsed )
        handle.close()
    #sequences[0].features[8].qualifiers['label'][0] = u'björn'
    return sequences

    #http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?
    #db=nuccore&
    #id=21614549&
    #strand=1&
    #seq_start=1&
    #seq_stop=100&
    #rettype=gb&
    #retmode=text

if __name__=="__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)

    #b= read(">a\naaa")

