#!/usr/bin/env python3
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

import copy              as _copy
import itertools         as _itertools

import sys               as _sys
import math              as _math

from Bio.Alphabet.IUPAC     import IUPACAmbiguousDNA as _IUPACAmbiguousDNA
from Bio.Seq                import Seq               as _Seq
from pydna.utils  import seguid                           as _seg
from pydna.utils  import rc                               as _rc
from pydna.common_sub_strings import common_sub_strings as _common_sub_strings

class Dseq(_Seq):
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

    If the ovhg parameter is psecified a crick strand also needs to be supplied,
    otherwise an exception is raised.

    >>> Dseq(watson="agt",ovhg=2)
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


    >>> Dseq("aaa","ttt")
    Dseq(-3)
    aaa
    ttt
    >>> Dseq("aaa","ttt",ovhg=0)
    Dseq(-3)
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

    Coercing to string

    >>> a=Dseq("tttcccc","aaacccc")
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

    '''

    def __init__(self,
                  watson,
                  crick         = None,
                  ovhg          = None,
                  linear        = None,
                  circular      = None,
                  alphabet      = _IUPACAmbiguousDNA() ):

        watson = "".join(watson.split())

        if ovhg is None:
            if crick is None:
                self.crick = _rc(watson)
                self._ovhg = 0
            else:
                crick = "".join(crick.split())

                olaps = _common_sub_strings(str(watson).lower(),
                                            str(_rc(crick).lower()),
                                            int(_math.log(len(watson))/_math.log(4)))
                try:
                    F,T,L = olaps[0]
                except IndexError:
                    raise Exception("Could not anneal the two strands! "
                                    "ovhg should be provided")
                ovhgs = [ol[1]-ol[0] for ol in olaps if ol[2]==L]
                if len(ovhgs)>1:
                    for o in ovhgs:
                        print(o)
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
        asn = ((-self._ovhg*" ") + str(_rc(self.crick)))

        self.todata = "".join([a.strip() or b.strip() for a,b in _itertools.zip_longest(sns,asn, fillvalue=" ")])
        self.dsdata = "".join([a for a, b in _itertools.zip_longest(sns,asn, fillvalue=" ") if a.lower()==b.lower()])

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

        _Seq.__init__(self, self.todata, alphabet)

    def mw(self):
        nts = ( self.watson + self.crick ).lower()
        # (A x 313.2) + (T x 304.2) + (C x 289.2) + (G x 329.2) + 79.0
        return sum( [313.2 for n in nts if n=="a"] +
                    [304.2 for n in nts if n=="t"] +
                    [289.2 for n in nts if n=="c"] +
                    [329.2 for n in nts if n=="g"] +
                    [308.9 for n in nts if n=="n"]) + 79



    def find(self, sub, start=0, end=_sys.maxsize):
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
                    other.ovhg == self._ovhg and
                    self.circular == other.circular)
                    # Also test for alphabet ?
        except AttributeError:
            same = False
        return same


    def fig(self):
        '''Returns a representation of the sequence, truncated if
       longer than 30 bp:

       Examples
       --------

       >>> from pydna.dseq import Dseq
       >>> a=Dseq("atcgcttactagcgtactgatcatctgact")
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
        '''Returns a Dseq object where watson and crick have switched
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

       '''
        ovhg = len(self.watson) - len(self.crick) + self._ovhg
        return Dseq(self.crick, self.watson, ovhg=ovhg, circular = self.circular)

    def reverse_complement(self):
        '''Alias of the rc method'''
        return self.rc()

    def looped(self):
        '''Returns a circularized Dseq object. This can only be done if the
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

       '''
        if self.circular:
            return self
        type5, sticky5 = self.five_prime_end()
        type3, sticky3 = self.three_prime_end()
        if type5 == type3 and str(sticky5) == str(_rc(sticky3)):
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

       >>> from pydna.dseq import Dseq
       >>> a=Dseq("catcgatc", circular=True)
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
            str(self_tail) == str(_rc(other_tail))):
            answer = Dseq(self.watson + other.watson,
                          other.crick + self.crick,
                          self._ovhg,)
        elif not self:
            answer = _copy.copy(other)
        elif not other:
            answer = _copy.copy(self)
        else:
            raise TypeError("sticky ends not compatible!")
        return answer

    def __mul__(self, number):
        if not isinstance(number, int):
            raise TypeError("TypeError: can't multiply Dseq by non-int of type {}".format(type(number)))
        if number<=0:
            return self.__class__("")
        new = _copy.copy(self)
        for i in range(number-1):
            new += self
        return new

    def _fill_in_five_prime(self, nucleotides):
        stuffer = ''
        type, se = self.five_prime_end()
        if type == "5'":
            for n in _rc(se):
                if n in nucleotides:
                    stuffer+=n
                else:
                    break
        return self.crick+stuffer, self._ovhg+len(stuffer)

    def _fill_in_three_prime(self, nucleotides):
        stuffer = ''
        type, se = self.three_prime_end()
        if type == "5'":
            for n in _rc(se):
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

        enzymes = [e for (p,e) in sorted([(enzyme.search(_Seq(frags[0].dsdata))[::-1], enzyme) for enzyme in enzymes], reverse=True) if p]

        if not enzymes:
            return [self,]

        for enzyme in enzymes:
            for frag in frags:

                if enzyme.search(_Seq(frag.dsdata)):

                    watson_fragments = [str(s) for s in enzyme.catalyze(_Seq(frag.watson+"N"))]
                    crick_fragments  = [str(s) for s in enzyme.catalyze(_Seq(frag.crick+"N" ))[::-1]]

                    watson_fragments[-1] = watson_fragments[-1][:-1]
                    crick_fragments[0]   = crick_fragments[0][:-1]

                    s = list(zip(watson_fragments, crick_fragments))

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

        >>> from pydna.dseq import Dseq
        >>> seq=Dseq("ggatccnnngaattc")
        >>> seq
        Dseq(-15)
        ggatccnnngaattc
        cctaggnnncttaag
        >>> from Bio.Restriction import BamHI,EcoRI
        >>> type(seq.cut(BamHI))
        <class 'list'>
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

        enzymes = [e for (p,e) in sorted([(enzyme.search(_Seq(frags[0].dsdata))[::-1], enzyme) for enzyme in enzymes], reverse=True) if p]


        if not enzymes:
            return []

        for enz in enzymes:
            for frag in frags:

                ws = [x-1 for x in enz.search(_Seq(frag.watson)+"N")] #, linear = frag.linear
                cs = [x-1 for x in enz.search(_Seq(frag.crick) +"N")] #, linear = frag.linear

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
        if self.ovhg==rc_ovhg==0:
            return _seg(min(self.watson, self.crick))
        if self.ovhg<rc_ovhg:
            w = self.watson
            c = self.crick
            o = self.ovhg
        elif self.ovhg>rc_ovhg:
            w = self.crick
            c = self.watson
            o = rc_ovhg
        elif self.ovhg==rc_ovhg:
            w, c = sorted((self.watson, self.crick))
            o = self.ovhg
        return _seg( str(o) + w + "|" + c)

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

if __name__=="__main__":
    import os        as _os
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
