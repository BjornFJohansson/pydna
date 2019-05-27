#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
# doctest: +NORMALIZE_WHITESPACE 
# doctest: +SKIP
"""This module provides the :class:`Amplicon` class for PCR simulation. This class is not meant to be use directly
but is used by the :mod:`amplify` module"""

import math       as _math
import itertools  as _itertools
import re         as _re
import textwrap   as _textwrap
import copy       as _copy
import logging    as _logging
_module_logger = _logging.getLogger("pydna."+__name__)

from Bio.SeqRecord                  import SeqRecord      as _SeqRecord
from Bio.SeqUtils                   import GC             as _GC
from Bio.SeqUtils.MeltingTemp       import Tm_NN          as _Tm_NN
from pydna.utils                         import rc             as _rc
from pydna.dseqrecord                    import Dseqrecord     as _Dseqrecord
from pydna._pretty                       import pretty_str     as _pretty_str
from pydna.tm                            import tmbresluc      as _tmbresluc


class Amplicon(_Dseqrecord):
    '''The Amplicon class holds information about a PCR reaction involving two
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

    saltc : float, optional
        saltc = monovalent cations (mM) (Na,K..)
        default value is 50mM
        This is used for Tm calculations.

    forward_primer_concentration : float, optional
        primer concentration (nM)
        default set to 1000nM = 1µM
        This is used for Tm calculations.

    rc : float, optional
        primer concentration (nM)
        default set to 1000nM = 1µM
        This is used for Tm calculations.
    '''

    def __init__(    self,
                     record,
                     *args,
                     template=None,
                     forward_primer=None,
                     reverse_primer=None,
                     saltc=50.0,
                     fprimerc=1000.0,
                     rprimerc=1000.0,
                     **kwargs):

        super().__init__(record, *args, **kwargs)
        self.template       = template
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.fprimerc       = fprimerc
        self.rprimerc       = rprimerc
        self.saltc          = saltc

    @classmethod
    def from_SeqRecord(cls, record, *args, path=None, **kwargs):
        obj = super().from_SeqRecord(record, *args, **kwargs)
        obj.path = path
        return obj

    def __getitem__(self, sl):
        answer = _copy.copy(self)
        answer.seq = answer.seq.__getitem__(sl)
        answer.seq.alphabet = self.seq.alphabet
        sr = _SeqRecord("n"*len(self))
        sr.features = self.features
        answer.features = _SeqRecord.__getitem__(sr, sl).features
        return answer

    def __repr__(self):
        '''returns a short string representation of the object'''
        return "Amplicon({})".format(self.__len__())

    def _repr_pretty_(self, p, cycle):
        p.text("Amplicon({})".format(self.__len__()))
            
    def _repr_html_(self):
        return "Amplicon({})".format(self.__len__())

    def reverse_complement(self):
        answer = type(self)(super().reverse_complement())
        answer.template       = self.template.rc()
        answer.forward_primer = self.reverse_primer
        answer.reverse_primer = self.forward_primer
        answer.fprimerc       = self.rprimerc
        answer.rprimerc       = self.fprimerc
        answer.saltc          = self.saltc
        return answer
    
    rc = reverse_complement

    def figure(self):
        '''
        This method returns a simple figure of the two primers binding to a part
        of the template.

        ::

         5gctactacacacgtactgactg3
          |||||||||||||||||||||| tm 52.6 (dbd) 58.3
         5gctactacacacgtactgactg...caagatagagtcagtaaccaca3
         3cgatgatgtgtgcatgactgac...gttctatctcagtcattggtgt5
                                   |||||||||||||||||||||| tm 49.1 (dbd) 57.7
                                  3gttctatctcagtcattggtgt5



        Returns
        -------
        figure:string
             A string containing a text representation of the primers
             annealing on the template (see example above).


        Notes
        -----
        tm in the figure above is the melting temperature (tm) for each primer calculated according to
        SantaLucia 1998 [#]_.

        dbd is the tm calculation for enzymes with dsDNA binding domains like Pfu-Sso7d [#]_. See [#]_
        for more information.

        References
        ----------

        .. [#] J. SantaLucia, “A Unified View of Polymer, Dumbbell, and Oligonucleotide DNA Nearest-neighbor Thermodynamics,” Proceedings of the National Academy of Sciences 95, no. 4 (1998): 1460.
        .. [#] M. Nørholm, “A Mutant Pfu DNA Polymerprimerase Designed for Advanced Uracil-excision DNA Engineering,” BMC Biotechnology 10, no. 1 (2010): 21, doi:10.1186/1472-6750-10-21.
        .. [#] http://www.thermoscientificbio.com/webtools/tmc/

        '''

        tmf = _Tm_NN(str(self.forward_primer.footprint),
                    dnac1=self.fprimerc,
                    Na=self.saltc)
        tmr = _Tm_NN(str(self.reverse_primer.footprint),
                    dnac1=self.fprimerc,
                    Na=self.saltc)

        tmf_dbd = _tmbresluc(str(self.forward_primer.footprint), primerc=self.fprimerc)
        tmr_dbd = _tmbresluc(str(self.reverse_primer.footprint), primerc=self.rprimerc)

        f =   '''
            {sp1}5{faz}...{raz}3
             {sp3}{rap} tm {tmr} (dbd) {tmr_dbd}
            {sp3}3{rp}5
            5{fp}3
             {fap:>{fplength}} tm {tmf} (dbd) {tmf_dbd}
            {sp2}3{fzc}...{rzc}5
            '''.format( fp       = self.forward_primer.seq,
                        fap      = "|"*len(self.forward_primer.footprint),
                        fplength = len(self.forward_primer.seq),
                        tmf      = round(tmf,1),
                        tmr      = round(tmr,1),
                        tmf_dbd  = round(tmf_dbd,1),
                        tmr_dbd  = round(tmr_dbd,1),
                        rp       = self.reverse_primer.seq[::-1],
                        rap      = "|"*len(self.reverse_primer.footprint),
                        rplength = len(self.reverse_primer.seq),
                        faz      = self.forward_primer.footprint,
                        raz      = self.reverse_primer.footprint.reverse_complement(),
                        fzc      = self.forward_primer.footprint.complement(),
                        rzc      = self.reverse_primer.footprint[::-1],
                        sp1      = " "*(len(self.forward_primer.seq)-len(self.forward_primer.footprint)),
                        sp2      = " "*(len(self.forward_primer.seq)-len(self.forward_primer.footprint)),
                        sp3      = " "*(3+len(self.forward_primer.seq))
                       )
        return _pretty_str(_textwrap.dedent(f).strip("\n"))


    def program(self):

        r'''Returns a string containing a text representation of a proposed
       PCR program using Taq or similar polymerase.

       ::

        Taq (rate 30 nt/s)
        Three-step|         30 cycles     |      |Tm formula: Biopython Tm_NN
        94.0°C    |94.0°C                 |      |SaltC 50mM
        __________|_____          72.0°C  |72.0°C|
        04min00s  |30s  \         ________|______|
                  |      \ 46.0°C/ 0min 1s|10min |
                  |       \_____/         |      |
                  |         30s           |      |4-8°C

       '''


        # Primer melting temperatures are calculated with the Tm_NN formula from
        # biopython
        # simple salt concentration correction is used and the template concentration
        # is ignored. dnac1 = primer concentration
        
        tmf = _Tm_NN(str(self.forward_primer.footprint),
                    dnac1=self.fprimerc,
                    Na=self.saltc)
        tmr = _Tm_NN(str(self.reverse_primer.footprint),
                    dnac1=self.fprimerc,
                    Na=self.saltc)

        # Ta calculation according to
        # Rychlik, Spencer, and Rhoads, 1990, Optimization of the anneal
        # ing temperature for DNA amplification in vitro
        # http://www.ncbi.nlm.nih.gov/pubmed/2243783
        # The formula described uses the length and GC content of the product and
        # salt concentration (monovalent cations).

        #GC_prod=GC(str(self.seq))

        tmp = 81.5 + 0.41*_GC(str(self.seq)) + 16.6*_math.log10(self.saltc/1000.0) - 675/len(self)
        tml = min(tmf,tmr)
        ta = 0.3*tml+0.7*tmp-14.9

        # Taq polymerase extension rate is set to 30 nt/s
        # see https://www.thermofisher.com/pt/en/home/life-science/pcr/pcr-enzymes-master-mixes/taq-dna-polymerase-enzymes/taq-dna-polymerase.html
        taq_extension_rate = 30  # seconds/kB PCR product length
        extension_time_taq = int(round(taq_extension_rate * len(self) / 1000)) # seconds
        f  = _textwrap.dedent(r'''
                                 Taq (rate {rate} nt/s) 35 cycles             |{size}bp
                                 95.0°C    |95.0°C                 |      |Tm formula: Biopython Tm_NN
                                 |_________|_____          72.0°C  |72.0°C|SaltC {saltc:2}mM
                                 | 03min00s|30s  \         ________|______|Primer1C {forward_primer_concentration:3}µM
                                 |         |      \ {ta}°C/{0:2}min{1:2}s| 5min |Primer2C {reverse_primer_concentration:3}µM
                                 |         |       \_____/         |      |GC {GC_prod}%
                                 |         |         30s           |      |4-12°C'''.format(rate=taq_extension_rate,
                                                                                            forward_primer_concentration=self.fprimerc/1000,
                                                                                            reverse_primer_concentration=self.rprimerc/1000,
                                                                                            ta=round(ta,1),
                                                                                            saltc=self.saltc,
                                                                                            *divmod(extension_time_taq,60),
                                                                                            size= len(self.seq),
                                                                                            GC_prod= int(self.gc()) ))
        return _pretty_str(f)


    def taq_program(self):
        return self.program()


    def dbd_program(self):
        r'''Returns a string containing a text representation of a proposed
       PCR program using a polymerase with a DNA binding domain such as Pfu-Sso7d.

       ::

        Pfu-Sso7d (rate 15s/kb)             |{size}bp
        Three-step|          30 cycles   |      |Tm formula: Pydna tmbresluc
        98.0°C    |98.0°C                |      |SaltC 50mM
        __________|_____          72.0°C |72.0°C|Primer1C   1µM
        00min30s  |10s  \ 61.0°C ________|______|Primer2C   1µM
                  |      \______/ 0min 0s|10min |
                  |        10s           |      |4-12°C

       '''
        PfuSso7d_extension_rate = 15                #seconds/kB PCR product
        extension_time_PfuSso7d = PfuSso7d_extension_rate * len(self) / 1000  # seconds

        # The program returned is eaither a two step or three step progrem
        # This depends on the tm and length of the primers in the
        # original instructions from finnzyme. These do not seem to be

        tmf_dbd = _tmbresluc(str(self.forward_primer.footprint), primerc=self.fprimerc)
        tmr_dbd = _tmbresluc(str(self.reverse_primer.footprint), primerc=self.rprimerc)

        # Ta calculation for enzymes with dsDNA binding domains like Pfu-Sso7d
        # https://www.finnzymes.fi/tm_determination.html

        length_of_f = len(self.forward_primer.footprint)
        length_of_r = len(self.reverse_primer.footprint)

        if (length_of_f>20 and length_of_r>20 and tmf_dbd>=69.0 and tmr_dbd>=69.0) or (tmf_dbd>=72.0 and tmr_dbd>=72.0):
            f=_textwrap.dedent( r'''
                                    Pfu-Sso7d (rate {rate}s/kb)
                                    Two-step|    30 cycles |      |{size}bp
                                    98.0°C  |98.0C         |      |Tm formula: Pydna tmbresluc
                                    _____ __|_____         |      |SaltC {saltc:2}mM
                                    00min30s|10s  \        |      |Primer1C {forward_primer_concentration:3}µM
                                            |      \ 72.0°C|72.0°C|Primer2C {reverse_primer_concentration:3}µM
                                            |       \______|______|GC {GC_prod}%
                                            |      {0:2}min{1:2}s|10min |4-12°C
                                 '''.format(rate = PfuSso7d_extension_rate,
                                            forward_primer_concentration = self.fprimerc/1000,
                                            reverse_primer_concentration = self.rprimerc/1000,
                                            saltc = self.saltc,
                                            *map(int,divmod(extension_time_PfuSso7d,60)),
                                            GC_prod= int(self.gc()),
                                            size = len(self.seq) ))
        else:

            if (length_of_f>20 and length_of_r>20):
                ta = min(tmf_dbd,tmr_dbd)+3
            else:
                ta = min(tmf_dbd,tmr_dbd)

            f=_textwrap.dedent( r'''
                                    Pfu-Sso7d (rate {rate}s/kb)                 |{size}bp
                                    Three-step|          30 cycles   |      |Tm formula: Pydna tmbresluc
                                    98.0°C    |98.0°C                |      |SaltC {saltc:2}mM
                                    __________|_____          72.0°C |72.0°C|Primer1C {forward_primer_concentration:3}µM
                                    00min30s  |10s  \ {ta:.1f}°C ________|______|Primer2C {reverse_primer_concentration:3}µM
                                              |      \______/{0:2}min{1:2}s|10min |GC {GC_prod}%
                                              |        10s           |      |4-12°C
                                 '''.format(rate = PfuSso7d_extension_rate,
                                            size= len(self.seq),
                                            ta   = round(ta),
                                            forward_primer_concentration   = self.fprimerc/1000,
                                            reverse_primer_concentration   = self.rprimerc/1000,
                                            saltc= self.saltc,
                                            GC_prod= int(self.gc()),
                                            *map(int, divmod(extension_time_PfuSso7d,60)) ))
        return _pretty_str(f)


    def pfu_sso7d_program(self):
        return self.dbd_program()




if __name__=="__main__":
    import os as _os
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"]=""
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"]=cached

