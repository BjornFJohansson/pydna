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

import textwrap   as _textwrap
import copy       as _copy
import logging    as _logging


_module_logger = _logging.getLogger("pydna."+__name__)


from Bio.SeqRecord             import SeqRecord   as _SeqRecord
from pydna.dseqrecord          import Dseqrecord  as _Dseqrecord
from pydna._pretty             import pretty_str  as _pretty_str
from pydna.tm                  import ta_default   as _ta_default
from pydna.tm                  import ta_dbd       as _ta_dbd


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


    '''


    def __init__(    self,
                     record,
                     *args,
                     template=None,
                     forward_primer=None,
                     reverse_primer=None,
                     ta_func =_ta_default,
                     ta_func_dbd =_ta_dbd,
                     **kwargs):
        

        super().__init__(record, *args)
        self.template            = template
        self.forward_primer      = forward_primer
        self.reverse_primer      = reverse_primer
        self.ta_func             = ta_func 
        self.ta_func_dbd         = ta_func_dbd  
        self.__dict__.update(kwargs)
        
        # https://medium.com/@chipiga86/circular-references-without-memory-leaks-and-destruction-of-objects-in-python-43da57915b8d


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
        return answer

    
    rc = reverse_complement

    
    def ta(self):
        return self.ta_func(self.forward_primer.footprint, self.reverse_primer.footprint, str(self.seq))
    
    
    def ta_dbd(self):
        return self.ta_func_dbd(self.forward_primer.footprint, self.reverse_primer.footprint, str(self.seq))


    def figure(self):
        '''
        This method returns a simple figure of the two primers binding to a part
        of the template.

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


        f =   '''
            {sp1}5{faz}...{raz}3
             {sp3}{rap}
            {sp3}3{rp}5
            5{fp}3
             {fap:>{fplength}}
            {sp2}3{fzc}...{rzc}5
            '''.format( fp       = self.forward_primer.seq,
                        fap      = "|"*len(self.forward_primer.footprint),
                        fplength = len(self.forward_primer.seq),
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

        r'''Returns a string containing a text representation of a suggested
       PCR program using Taq or similar polymerase.

       ::

        |95°C|95°C               |    |tmf:59.5
        |____|_____          72°C|72°C|tmr:59.7
        |5min|30s  \ 59.1°C _____|____|30s/kb
        |    |      \______/ 0:32|5min|GC 51%
        |    |       30s         |    |1051bp

       '''


        # Taq polymerase extension rate is set to 30 nt/s
        # see https://www.thermofisher.com/pt/en/home/life-science/pcr/pcr-enzymes-master-mixes/taq-dna-polymerase-enzymes/taq-dna-polymerase.html
        taq_extension_rate = 30  # seconds/kB PCR product length
        extension_time_taq = int(round(taq_extension_rate * len(self) / 1000)) # seconds

        f=_textwrap.dedent(    r'''
                                |95°C|95°C               |    |tmf:{tmf:.1f}
                                |____|_____          72°C|72°C|tmr:{tmr:.1f}
                                |5min|30s  \ {ta:.1f}°C _____|____|{rate}s/kb
                                |    |      \______/{0:2}:{1:2}|5min|GC {GC_prod}%
                                |    |       30s         |    |{size}bp
                                '''[1:-1].format(rate = taq_extension_rate,
                                        size= len(self.seq),
                                        ta=round(self.ta(),1),
                                        tmf=self.forward_primer.tm(),
                                        tmr=self.reverse_primer.tm(),
                                        GC_prod= int(self.gc()),
                                        *map(int,divmod(extension_time_taq,60)))) 
                                                                                                
                                                                                          
        return _pretty_str(f)


    taq_program = program


    def dbd_program(self):
        r'''Returns a string containing a text representation of a suggested
       PCR program using a polymerase with a DNA binding domain such as Pfu-Sso7d.

       ::

        |98°C|98°C               |    |tmf:53.8
        |____|_____          72°C|72°C|tmr:54.8
        |30s |10s  \ 57.0°C _____|____|15s/kb
        |    |      \______/ 0:15|5min|GC 51%
        |    |       10s         |    |1051bp
        
        
        |98°C|98°C      |    |tmf:82.5
        |____|____      |    |tmr:84.4
        |30s |10s \ 72°C|72°C|15s/kb
        |    |     \____|____|GC 52%
        |    |      3:45|5min|15058bp   

       '''
        PfuSso7d_extension_rate = 15                #seconds/kB PCR product
        extension_time_PfuSso7d = PfuSso7d_extension_rate * len(self) / 1000  # seconds


        # The program returned is eaither a two step or three step progrem
        # This depends on the tm and length of the primers in the
        # original instructions from finnzyme. These do not seem to be

        # Ta calculation for enzymes with dsDNA binding domains like phusion or Pfu-Sso7d
        # https://www.finnzymes.fi/tm_determination.html
        
        
        tmf = self.forward_primer.tm_dbd()
        tmr = self.reverse_primer.tm_dbd()


        if (tmf>=69.0 and tmr>=69.0):
            f=_textwrap.dedent(    r'''
                                    |98°C|98°C      |    |tmf:{tmf:.1f}
                                    |____|____      |    |tmr:{tmr:.1f}
                                    |30s |10s \ 72°C|72°C|{rate}s/kb
                                    |    |     \____|____|GC {GC_prod}%
                                    |    |     {0:2}:{1:2}|5min|{size}bp
                                    '''[1:-1].format(rate=PfuSso7d_extension_rate,
                                                     tmf=tmf,
                                                     tmr=tmr,
                                                     GC_prod=int(self.gc()),
                                                     size=len(self.seq),
                                                     *map(int,divmod(extension_time_PfuSso7d,60)),))   
        else:
            f=_textwrap.dedent(    r'''
                                    |98°C|98°C               |    |tmf:{tmf:.1f}
                                    |____|_____          72°C|72°C|tmr:{tmr:.1f}
                                    |30s |10s  \ {ta:.1f}°C _____|____|{rate}s/kb
                                    |    |      \______/{0:2}:{1:2}|5min|GC {GC_prod}%
                                    |    |       10s         |    |{size}bp
                                    '''[1:-1].format(rate = PfuSso7d_extension_rate,
                                            size= len(self.seq),
                                            ta   = round(self.ta_dbd()),
                                            tmf=tmf,
                                            tmr=tmr,
                                            GC_prod= int(self.gc()),
                                            *map(int,divmod(extension_time_PfuSso7d,60))) )      
        return _pretty_str(f)


    pfu_sso7d_program = dbd_program


if __name__=="__main__":
    import os as _os
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"]=""
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"]=cached

