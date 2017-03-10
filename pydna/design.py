#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2013, 2014 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''
This module contain functions for primer design.

'''
import warnings

import math                                       as _math
from operator import itemgetter                   as _itemgetter
import os                                         as _os
import copy                                       as _copy
from Bio.Alphabet import Alphabet                 as _Alphabet
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA  as _IUPACAmbiguousDNA
from Bio.Seq import Seq                           as _Seq
from pydna.amplify import Anneal                  as _Anneal
from pydna.amplify import pcr                     as _pcr
from pydna.tm import tmbresluc                    as _tmbresluc
from pydna.dseqrecord import Dseqrecord           as _Dseqrecord
from pydna._pretty import pretty_str              as _pretty_str
from pydna.primer    import Primer                as _Primer

import logging    as _logging
_module_logger = _logging.getLogger("pydna."+__name__)


def print_primer_pair(*args,**kwargs):
    from pydna import _PydnaDeprecationWarning
    warnings.warn("This function has been deprecated."
                  "consider using the primer_design"
                  "function instead", _PydnaDeprecationWarning)
    f,r = cloning_primers(*args,**kwargs)
    return _pretty_str("\n"+f.format("fasta")+"\n"+r.format("fasta") + "\n")

def cloning_primers( template,
                     minlength=16,
                     maxlength=29,
                     fp=None,
                     rp=None,
                     fp_tail='',
                     rp_tail='',
                     target_tm=55.0,
                     fprimerc=1000.0,
                     rprimerc=1000.0,
                     saltc=50.0,
                     formula = _tmbresluc):
    from pydna import _PydnaDeprecationWarning
    warnings.warn("This function has been deprecated."
                  "consider using the primer_design"
                  "function instead", _PydnaDeprecationWarning)

    '''
    **Do not use this function, use pydna.design.primer_design instead**
    **This function will be deprecated and removed in a future version of pydna**
    **This can be discussed in the google group https://groups.google.com/forum/#!forum/pydna ***
    
    This function can design primers for PCR amplification of a given sequence.
    This function accepts a Dseqrecord object containing the template sequence and
    returns a tuple cntaining two ::mod`Bio.SeqRecord.SeqRecord` objects describing
    the primers.

    Primer tails can optionally be given in the form of strings.

    An predesigned primer can be given, either the forward or reverse primers. In this
    case this function tries to design a primer with a Tm to match the given primer.


    Parameters
    ----------

    template : Dseqrecord
        a Dseqrecord object. The only required argument.

    minlength : int, optional
        Minimum length of the annealing part of the primer.

    maxlength : int, optional
        Maximum length (including tail) for designed primers.

    fp, rp : SeqRecord, optional
        optional Biopython SeqRecord objects containing one primer each.

    fp_tail, rp_tail : string, optional
        optional tails to be added to the forwars or reverse primers.

    target_tm : float, optional
        target tm for the primers, set to 55°C by default.

    fprimerc : float, optional
        Concentration of forward primer in nM, set to 1000.0 nM by default.

    rprimerc : float, optional
        Concentration of reverse primer in nM, set to 1000.0 nM by default.

    saltc  : float, optional
        Salt concentration (monovalet cations) :mod:`tmbresluc` set to 50.0 mM by default

    formula : function
        formula used for tm calculation
        this is the name of a function.
        built in options are:

        1. :func:`pydna.amplify.tmbresluc` (default)
        2. :func:`pydna.amplify.basictm`
        3. :func:`pydna.amplify.tmstaluc98`
        4. :func:`pydna.amplify.tmbreslauer86`

        These functions are imported from the :mod:`pydna.amplify` module, but can be
        substituted for some other custom made function.

    Returns
    -------
    fp, rp : tuple
        fp is a :mod:Bio.SeqRecord object describing the forward primer
        rp is a :mod:Bio.SeqRecord object describing the reverse primer



    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.design import cloning_primers
    >>> from pydna.amplify import pcr
    >>> t=Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    >>> t
    Dseqrecord(-64)
    >>> pf,pr = cloning_primers(t)
    >>> pf
    fw64 17-mer:5'-atgactgctaacccttc-3'
    >>> pr
    rv64 17-mer:5'-catcgtaagtttcgaac-3'
    >>> pcr_prod = pcr(pf, pr, t)
    >>> pcr_prod
    Amplicon(64)
    >>>
    >>> print(pcr_prod.figure())
    5atgactgctaacccttc...gttcgaaacttacgatg3
                         ||||||||||||||||| tm 49.0 (dbd) 52.9
                        3caagctttgaatgctac5
    5atgactgctaacccttc3
     ||||||||||||||||| tm 51.6 (dbd) 54.0
    3tactgacgattgggaag...caagctttgaatgctac5
    >>> pf,pr = cloning_primers(t, fp_tail="GGATCC", rp_tail="GAATTC")
    >>> pf
    fw64 23-mer:5'-GGATCCatgactgct..ttc-3'
    >>> pr
    rv64 23-mer:5'-GAATTCcatcgtaag..aac-3'
    >>> pcr_prod = pcr(pf, pr, t)
    >>> print(pcr_prod.figure())
          5atgactgctaacccttc...gttcgaaacttacgatg3
                               ||||||||||||||||| tm 49.0 (dbd) 52.9
                              3caagctttgaatgctacCTTAAG5
    5GGATCCatgactgctaacccttc3
           ||||||||||||||||| tm 51.6 (dbd) 54.0
          3tactgacgattgggaag...caagctttgaatgctac5
    >>> print(pcr_prod.seq)
    GGATCCatgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatgGAATTC
    >>>
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> pf = SeqRecord(Seq("atgactgctaacccttccttggtgttg"))
    >>> pf,pr = cloning_primers(t, fp = pf, fp_tail="GGATCC", rp_tail="GAATTC")
    >>> pf
    SeqRecord(seq=Seq('GGATCCatgactgctaacccttccttggtgttg', Alphabet()), id='fw64', name='fw64', description='fw64 id?', dbxrefs=[])
    >>> pr
    rv64 34-mer:5'-GAATTCcatcgtaag..gtc-3'
    >>> ampl = pcr(pf,pr,t)
    >>> print(ampl.figure())
          5atgactgctaacccttccttggtgttg...gacgacatttcgttcgaaacttacgatg3
                                         |||||||||||||||||||||||||||| tm 61.7 (dbd) 72.2
                                        3ctgctgtaaagcaagctttgaatgctacCTTAAG5
    5GGATCCatgactgctaacccttccttggtgttg3
           ||||||||||||||||||||||||||| tm 63.7 (dbd) 72.3
          3tactgacgattgggaaggaaccacaac...ctgctgtaaagcaagctttgaatgctac5
    >>>


    '''

    if fp and not rp:
        fp = _Primer(_Seq(fp_tail, _IUPACAmbiguousDNA())) + fp
        p  = _Anneal([fp], template).forward_primers.pop()
        fp = _Primer(p.footprint)
        fp_tail = _Primer(p.tail)
        rp = _Primer(_Seq(str(template[-(maxlength*3-len(rp_tail)):].reverse_complement().seq), _IUPACAmbiguousDNA()))
        target_tm = formula(str(fp.seq).upper(), primerc=fprimerc, saltc=saltc)
    elif not fp and rp:
        rp = _Primer(_Seq(rp_tail, _IUPACAmbiguousDNA())) + rp
        p =  _Anneal([rp], template).reverse_primers.pop()
        rp = _Primer(p.footprint)
        rp_tail = _Primer(p.tail)
        fp = _Primer(_Seq(str(template[:maxlength*3-len(fp_tail)].seq), _IUPACAmbiguousDNA()))
        target_tm = formula(str(rp.seq).upper(), primerc=rprimerc, saltc=saltc)
    elif not fp and not rp:
        fp = _Primer(_Seq(str(template[:maxlength*3-len(fp_tail)].seq), _IUPACAmbiguousDNA()))
        rp = _Primer(_Seq(str(template[-maxlength*3+len(rp_tail):].reverse_complement().seq), _IUPACAmbiguousDNA()))
    else:
        raise Exception("Specify maximum one of the two primers, not both.")

    fp.concentration = fprimerc
    rp.concentration = rprimerc
    
    lowtm, hightm = sorted( [( formula(str(fp.seq), fprimerc, saltc), fp, "f" ),
                             ( formula(str(rp.seq), rprimerc, saltc), rp, "r" ) ] , key=_itemgetter(0))
                          
    while lowtm[0] > target_tm and len(lowtm[1])>minlength:
        shorter = lowtm[1][:-1]
        tm      = formula(str(shorter.seq).upper(), primerc=fprimerc, saltc=saltc)
        lowtm   = (tm, shorter, lowtm[2])

    while hightm[0] > lowtm[0] + 2.0 and len(hightm[1])>minlength:
        shorter = hightm[1][:-1]
        tm = formula(str(shorter.seq).upper(), primerc = rprimerc, saltc = saltc)
        hightm = (tm, shorter, hightm[2])

    fp, rp = sorted((lowtm, hightm), key=_itemgetter(2))

    fp = fp_tail + fp[1]
    rp = rp_tail + rp[1]

    #fp.description = "fw{}".format(len(template))
    #rp.description = "rv{}".format(len(template))

    #fp.name = "fw{}".format(len(template))[:15]
    #rp.name = "rv{}".format(len(template))[:15]

    fp.description = "fw{}".format(len(template))+' '+template.accession
    rp.description = "rv{}".format(len(template))+' '+template.accession

    fp.id = "fw{}".format(len(template))
    rp.id = "rv{}".format(len(template))

    fp.name = fp.id
    rp.name = rp.id

    if fp.seq.alphabet == _Alphabet():
        fp.seq.alphabet = _IUPACAmbiguousDNA()
    if rp.seq.alphabet == _Alphabet():
        rp.seq.alphabet = _IUPACAmbiguousDNA()

    return fp, rp

def primer_design(    template,
                      fp=None,
                      rp=None,
                      target_tm=55.0,
                      fprimerc=1000.0,  # nM
                      rprimerc=1000.0,  # nM
                      saltc=50.0,
                      formula = _tmbresluc):

    '''This function designs a forward primer and a reverse primer for PCR amplification 
    of a given template sequence.
    
    The template argument is a Dseqrecord object or equivalent containing the template sequence.
    
    The optional fp and rp arguments can contain an existing primer for the sequence (either the forward or reverse primer).
    One or the other primers can be specified, not both (since then there is nothing to design!, use the pydna.amplify.pcr function instead).
    
    If one ofthe primers is given, the other primer is designed to match in terms of Tm.
    If both primers are designed, they will be designed to target_tm
    
    fprimerc, rprimerc and saltc are formward and reverse primer concentration (nM). Saltc is the salt concentration. 
    These arguments might affect how Tm is calculated.
    
    formula is a function that can take at least three arguments f( str, primerc=float, saltc=float).
    There are several of these in the pydna.tm module.
    
    The function returns a pydna.amplicon.Amplicon class instance. This object has 
    the object.forward_primer and object.reverse_primer properties which contain the designed primers.
    
    
    Parameters
    ----------

    template : pydna.dseqrecord.Dseqrecord
        a Dseqrecord object. The only required argument.

    fp, rp : pydna.primer.Primer, optional
        optional pydna.primer.Primer objects containing one primer each.

    target_tm : float, optional
        target tm for the primers, set to 55°C by default.

    fprimerc : float, optional
        Concentration of forward primer in nM, set to 1000.0 nM by default.

    rprimerc : float, optional
        Concentration of reverse primer in nM, set to 1000.0 nM by default.

    saltc  : float, optional
        Salt concentration (monovalet cations) :mod:`tmbresluc` set to 50.0 mM by default

    formula : function
        formula used for tm calculation
        this is the name of a function.
        built in options are:

        1. :func:`pydna.amplify.tmbresluc` (default)
        2. :func:`pydna.amplify.basictm`
        3. :func:`pydna.amplify.tmstaluc98`
        4. :func:`pydna.amplify.tmbreslauer86`

        These functions are imported from the :mod:`pydna.amplify` module, but can be
        substituted for some other custom made function.

    Returns
    -------
    result : Amplicon

    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> t=Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    >>> t
    Dseqrecord(-64)
    >>> from pydna.design import primer_design
    >>> ampl = primer_design(t)
    >>> ampl
    Amplicon(64)
    >>> ampl.forward_primer
    fw64 18-mer:5'-atgactgctaacccttcc-3'
    >>> ampl.reverse_primer
    rv64 19-mer:5'-catcgtaagtttcgaacga-3'
    >>> print(ampl.figure())
    5atgactgctaacccttcc...tcgttcgaaacttacgatg3
                          ||||||||||||||||||| tm 53.8 (dbd) 60.6
                         3agcaagctttgaatgctac5
    5atgactgctaacccttcc3
     |||||||||||||||||| tm 54.4 (dbd) 58.4
    3tactgacgattgggaagg...agcaagctttgaatgctac5
    >>> pf = "GGATCC" + ampl.forward_primer
    >>> pr = "GGATCC" + ampl.reverse_primer  
    >>> pf
    fw64 24-mer:5'-GGATCCatgactgct..tcc-3'
    >>> pr
    rv64 25-mer:5'-GGATCCcatcgtaag..cga-3'
    >>> from pydna.amplify import pcr
    >>> pcr_prod = pcr(pf, pr, t)
    >>> print(pcr_prod.figure())
          5atgactgctaacccttcc...tcgttcgaaacttacgatg3
                                ||||||||||||||||||| tm 53.8 (dbd) 60.6
                               3agcaagctttgaatgctacCCTAGG5
    5GGATCCatgactgctaacccttcc3
           |||||||||||||||||| tm 54.4 (dbd) 58.4
          3tactgacgattgggaagg...agcaagctttgaatgctac5
    >>> print(pcr_prod.seq)
    GGATCCatgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatgGGATCC
    >>> from pydna.primer import Primer
    >>> pf = Primer("atgactgctaacccttccttggtgttg", id="myprimer")
    >>> ampl = primer_design(t, fp = pf)
    >>> ampl.forward_primer
    myprimer 27-mer:5'-atgactgctaaccct..ttg-3'
    >>> ampl.reverse_primer
    rv64 28-mer:5'-catcgtaagtttcga..gtc-3'

    '''
    
    def design(target_tm, template):
        ''' returns a string '''
        tmp=0
        length=0
        while tmp<target_tm:
            length+=1
            p = str(template.seq[:length])
            tmp = formula(p.upper())
        ps = p[:-1]
        tmps = formula(str(ps).upper())
        _module_logger.debug(((p,   tmp),(ps, tmps)))
        return min( ( abs(target_tm-tmp), p), (abs(target_tm-tmps), ps) )[1]
    
    if fp and not rp:
        fp  = _Anneal((fp,), template).forward_primers.pop()
        target_tm = formula( str(fp.footprint), primerc=fprimerc, saltc=saltc)
        _module_logger.debug("forward primer given, design reverse primer:")
        rp = _Primer(design(target_tm, template.rc()))
    elif not fp and rp:
        rp =  _Anneal([rp], template).reverse_primers.pop()
        target_tm = formula( str(rp.footprint), primerc=rprimerc, saltc=saltc)
        _module_logger.debug("reverse primer given, design forward primer:")
        fp = _Primer(design(target_tm, template))
    elif not fp and not rp:
        _module_logger.debug("no primer given, design forward primer:")
        fp = _Primer((design(target_tm, template)))
        target_tm = formula( str(fp.seq), primerc=fprimerc, saltc=saltc)
        _module_logger.debug("no primer given, design reverse primer:")
        rp = _Primer(design(target_tm, template.rc()))
    else:
        raise Exception("Specify maximum one of the two primers.")

    ampl = _Anneal( (fp, rp), template)
    
    prod = ampl.products[0]
    
    prod.forward_primer.concentration = fprimerc
    prod.reverse_primer.concentration = rprimerc

    if prod.forward_primer.id == "<unknown id>":
        prod.forward_primer.id = "fw{}".format(len(template))
        
    if prod.reverse_primer.id == "<unknown id>":
        prod.reverse_primer.id = "rv{}".format(len(template))

    if prod.forward_primer.name == "<unknown name>":
        prod.forward_primer.name = "fw{}".format(len(template))
        
    if prod.reverse_primer.name == "<unknown name>":
        prod.reverse_primer.name = "rv{}".format(len(template))

    prod.forward_primer.description = prod.forward_primer.id+' '+template.accession
    prod.reverse_primer.description = prod.reverse_primer.id+' '+template.accession

    return prod


def integration_primers( up,
                         cas,
                         dn,
                         uplink     = _Dseqrecord(''),
                         dnlink     = _Dseqrecord(''),
                         minlength  = 16,
                         maxlength  = 80,
                         min_olap   = 50,
                         target_tm  = 55.0,
                         fprimerc   = 1000.0,
                         rprimerc   = 1000.0,
                         saltc      = 50.0,
                         formula    = _tmbresluc):

    fp_tail = str(up[-min_olap:].seq) + str(uplink.seq)
    rp_tail = str(dn[:min_olap].rc().seq) + str(dnlink.rc().seq)
    from pydna import _PydnaDeprecationWarning
    warnings.warn("This function has been deprecated."
                  "consider using the primer_design"
                  "function instead", _PydnaDeprecationWarning)
    
    '''
    **Do not use this function, use pydna.design.assembly_fragments instead**
    **This function will be deprecated and removed in a future version of pydna**
    **This can be discussed in the google group https://groups.google.com/forum/#!forum/pydna ***
    '''
    return cloning_primers( cas,
                            minlength=minlength,
                            maxlength=maxlength,
                            fp_tail=fp_tail,
                            rp_tail=rp_tail,
                            target_tm=target_tm,
                            fprimerc=fprimerc,
                            rprimerc=rprimerc,
                            saltc=saltc,
                            formula = formula)


def assembly_primers(templates,
                     circular   = False,
                     vector     = None,
                     maxlink    = 40,
                     minlength  = 16,
                     maxlength  = 50,
                     min_olap   = 35,
                     target_tm  = 55.0,
                     fprimerc   = 1000.0,
                     rprimerc   = 1000.0,
                     saltc      = 50.0,
                     formula    = _tmbresluc):
    from pydna import _PydnaDeprecationWarning
    warnings.warn("This function has been deprecated."
              "consider using the primer_design"
              "function instead", _PydnaDeprecationWarning)


    '''**Do not use this function, use pydna.design.assembly_fragments instead**
    **This function will be deprecated and removed in a future version of pydna**
    **This can be discussed in the google group https://groups.google.com/forum/#!forum/pydna *** 
    
    This function return primer pairs that are useful for fusion of DNA sequences given in template.
    Given two sequences that we wish to fuse (a and b) to form fragment c.

    ::

       _________ a _________           __________ b ________
      /                     \\         /                     \\
      agcctatcatcttggtctctgca   <-->  TTTATATCGCATGACTCTTCTTT
      |||||||||||||||||||||||         |||||||||||||||||||||||
      tcggatagtagaaccagagacgt   <-->  AAATATAGCGTACTGAGAAGAAA


           agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
           ||||||||||||||||||||||||||||||||||||||||||||||
           tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA
           \\___________________ c ______________________/


    We can design tailed primers to fuse a and b by fusion PCR, Gibson assembly or
    in-vivo homologous recombination. The basic requirements for the primers for
    the three techniques are the same.


    Design tailed primers incorporating a part of the next or previous fragment to be assembled.

    ::


      agcctatcatcttggtctctgca
      |||||||||||||||||||||||
                      gagacgtAAATATA

      |||||||||||||||||||||||
      tcggatagtagaaccagagacgt


                             TTTATATCGCATGACTCTTCTTT
                             |||||||||||||||||||||||

                      ctctgcaTTTATAT
                             |||||||||||||||||||||||
                             AAATATAGCGTACTGAGAAGAAA

    PCR products with flanking sequences are formed in the PCR process.

    ::

      agcctatcatcttggtctctgcaTTTATAT
      ||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATA
                      \\____________/

                         identical
                         sequences
                       ____________
                      /            \\
                      ctctgcaTTTATATCGCATGACTCTTCTTT
                      ||||||||||||||||||||||||||||||
                      gagacgtAAATATAGCGTACTGAGAAGAAA

    The fragments can be fused by any of the techniques mentioned earlier.

    ::

      agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
      ||||||||||||||||||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA




    Parameters
    ----------

    templates : list of Dseqrecord
        list Dseqrecord object for which fusion primers should be constructed.

    minlength : int, optional
        Minimum length of the annealing part of the primer.

    maxlength : int, optional
        Maximum length (including tail) for designed primers.

    tot_length : int, optional
        Maximum total length of a the primers

    target_tm : float, optional
        target tm for the primers

    primerc : float, optional
        Concentration of each primer in nM, set to 1000.0 nM by default

    saltc  : float, optional
        Salt concentration (monovalet cations) :mod:`tmbresluc` set to 50.0 mM by default

    formula : function
        formula used for tm calculation
        this is the name of a function.
        built in options are:

        1. :func:`pydna.tm.tmbresluc` (default)
        2. :func:`pydna.tm.basictm`
        3. :func:`pydna.tm.tmstaluc98`
        4. :func:`pydna.tm.tmbreslauer86`

        These functions are imported from the :mod:`pydna.tm` module, but can be
        substituted for some other custom made function.

    Returns
    -------
    primer_pairs : list of tuples of :mod:`Bio.Seqrecord` objects

        ::

          [(forward_primer_1, reverse_primer_1),
           (forward_primer_2, reverse_primer_2), ...]


    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.assembly import Assembly
    >>> from pydna.design import assembly_primers
    >>> from pydna.amplify import pcr
    >>> a=Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    >>> b=Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc")
    >>> c=Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt")
    >>> primer_pairs = assembly_primers([a,b,c], circular = True)
    >>> p=[]
    >>> for t, (f,r) in zip([a,b,c], primer_pairs): p.append(pcr(f,r,t))
    >>> p
    [Amplicon(100), Amplicon(101), Amplicon(102)]
    >>> assemblyobj = Assembly(p)
    >>> assemblyobj
    Assembly:
    Sequences........................: [100] [101] [102]
    Sequences with shared homologies.: [100] [101] [102]
    Homology limit (bp)..............: 25
    Number of overlaps...............: 3
    Nodes in graph(incl. 5' & 3')....: 5
    Only terminal overlaps...........: No
    Circular products................: [195]
    Linear products..................: [231] [231] [231] [167] [166] [165] [36] [36] [36]
    >>> assemblyobj.linear_products
    [Contig(-231), Contig(-231), Contig(-231), Contig(-167), Contig(-166), Contig(-165), Contig(-36), Contig(-36), Contig(-36)]
    >>> assemblyobj.circular_products[0].cseguid()
    'V3Mi8zilejgyoH833UbjJOtDMbc'
    >>> (a+b+c).looped().cseguid()
    'V3Mi8zilejgyoH833UbjJOtDMbc'
    >>> print(assemblyobj.circular_products[0].small_fig())
     -|100bp_PCR_prod|36
    |                 \\/
    |                 /\\
    |                 36|101bp_PCR_prod|36
    |                                   \\/
    |                                   /\\
    |                                   36|102bp_PCR_prod|36
    |                                                     \\/
    |                                                     /\\
    |                                                     36-
    |                                                        |
     --------------------------------------------------------
    >>>





    '''

    if vector:
        circular = False

    if not hasattr(templates, '__iter__'):
        raise Exception("first argument has to be an iterable")

    tail_length =  int(_math.ceil(float(min_olap)/2))

    if circular:
        templates = tuple(templates) + (templates[0],)

    linkers = []
    newtemplates=[]

    for t in templates:
        if len(t)<=maxlink:
            linker = t
        else:
            linker= _Dseqrecord('')
            newtemplates.append(t)
        linkers.append(linker)

    newtails = []
    linkers = linkers[1:]+[linkers[0]]

    for t,n,l in zip(newtemplates, newtemplates[1:], linkers):
        length = tail_length - len(l)//2
        newtails.extend( (str(n[:length].rc().seq)+str(l.rc().seq),
                          str(t.seq[-length:])+str(l.seq) ))

    if not circular:
        newtails = [""] + newtails + [""]
    else:
        newtails = [newtails[-1]] + newtails[:-1]

    tails =  list(zip(newtails[::2], newtails[1::2]))

    primer_pairs = []

    for template, (fp_tail, rp_tail) in zip(newtemplates, tails):

        fp, rp = cloning_primers(   template,
                                    minlength  = minlength,
                                    maxlength  = maxlength,
                                    fp_tail    = fp_tail,
                                    rp_tail    = rp_tail,
                                    target_tm  = target_tm,
                                    fprimerc   = fprimerc,
                                    rprimerc   = rprimerc,
                                    saltc      = saltc,
                                    formula    = formula)

        #fp = _Primer(_Seq(fp_tail, IUPACAmbiguousDNA())) + fp
        #rp = _Primer(_Seq(fp_tail, IUPACAmbiguousDNA())) + rp

        primer_pairs.append((fp, rp))

    if vector:
        fake_cas = templates[0]+templates[-1]
        pa, pb   = integration_primers( vector,
                                        fake_cas,
                                        vector,
                                        min_olap=min_olap)

        primer_pairs[0]  = (pa, primer_pairs[0][1])
        primer_pairs[-1] = (primer_pairs[-1][0], pb)

    return primer_pairs

def assembly_fragments(f, overlap=35, maxlink=40):
    
    '''This function return a list of :mod:`pydna.amplicon.Amplicon` objects where 
    primers have been modified with tails so that the fragments can be fused in 
    the order they appear in the list by for example Gibson assembly or homologous 
    recombination.
    
    Given that we have two linear :mod:`pydna.amplicon.Amplicon` objects a and b 
    
    we can modify the reverse primer of a and forward primer of b with tails to allow 
    fusion fusion PCR, Gibson assembly or in-vivo homologous recombination.
    The basic requirements for the primers for the three techniques are the same.

    ::
        
                                <-->

       _________ a _________           __________ b ________
      /                     \\         /                     \\
      agcctatcatcttggtctctgca         TTTATATCGCATGACTCTTCTTT
      |||||||||||||||||||||||         |||||||||||||||||||||||
                       <gacgt                          <AGAAA
      agcct>                          TTTAT>
      |||||||||||||||||||||||         |||||||||||||||||||||||
      tcggatagtagaaccagagacgt         AAATATAGCGTACTGAGAAGAAA


           agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
           ||||||||||||||||||||||||||||||||||||||||||||||
           tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA
           \\___________________ c ______________________/


    Design tailed primers incorporating a part of the next or previous fragment to be assembled.

    ::


      agcctatcatcttggtctctgca
      |||||||||||||||||||||||
                      gagacgtAAATATA

      |||||||||||||||||||||||
      tcggatagtagaaccagagacgt


                             TTTATATCGCATGACTCTTCTTT
                             |||||||||||||||||||||||

                      ctctgcaTTTATAT
                             |||||||||||||||||||||||
                             AAATATAGCGTACTGAGAAGAAA

    PCR products with flanking sequences are formed in the PCR process.

    ::

      agcctatcatcttggtctctgcaTTTATAT
      ||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATA
                      \\____________/

                         identical
                         sequences
                       ____________
                      /            \\
                      ctctgcaTTTATATCGCATGACTCTTCTTT
                      ||||||||||||||||||||||||||||||
                      gagacgtAAATATAGCGTACTGAGAAGAAA

    The fragments can be fused by any of the techniques mentioned earlier to form c:

    ::

      agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
      ||||||||||||||||||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA


    The first argument of this function is a list of sequence objects containing 
    Amplicons and other similar objects.
    
    **At least every second sequence object needs to be an Amplicon**
    
    This rule exists because if a sequence object is that is not a PCR product
    is to be fused with another fragment, that other fragment needs to be an Amplicon
    so that the primer of the other object can be modified to include the whole stretch
    of sequence homology needed for the fusion. See the example below where a is a 
    non-amplicon (a linear plasmid  vector for instance)

    ::

       _________ a _________           __________ b ________
      /                     \\         /                     \\
      agcctatcatcttggtctctgca   <-->  TTTATATCGCATGACTCTTCTTT
      |||||||||||||||||||||||         |||||||||||||||||||||||
      tcggatagtagaaccagagacgt                          <AGAAA
                                      TTTAT>
                                      |||||||||||||||||||||||
                                <-->  AAATATAGCGTACTGAGAAGAAA


           agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
           ||||||||||||||||||||||||||||||||||||||||||||||
           tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA
           \\___________________ c ______________________/


    In this case only the forward primer of b is fitted with a tail with a part a:

    ::


      agcctatcatcttggtctctgca
      |||||||||||||||||||||||
      tcggatagtagaaccagagacgt


                             TTTATATCGCATGACTCTTCTTT
                             |||||||||||||||||||||||
                                              <AGAAA
               tcttggtctctgcaTTTATAT
                             |||||||||||||||||||||||
                             AAATATAGCGTACTGAGAAGAAA

    PCR products with flanking sequences are formed in the PCR process.

    ::

      agcctatcatcttggtctctgcaTTTATAT
      ||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATA
                      \\____________/

                         identical
                         sequences
                       ____________
                      /            \\
                      ctctgcaTTTATATCGCATGACTCTTCTTT
                      ||||||||||||||||||||||||||||||
                      gagacgtAAATATAGCGTACTGAGAAGAAA

    The fragments can be fused by for example Gibson assembly:

    ::

      agcctatcatcttggtctctgcaTTTATAT
      ||||||||||||||||||||||||||||||
      tcggatagtagaacca
                      
                                   TCGCATGACTCTTCTTT
                      ||||||||||||||||||||||||||||||
                      gagacgtAAATATAGCGTACTGAGAAGAAA 
                      
    to form c:

    ::
                
      agcctatcatcttggtctctgcaTTTATATCGCATGACTCTTCTTT
      ||||||||||||||||||||||||||||||||||||||||||||||
      tcggatagtagaaccagagacgtAAATATAGCGTACTGAGAAGAAA
      

    The first argument of this function is a list of sequence objects containing 
    Amplicons and other similar objects.
    
    The overlap argument controls how many base pairs of overlap required between 
    adjacent sequence fragments. In the junction between Amplicons, tails with the 
    length of about half of this value is added to the two primers
    closest to the junction.
    
    ::
        
            >       <
            Amplicon1
                     Amplicon2
                     >       <
                     
                     ⇣

            >       <-
            Amplicon1
                     Amplicon2
                    ->       <                     
                     
    In the case of an Amplicon adjacent to a Dseqrecord object, the tail will 
    be twice as long (1*overlap) since the 
    recombining sequence is present entirely on this primer:
        
    ::
        
            Dseqrecd1
                     Amplicon1
                     >       <
                     
                     ⇣

            Dseqrecd1
                     Amplicon1
                   -->       <
    
    Note that if the sequence of DNA fragments starts or stops with an Amplicon, 
    the very first and very last prinmer will not be modified i.e. assembles are 
    always assumed to be linear. There are simple tricks around that for circular
    assemblies depicted in the last two examples below.
    
    The maxlink arguments controls the cut off length for sequences that will be
    synhtesized by adding them to primers for the adjacent fragment(s). The 
    argument list may contain short spacers (such as spacers between fusion proteins).
    

    ::

        Example 1: Linear assembly of PCR products (pydna.amplicon.Amplicon class objects) ------
        
        
        >       <         >       <
        Amplicon1         Amplicon3
                 Amplicon2         Amplicon4
                 >       <         >       <
        
                             ⇣
                             pydna.design.assembly_fragments
                             ⇣ 
        
        >       <-       ->       <-                       pydna.assembly.Assembly
        Amplicon1         Amplicon3                         
                 Amplicon2         Amplicon4     ➤  Amplicon1Amplicon2Amplicon3Amplicon4
                ->       <-       ->       <
        
        
        Example 2: Linear assembly of alternating Amplicons and other fragments
        
        
        >       <         >       <
        Amplicon1         Amplicon2
                 Dseqrecd1         Dseqrecd2
                      
                             ⇣
                             pydna.design.assembly_fragments
                             ⇣ 
                          
        >       <--     -->       <--                     pydna.assembly.Assembly
        Amplicon1         Amplicon2
                 Dseqrecd1         Dseqrecd2     ➤  Amplicon1Dseqrecd1Amplicon2Dseqrecd2
        
        
        Example 3: Linear assembly of alternating Amplicons and other fragments
        
        
        Dseqrecd1         Dseqrecd2
                 Amplicon1         Amplicon2
                 >       <       -->       <
        
                             ⇣
                     pydna.design.assembly_fragments
                             ⇣
                                                          pydna.assembly.Assembly
        Dseqrecd1         Dseqrecd2
                 Amplicon1         Amplicon2     ➤  Dseqrecd1Amplicon1Dseqrecd2Amplicon2
               -->       <--     -->       <
        
        
        Example 4: Circular assembly of alternating Amplicons and other fragments
        
                         ->       <==
        Dseqrecd1         Amplicon2
                 Amplicon1         Dseqrecd1
               -->       <-
                             ⇣
                             pydna.design.assembly_fragments
                             ⇣ 
                                                           pydna.assembly.Assembly
                         ->       <==
        Dseqrecd1         Amplicon2                    -Dseqrecd1Amplicon1Amplicon2-  
                 Amplicon1                       ➤    |                             |
               -->       <-                            -----------------------------
        
        ------ Example 5: Circular assembly of Amplicons
        
        >       <         >       <
        Amplicon1         Amplicon3
                 Amplicon2         Amplicon1
                 >       <         >       <
        
                             ⇣
                             pydna.design.assembly_fragments
                             ⇣ 
        
        >       <=       ->       <-        
        Amplicon1         Amplicon3                  
                 Amplicon2         Amplicon1
                ->       <-       +>       <
        
                             ⇣
                     make new Amplicon using the Amplicon1.template and 
                     the last fwd primer and the first rev primer.
                             ⇣
                                                           pydna.assembly.Assembly
        +>       <=       ->       <-        
         Amplicon1         Amplicon3                  -Amplicon1Amplicon2Amplicon3-
                  Amplicon2                      ➤   |                             |
                 ->       <-                          -----------------------------
        
        


    Parameters
    ----------

    f : list of :mod:`pydna.amplicon.Amplicon` and other Dseqrecord like objects
        list Amplicon and Dseqrecord object for which fusion primers should be constructed.

    overlap : int, optional
        Minimum length of the annealing part of the primer.

    maxlink : int, optional
        Maximum length of spacers that will be included in tails for designed primers.

    Returns
    -------
    seqs : list of :mod:`pydna.amplicon.Amplicon` and other Dseqrecord like objects :mod:`pydna.amplicon.Amplicon` objects

        ::

          [Amplicon1,
           Amplicon2, ...]


    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord    
    >>> from pydna.design import primer_design
    >>> a=primer_design(Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg"))
    >>> b=primer_design(Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc"))
    >>> c=primer_design(Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt"))
    >>> from pydna.design import assembly_fragments
    >>> # We would like a circular recombination, so the first sequence has to be repeated
    >>> fa1,fb,fc,fa2 = assembly_fragments([a,b,c,a])
    >>> # Since all fragments are Amplicons, we need to extract the rp of the 1st and fp of the last fragments.
    >>> from pydna.amplify import pcr
    >>> fa = pcr(fa2.forward_primer, fa1.reverse_primer, a)
    >>> [fa,fb,fc]
    [Amplicon(100), Amplicon(101), Amplicon(102)]
    >>> from pydna.assembly import Assembly
    >>> assemblyobj = Assembly([fa,fb,fc])
    >>> assemblyobj
    Assembly:
    Sequences........................: [100] [101] [102]
    Sequences with shared homologies.: [100] [101] [102]
    Homology limit (bp)..............: 25
    Number of overlaps...............: 3
    Nodes in graph(incl. 5' & 3')....: 5
    Only terminal overlaps...........: No
    Circular products................: [195]
    Linear products..................: [231] [231] [231] [167] [166] [165] [36] [36] [36]
    >>> assemblyobj.linear_products
    [Contig(-231), Contig(-231), Contig(-231), Contig(-167), Contig(-166), Contig(-165), Contig(-36), Contig(-36), Contig(-36)]
    >>> assemblyobj.circular_products[0].cseguid()
    'V3Mi8zilejgyoH833UbjJOtDMbc'
    >>> (a+b+c).looped().cseguid()
    'V3Mi8zilejgyoH833UbjJOtDMbc'
    >>> print(assemblyobj.circular_products[0].small_fig())
     -|100bp_PCR_prod|36
    |                 \\/
    |                 /\\
    |                 36|101bp_PCR_prod|36
    |                                   \\/
    |                                   /\\
    |                                   36|102bp_PCR_prod|36
    |                                                     \\/
    |                                                     /\\
    |                                                     36-
    |                                                        |
     --------------------------------------------------------
    >>>

    '''
    # sanity check for arguments
    nf = [item for item in f if len(item)>maxlink]
    if not all(hasattr(i[0],"template") or hasattr(i[1],"template") for i in zip(nf,nf[1:])):
        raise Exception("Every second fragment larger than maxlink has to be an Amplicon object")
    
    tail_length = _math.ceil(overlap/2)
    f = [_copy.copy(f) for f in f]    
    if len(f[0])<=maxlink:
        f[1].forward_primer = f[0].seq.todata + f[1].forward_primer
        f=f[1:]
    if len(f[-1])<=maxlink:
        f[-2].reverse_primer = f[-1].seq.rc().todata + f[-2].reverse_primer
        f=f[:-1]    
    empty = _Dseqrecord("")    
    for i in range(len(f)-1):                                  
        if len(f[i+1])<=maxlink:  # f[i+1] is smaller than maxlink
            if hasattr(f[i], "template") and hasattr(f[i+2], "template"):
                lnk = str(f[i].seq[(-tail_length+len(f[i+1])//2 if len(f[i+1])//2<overlap else len(f[i+1])): ])
                f[i+2].forward_primer = lnk + f[i+1].seq.todata + f[i+2].forward_primer
                lnk = str(f[i+2].seq[ :(tail_length-len(f[i+1])//2 if len(f[i+1])//2<overlap else 0)].rc())
                f[i].reverse_primer = lnk + f[i+1].seq.rc().todata + f[i].reverse_primer
            elif hasattr(f[i] , "template"):
                lnk = str(f[i+2].seq[:overlap].rc())
                f[i].reverse_primer = lnk + f[i+1].seq.rc().todata + f[i].reverse_primer
            elif hasattr(f[i+2] , "template"):
               lnk = str(f[i].seq[-overlap:])
               f[i+2].forward_primer = lnk + f[i+1].seq.todata + f[i+2].forward_primer
            f[i+1]=empty
        else:                    # f[i+1] is larger than maxlink
            if hasattr(f[i], "template") and hasattr(f[i+1], "template"):
                lnk = str(f[i].seq[-tail_length:])            
                f[i+1].forward_primer = lnk + f[i+1].forward_primer
                lnk = str(f[i+1].seq[:tail_length].rc())
                f[i].reverse_primer = lnk + f[i].reverse_primer            
            elif hasattr(f[i] , "template"):
                lnk = str(f[i+1].seq[:overlap].rc())
                f[i].reverse_primer = lnk + f[i].reverse_primer                
            elif hasattr(f[i+1] , "template"):
                lnk = str(f[i].seq[-overlap:])
                f[i+1].forward_primer = lnk + f[i+1].forward_primer
 
    f = [item for item in f if len(item)]
    return [_pcr(p.forward_primer,p.reverse_primer,p.template) if hasattr(p, "template") else p for p in f]

if __name__=="__main__":
    import os as _os
    cache = _os.getenv("pydna_cache", "nocache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
