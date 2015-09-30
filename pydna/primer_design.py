#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2013, 2014 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''
This module contain functions for primer design.

'''



import math
from operator import itemgetter
from Bio.Alphabet import Alphabet
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydna.amplify import Anneal, tmbresluc, Dseqrecord
from pydna.dsdna import parse
from pretty import pretty_str

class Primer(SeqRecord):

    def __init__(self, seq = Seq(""), *args, **kwargs):
        if seq.alphabet == Alphabet():
            seq.alphabet = IUPACAmbiguousDNA()
        super(Primer, self).__init__(seq, *args, **kwargs)

    def __repr__(self):
        return self.__class__.__name__ \
         + "(seq=%s, id=%s, name=%s, description=%s, dbxrefs=%s)" \
         % tuple(map(repr, (self.seq, self.id, self.name,
                            self.description, self.dbxrefs)))

    def __add__(self, other):
        new = super(Primer, self).__add__(other)
        return Primer(new.seq)

    def __radd__(self, other):
        new = super(Primer, self).__radd__(other)
        return Primer(new.seq)

def print_primer_pair(*args,**kwargs):
    f,r = cloning_primers(*args,**kwargs)
    return pretty_str("\n"+f.format("fasta")+"\n"+r.format("fasta") + "\n")

def cloning_primers(template,
                    minlength=16,
                    maxlength=29,
                    fp=None,
                    rp=None,
                    fp_tail='',
                    rp_tail='',
                    target_tm=55.0,
                    primerc = 1000.0,
                    saltc=50.0,
                    formula = tmbresluc,
                    path = u""):

    '''This function can design primers for PCR amplification of a given sequence.
    This function accepts a Dseqrecord object containing the template sequence and
    returns a tuple cntaining two ::mod`Bio.SeqRecord.SeqRecord` objects describing
    the primers.

    Primer tails can optionally be given in the form of strings.

    An predesigned primer can be given, either the forward or reverse primers. In this
    case this function tries to design a primer with a Tm to match the given primer.


    Parameters
    ----------

    template : Dseqrecord
        a Dseqrecord object.

    minlength : int, optional
        Minimum length of the annealing part of the primer

    maxlength : int, optional
        Maximum length (including tail) for designed primers.

    fp, rp : SeqRecord, optional
        optional Biopython SeqRecord objects containing one primer each.

    fp_tail, rp_tail : string, optional
        optional tails to be added to the forwars or reverse primers

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

        1. :func:`pydna.amplify.tmbresluc` (default)
        2. :func:`pydna.amplify.basictm`
        3. :func:`pydna.amplify.tmstaluc98`
        4. :func:`pydna.amplify.tmbreslauer86`

        These functions are imported from the :mod:`pydna.amplify` module, but can be
        substituted for some other custom made function.

    path : unicode, optional
        This variable can be set to a path to a text file, which will be created
        if it does not exist.
        This file (if it exists) will be parsed for sequences in fasta or
        genbank format and a Biopython SeqRecord object will be created for
        each sequence.

        If a SeqRecord object is found with the same description as any of the
        primers designed, the SeqRecord object parsed from the file will be
        returned by this function instead of the newly designed primer.

        If no sequence with the same description can be found, the primer(s)
        will be appended to the file in fasta format.


    Returns
    -------
    fp, rp : tuple
        fp is a :mod:Bio.SeqRecord object describing the forward primer
        rp is a :mod:Bio.SeqRecord object describing the reverse primer



    Examples
    --------

    >>> import pydna
    >>> t=pydna.Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    >>> t
    Dseqrecord(-64)
    >>> pf,pr = pydna.cloning_primers(t)
    >>> pf
    Primer(seq=Seq('atgactgctaacccttc', IUPACAmbiguousDNA()), id='fw64', name='fw64', description='fw64 -', dbxrefs=[])
    >>> pr
    Primer(seq=Seq('catcgtaagtttcgaac', IUPACAmbiguousDNA()), id='rv64', name='rv64', description='rv64 -', dbxrefs=[])
    >>> pcr_prod = pydna.pcr(pf, pr, t)
    >>> pcr_prod
    Amplicon(64)
    >>>
    >>> print(pcr_prod.figure())
    5atgactgctaacccttc...gttcgaaacttacgatg3
                         ||||||||||||||||| tm 42.4 (dbd) 52.9
                        3caagctttgaatgctac5
    5atgactgctaacccttc3
     ||||||||||||||||| tm 44.5 (dbd) 54.0
    3tactgacgattgggaag...caagctttgaatgctac5
    >>> pf,pr = pydna.cloning_primers(t, fp_tail="GGATCC", rp_tail="GAATTC")
    >>> pf
    Primer(seq=Seq('GGATCCatgactgctaacccttc', IUPACAmbiguousDNA()), id='fw64', name='fw64', description='fw64 -', dbxrefs=[])
    >>> pr
    Primer(seq=Seq('GAATTCcatcgtaagtttcgaac', IUPACAmbiguousDNA()), id='rv64', name='rv64', description='rv64 -', dbxrefs=[])
    >>> pcr_prod = pydna.pcr(pf, pr, t)
    >>> print(pcr_prod.figure())
          5atgactgctaacccttc...gttcgaaacttacgatg3
                               ||||||||||||||||| tm 42.4 (dbd) 52.9
                              3caagctttgaatgctacCTTAAG5
    5GGATCCatgactgctaacccttc3
           ||||||||||||||||| tm 44.5 (dbd) 54.0
          3tactgacgattgggaag...caagctttgaatgctac5
    >>> print(pcr_prod.seq)
    GGATCCatgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatgGAATTC
    >>>
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> pf = SeqRecord(Seq("atgactgctaacccttccttggtgttg"))
    >>> pf,pr = pydna.cloning_primers(t, fp = pf, fp_tail="GGATCC", rp_tail="GAATTC")
    >>> pf
    Primer(seq=Seq('GGATCCatgactgctaacccttccttggtgttg', Alphabet()), id='fw64', name='fw64', description='fw64 -', dbxrefs=[])
    >>> pr
    Primer(seq=Seq('GAATTCcatcgtaagtttcgaacgaaatgtcgtc', IUPACAmbiguousDNA()), id='rv64', name='rv64', description='rv64 -', dbxrefs=[])
    >>> ampl = pydna.pcr(pf,pr,t)
    >>> print(ampl.figure())
          5atgactgctaacccttccttggtgttg...gacgacatttcgttcgaaacttacgatg3
                                         |||||||||||||||||||||||||||| tm 57.5 (dbd) 72.2
                                        3ctgctgtaaagcaagctttgaatgctacCTTAAG5
    5GGATCCatgactgctaacccttccttggtgttg3
           ||||||||||||||||||||||||||| tm 59.0 (dbd) 72.3
          3tactgacgattgggaaggaaccacaac...ctgctgtaaagcaagctttgaatgctac5
    >>>


    '''

    if fp and not rp:
        fp = Primer(Seq(fp_tail, IUPACAmbiguousDNA())) + fp
        p  = Anneal([fp], template).fwd_primers.pop()
        fp = Primer(p.footprint)
        fp_tail = Primer(p.tail)
        rp = Primer(Seq(str(template[-(maxlength*3-len(rp_tail)):].reverse_complement().seq), IUPACAmbiguousDNA()))
        target_tm = formula(str(fp.seq).upper(), primerc=primerc, saltc=saltc)
    elif not fp and rp:
        rp = Primer(Seq(rp_tail, IUPACAmbiguousDNA())) + rp
        p =  Anneal([rp], template).rev_primers.pop()
        rp = Primer(p.footprint)
        rp_tail = Primer(p.tail)
        fp = Primer(Seq(str(template[:maxlength*3-len(fp_tail)].seq), IUPACAmbiguousDNA()))
        target_tm = formula(str(rp.seq).upper(), primerc=primerc, saltc=saltc)
    elif not fp and not rp:
        fp = Primer(Seq(str(template[:maxlength-len(fp_tail)].seq), IUPACAmbiguousDNA()))
        rp = Primer(Seq(str(template[-maxlength+len(rp_tail):].reverse_complement().seq), IUPACAmbiguousDNA()))
    else:
        raise Exception("Specify maximum one of the two primers, not both.")

    lowtm, hightm = sorted( [( formula(str(fp.seq), primerc, saltc), fp, "f" ),
                             ( formula(str(rp.seq), primerc, saltc), rp, "r" ) ] )

    while lowtm[0] > target_tm and len(lowtm[1])>minlength:
        shorter = lowtm[1][:-1]
        tm      = formula(str(shorter.seq).upper(), primerc=primerc, saltc=saltc)
        lowtm   = (tm, shorter, lowtm[2])

    while hightm[0] > lowtm[0] + 2.0 and len(hightm[1])>minlength:
        shorter = hightm[1][:-1]
        tm = formula(str(shorter.seq).upper(), primerc = primerc, saltc = saltc)
        hightm = (tm, shorter, hightm[2])

    fp, rp = sorted((lowtm, hightm), key=itemgetter(2))

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

    if fp.seq.alphabet == Alphabet():
        fp.seq.alphabet = IUPACAmbiguousDNA()
    if rp.seq.alphabet == Alphabet():
        rp.seq.alphabet = IUPACAmbiguousDNA()

#    If the path argument is supplied primers will be written to a file with that
#    path. If the file does not exist, it will be created and both primers will be
#    written to it. If the file exists, the file will be parsed for sequences in
#    fasta or genbank into .

    if path:
        try:
            with open(path, 'rU') as f: raw = f.read()
        except IOError:
            raw = u""
            with open(path, 'w') as f:
                f.write(fp.format("fasta"))
                f.write(rp.format("fasta"))
        else:
            primer_dict = {x.description:x for x in parse(raw, ds=False)}
            try:
                fp = primer_dict[fp.description]
            except KeyError:
                with open(path, 'a') as f: f.write(u"\n"+fp.format("fasta").strip())
            try:
                rp = primer_dict[rp.description]
            except KeyError:
                with open(path, 'a') as f: f.write(u"\n"+rp.format("fasta").strip())
    return fp, rp

def integration_primers( up,
                         cas,
                         dn,
                         uplink     = Dseqrecord(''),
                         dnlink     = Dseqrecord(''),
                         minlength  = 16,
                         min_olap   = 50,
                         target_tm  = 55.0,
                         primerc    = 1000.0,
                         saltc      = 50.0,
                         formula    = tmbresluc,
                         path       = u""):

    fp_tail = str(up[-min_olap:].seq) + str(uplink.seq)
    rp_tail = str(dn[:min_olap].rc().seq) + str(dnlink.rc().seq)

    maxlength  = minlength + max(len(fp_tail), len(rp_tail))

    return cloning_primers( cas,
                            minlength=minlength,
                            maxlength=maxlength,
                            fp_tail=fp_tail,
                            rp_tail=rp_tail,
                            target_tm=target_tm,
                            primerc = primerc,
                            saltc=saltc,
                            formula = formula,
                            path = path)


def assembly_primers(templates,
                     circular   = False,
                     vector     = None,
                     maxlink    = 40,
                     minlength  = 16,
                     maxlength  = 50,
                     min_olap   = 35,
                     target_tm  = 55.0,
                     primerc    = 1000.0,
                     saltc      = 50.0,
                     formula    = tmbresluc,
                     path       = u""):


    '''This function return primer pairs that are useful for fusion of DNA sequences given in template.
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

        1. :func:`pydna.amplify.tmbresluc` (default)
        2. :func:`pydna.amplify.basictm`
        3. :func:`pydna.amplify.tmstaluc98`
        4. :func:`pydna.amplify.tmbreslauer86`

        These functions are imported from the :mod:`pydna.amplify` module, but can be
        substituted for some other custom made function.

    Returns
    -------
    primer_pairs : list of tuples of :mod:`Bio.Seqrecord` objects

        ::

          [(forward_primer_1, reverse_primer_1),
           (forward_primer_2, reverse_primer_2), ...]


    Examples
    --------

    >>> import pydna
    >>> a=pydna.Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    >>> b=pydna.Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc")
    >>> c=pydna.Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt")
    >>> primer_pairs = pydna.assembly_primers([a,b,c], circular = True)
    >>> p=[]
    >>> for t, (f,r) in zip([a,b,c], primer_pairs): p.append(pydna.pcr(f,r,t))
    >>> p
    [Amplicon(100), Amplicon(101), Amplicon(102)]
    >>> assemblyobj = pydna.Assembly(p)
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
        raise Exception("argument has to be an iterable")

    tail_length =  int(math.ceil(float(min_olap)/2))

    if circular:
        templates = tuple(templates) + (templates[0],)

    linkers = []
    newtemplates=[]

    for t in templates:
        if len(t)<=maxlink:
            linker = t
        else:
            linker= Dseqrecord('')
            newtemplates.append(t)
        linkers.append(linker)

    newtails = []
    linkers = linkers[1:]+[linkers[0]]

    for t,n,l in zip(newtemplates, newtemplates[1:], linkers):
        length = tail_length - len(l)/2
        newtails.extend( (str(n[:length].rc().seq)+str(l.rc().seq),
                          str(t.seq[-length:])+str(l.seq) ))

    if not circular:
        newtails = [""] + newtails + [""]
    else:
        pass
        newtails = [newtails[-1]] + newtails[:-1]

    tails =  zip(newtails[::2], newtails[1::2])

    primer_pairs = []

    for template, (fp_tail, rp_tail) in zip(newtemplates, tails):

        fp, rp = cloning_primers(   template,
                                    minlength  = minlength,
                                    maxlength  = maxlength,
                                    fp_tail = fp_tail,
                                    rp_tail = rp_tail,
                                    target_tm  = target_tm,
                                    primerc    = primerc,
                                    saltc      = saltc,
                                    formula    = formula,
                                    path       = path)

        #fp = Primer(Seq(fp_tail, IUPACAmbiguousDNA())) + fp
        #rp = Primer(Seq(fp_tail, IUPACAmbiguousDNA())) + rp

        primer_pairs.append((fp, rp))

    if vector:
        fake_cas = templates[0]+templates[-1]
        pa, pb = integration_primers( vector,fake_cas,vector, min_olap=min_olap)
        primer_pairs[0]  = (pa, primer_pairs[0][1])
        primer_pairs[-1] = (primer_pairs[-1][0], pb)

    return primer_pairs


if __name__=="__main__":
    import doctest
    doctest.testmod()