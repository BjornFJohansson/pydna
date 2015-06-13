#!/usr/bin/env python
# -*- coding: utf-8 -*-
# doctest: +NORMALIZE_WHITESPACE
# doctest: +SKIP
'''This module provides functions for PCR. Primers with 5' tails as well as inverse PCR on
circular templates are handled correctly.

'''

import cPickle
import shelve

import itertools
import math
import re
import textwrap
import copy
import operator
import os


from math                           import log10
from Bio.Seq                        import Seq
from Bio.Alphabet.IUPAC             import ambiguous_dna
from Bio.SeqRecord                  import SeqRecord
from Bio.SeqUtils                   import GC
from Bio.SeqUtils.CheckSum          import seguid
from Bio.SeqUtils.MeltingTemp       import Tm_staluc
from Bio.SeqFeature                 import SeqFeature
from Bio.SeqFeature                 import CompoundLocation
from Bio.SeqFeature                 import FeatureLocation
from pydna.dsdna                    import rc
from pydna.dsdna                    import Dseqrecord
from pydna.pretty                   import pretty_str, pretty_unicode


def _annealing_positions(primer, template, limit=15):
    '''Finds the annealing position(s) for a primer on a template where the
    primer anneals perfectly with at least limit nucleotides in the 3' part.

    start is a position (integer)
    footprint1 and tail1 are strings.

    ::

        <- - - - - - - - - - template - - - - - - - - - - - - - >

        <------------- start ---->
     5'-...gctactacacacgtactgactgcctccaagatagagtcagtaaccacactcgat...3'
           ||||||||||||||||||||||||||||||||||||||||||||||||
                                  3'-gttctatctcagtcattggtgtATAGTG-5'

                                                        <tail>
                                     <---footprint----->
                                     <--------- primer ------>

    Parameters
    ----------
    primer : string
        The primer sequence 5'-3'

    template : string
        The template sequence 5'-3'

    limit : int = 15, optional
        footprint needs to be at least of length limit.

    Returns
    -------
    describe : list of tuples (int, string, string)
        [ (start1, footprint1, tail1), (start2, footprint2, tail2),..., ]
    '''

    if len(primer)<limit:
        return []
    prc = rc(primer)
    head = prc[:limit].upper()

    table = {"R":"(A|G)",
             "Y":"(C|T)",
             "S":"(G|C)",
             "W":"(A|T)",
             "K":"(G|T)",
             "M":"(A|C)",
             "B":"(C|G|T)",
             "D":"(A|G|T)",
             "H":"(A|C|T)",
             "V":"(A|C|G)",
             "N":"(A|G|C|T)"}

    for key in table:
        head=head.replace(key, table[key])

    positions = [m.start() for m in re.finditer('(?={})'.format(head), template, re.I)]

    if positions:
        tail = prc[limit:]
        length = len(tail)
        results = []
        for match_start in positions:
            tm = template[match_start+limit:match_start+limit+length]
            footprint = rc(template[match_start:match_start+limit]+"".join([b for a,b in itertools.takewhile(lambda x: x[0].lower()==x[1].lower(), zip(tail, tm))]))
            results.append((match_start, footprint, primer[: len(primer) - len(footprint) ]))
        return results
    return []

class Primer(SeqRecord):
    '''This class holds information about a primer with on a template '''

    def __init__(self,
                 seq_obj,
                 position=None,
                 footprint=None,
                 tail=None):

        self.position  = position
        self.footprint = footprint
        self.tail      = tail

        seq_obj.seq.alphabet = ambiguous_dna

        SeqRecord.__init__(self,
                           seq                = seq_obj.seq,
                           id                 = seq_obj.id,
                           name               = seq_obj.name,
                           description        = seq_obj.description,
                           dbxrefs            = seq_obj.dbxrefs,
                           features           = seq_obj.features,
                           annotations        = seq_obj.annotations,
                           letter_annotations = seq_obj.letter_annotations)


class Amplicon(Dseqrecord):
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
                     template=None,
                     forward_primer=None,
                     reverse_primer=None,
                     saltc=None,
                     forward_primer_concentration=None,
                     reverse_primer_concentration=None,
                     *args,
                     **kwargs):

        #Dseqrecord.__init__(self,record,*args,**kwargs)
        super(Amplicon, self).__init__(record, *args, **kwargs)
        self.template = template
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.fwd_primer = self.forward_primer
        self.rev_primer = self.reverse_primer
        self.forward_primer_concentration = forward_primer_concentration
        self.reverse_primer_concentration = reverse_primer_concentration
        self.saltc = saltc

        self.tmf = Tm_staluc(str(self.forward_primer.footprint),dnac=50, saltc=self.saltc)
        self.tmr = Tm_staluc(str(self.reverse_primer.footprint),dnac=50, saltc=self.saltc)
        self.tmf_dbd = tmbresluc(str(self.forward_primer.footprint),primerc=self.forward_primer_concentration)
        self.tmr_dbd = tmbresluc(str(self.reverse_primer.footprint),primerc=self.reverse_primer_concentration)


    def __getitem__(self, sl):
        answer = copy.copy(self)
        answer.seq = answer.seq.__getitem__(sl)
        answer.seq.alphabet = self.seq.alphabet
        sr = SeqRecord("n"*len(self))
        sr.features = self.features
        answer.features = SeqRecord.__getitem__(sr, sl).features
        return answer

    def __repr__(self):
        '''returns a short string representation of the object'''
        return "Amplicon({})".format(self.__len__())

    def flankup(self, flankuplength=50):
        '''Returns a Dseqrecord object containing flankuplength bases upstream of the forward primer footprint,
       Truncated if the template is not long enough.

       ::

        <--- flankup --->

                  5TAATAAactactgactatct3
                         ||||||||||||||
        acgcattcagctactgtactactgactatctatcg

       '''
        return self.template.seq[self.forward_primer.position-flankuplength-len(self.forward_primer.footprint):self.forward_primer.position-len(self.forward_primer.footprint)]

    def flankdn(self, flankdnlength=50):
        '''Returns a Dseqrecord object containing flankdnlength bases downstream of the reverse primer footprint.
       Truncated if the template is not long enough.

       ::

                                       <---- flankdn ------>

                        3actactgactatctTAATAA5
                         ||||||||||||||
        acgcattcagctactgtactactgactatctatcgtacatgtactatcgtat


       '''
        return self.template.seq[self.reverse_primer.position+len(self.reverse_primer.footprint):self.reverse_primer.position+flankdnlength+len(self.reverse_primer.footprint)]



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
                        tmf      = round(self.tmf,1),
                        tmr      = round(self.tmr,1),
                        tmf_dbd  = round(self.tmf_dbd,1),
                        tmr_dbd  = round(self.tmr_dbd,1),
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
        return pretty_str(textwrap.dedent(f).strip("\n"))

    def program(self):
        '''Returns a string containing a text representation of a proposed
       PCR program using Taq or similar polymerase.

       ::

        Taq (rate 30 nt/s)
        Three-step|         30 cycles     |      |SantaLucia 1998
        94.0°C    |94.0°C                 |      |SaltC 50mM
        __________|_____          72.0°C  |72.0°C|
        04min00s  |30s  \         ________|______|
                  |      \ 46.0°C/ 0min 1s|10min |
                  |       \_____/         |      |
                  |         30s           |      |4-8°C

       '''

        # Ta calculation according to
        # Rychlik, Spencer, and Rhoads, 1990, Optimization of the anneal
        # ing temperature for DNA amplification in vitro
        # http://www.ncbi.nlm.nih.gov/pubmed/2003928
        GC_prod=GC(str(self.seq))
        tml = min(self.tmf,self.tmr)
        #print GC(str(self.product.seq)), self.saltc/1000.0, len(self.product)
        tmp = 81.5 + 0.41*GC(str(self.seq)) + 16.6*log10(self.saltc/1000.0) - 675/len(self)
        ta = 0.3*tml+0.7*tmp-14.9
        # Fermentas recombinant taq
        taq_extension_rate = 30  # seconds/kB PCR product length
        extension_time_taq = taq_extension_rate * len(self) / 1000 # seconds
        f  = textwrap.dedent(u'''
                                 Taq (rate {rate} nt/s) 35 cycles             |{size}bp
                                 95.0°C    |95.0°C                 |      |SantaLucia 1998
                                 |_________|_____          72.0°C  |72.0°C|SaltC {saltc:2}mM
                                 | 03min00s|30s  \         ________|______|
                                 |         |      \ {ta}°C/{0:2}min{1:2}s| 5min |
                                 |         |       \_____/         |      |
                                 |         |         30s           |      |4-12°C'''.format(rate      = taq_extension_rate,
                                           ta      = math.ceil(ta),
                                           saltc   = self.saltc,
                                            *divmod(extension_time_taq,60),
                                            size = len(self.seq)))

        return pretty_unicode(f)

    def taq_program(self):
        return self.program()

    def dbd_program(self):
        '''Returns a string containing a text representation of a proposed
       PCR program using a polymerase with a DNA binding domain such as Pfu-Sso7d.

       ::

        Pfu-Sso7d (rate 15s/kb)
        Three-step|          30 cycles   |      |Breslauer1986,SantaLucia1998
        98.0°C    |98.0°C                |      |SaltC 50mM
        __________|_____          72.0°C |72.0°C|Primer1C   1µM
        00min30s  |10s  \ 61.0°C ________|______|Primer2C   1µM
                  |      \______/ 0min 0s|10min |
                  |        10s           |      |4-8°C

       '''
        PfuSso7d_extension_rate = 15 #seconds/kB PCR product
        extension_time_PfuSso7d = PfuSso7d_extension_rate * len(self) / 1000  # seconds

        # Ta calculation for enzymes with dsDNA binding domains like Pfu-Sso7d
        # https://www.finnzymes.fi/tm_determination.html

        length_of_f = len(self.forward_primer.footprint)
        length_of_r = len(self.reverse_primer.footprint)

        if (length_of_f>20 and length_of_r>20 and self.tmf_dbd>=69.0 and self.tmr_dbd>=69.0) or (self.tmf_dbd>=72.0 and self.tmr_dbd>=72.0):
            f=textwrap.dedent(  '''
                                    Pfu-Sso7d (rate {rate}s/kb)
                                    Two-step|    30 cycles |      |{size}bp
                                    98.0°C  |98.0C         |      |Breslauer1986,SantaLucia1998
                                    _____ __|_____         |      |SaltC {saltc:2}mM
                                    00min30s|10s  \  72.0°C|72.0°C|Primer1C {forward_primer_concentration:3}µM
                                            |      \_______|______|Primer2C {reverse_primer_concentration:3}µM
                                            |      {0:2}min{1:2}s|10min |4-8°C
                                 '''.format(rate = PfuSso7d_extension_rate,
                                            forward_primer_concentration = self.forward_primer_concentration,
                                            reverse_primer_concentration = self.rc,
                                            saltc = self.saltc,
                                            *divmod(extension_time_PfuSso7d,60),
                                            size = len(self.seq)))
        else:

            if (length_of_f>20 and length_of_r>20):
                ta = min(self.tmf_dbd,self.tmr_dbd)+3
            else:
                ta = min(self.tmf_dbd,self.tmr_dbd)

            f=textwrap.dedent(  '''
                                    Pfu-Sso7d (rate {rate}s/kb)
                                    Three-step|          30 cycles   |      |Breslauer1986,SantaLucia1998
                                    98.0°C    |98.0°C                |      |SaltC {saltc:2}mM
                                    __________|_____          72.0°C |72.0°C|Primer1C {forward_primer_concentration:3}µM
                                    00min30s  |10s  \ {ta}°C ________|______|Primer2C {reverse_primer_concentration:3}µM
                                              |      \______/{0:2}min{1:2}s|10min |
                                              |        10s           |      |4-8°C
                                 '''.format(rate = PfuSso7d_extension_rate,
                                            ta   = math.ceil(ta),
                                            forward_primer_concentration   = self.forward_primer_concentration/1000,
                                            reverse_primer_concentration   = self.reverse_primer_concentration/1000,
                                            saltc= self.saltc,
                                            *divmod(extension_time_PfuSso7d,60)))
        return pretty_str(f)



class Anneal(object):
    u'''

    Parameters
    ----------
    primers : iterable containing SeqRecord objects
        Primer sequences 5'-3'.

    template : Dseqrecord object
        The template sequence 5'-3'.

    limit : int, optional
        limit length of the annealing part of the primers.

    Attributes
    ----------
    products: list
        A list of Amplicon objects, one for each primer pair that may form a PCR product.


    Examples
    --------
    >>> import pydna
    >>> template = pydna.Dseqrecord("tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
    >>> from Bio.SeqRecord import SeqRecord
    >>> p1 = pydna.read(">p1\\ntacactcaccgtctatcattatc", ds = False)
    >>> p2 = pydna.read(">p2\\ngtgctatcagatgatacagtcg", ds = False)
    >>> ann = pydna.Anneal((p1, p2), template)
    >>> print ann.report()
    Template na 51 nt linear:
    Primer p1 anneals forward at position 23
    <BLANKLINE>
    Primer p2 anneals reverse at position 29
    >>> ann.products
    [Amplicon(51)]
    >>> amplicon_list = ann.products
    >>> amplicon = amplicon_list.pop()
    >>> amplicon
    Amplicon(51)
    >>> print amplicon.figure()
    5tacactcaccgtctatcattatc...cgactgtatcatctgatagcac3
                               |||||||||||||||||||||| tm 50.6 (dbd) 60.5
                              3gctgacatagtagactatcgtg5
    5tacactcaccgtctatcattatc3
     ||||||||||||||||||||||| tm 49.4 (dbd) 58.8
    3atgtgagtggcagatagtaatag...gctgacatagtagactatcgtg5
    >>> amplicon.annotations['date'] = '02-FEB-2013'   # Set the date for this example to pass the doctest
    >>> print amplicon
    Dseqrecord
    circular: False
    size: 51
    ID: 51bp U96-TO06Y6pFs74SQx8M1IVTBiY
    Name: 51bp_PCR_prod
    Description: Product_p1_p2
    Number of features: 4
    /date=02-FEB-2013
    Dseq(-51)
    taca..gcac
    atgt..cgtg
    >>> print amplicon.program()
    <BLANKLINE>
    Taq (rate 30 nt/s) 35 cycles             |51bp
    95.0°C    |95.0°C                 |      |SantaLucia 1998
    |_________|_____          72.0°C  |72.0°C|SaltC 50mM
    | 03min00s|30s  \         ________|______|
    |         |      \ 44.0°C/ 0min 1s| 5min |
    |         |       \_____/         |      |
    |         |         30s           |      |4-12°C
    >>>

    '''
    def __init__( self,
                  primers,
                  template,
                  limit=13):

        refresh = False
        cached  = None

        primers = [p for p in primers if p.seq]

        key = str(template.seguid()) + "|".join(sorted([seguid(p.seq) for p in primers]))+str(limit)

        if os.environ["pydna_cache"] in ("compare", "cached"):
            cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "amplify"), protocol=cPickle.HIGHEST_PROTOCOL, writeback=False)
            try:
                cached = cache[key]
            except:
                if os.environ["pydna_cache"] == "compare":
                    raise Exception("no result for this key!")
                else:
                    refresh = True
            cache.close()

        if refresh or os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):

            self.key = key
            self.template = copy.deepcopy(template)

            self.limit = limit
            self._products = None

            self.fwd_primers = []
            self.rev_primers = []

            twl = len(self.template.seq.watson)
            tcl = len(self.template.seq.crick)

            if self.template.linear:
                tw = self.template.seq.watson
                tc = self.template.seq.crick
            else:
                tw = self.template.seq.watson+self.template.seq.watson
                tc = self.template.seq.crick +self.template.seq.crick

            for p in primers:
                self.fwd_primers.extend((Primer(p,
                                                tcl-pos - min(self.template.seq.ovhg, 0),
                                                Seq(fp), Seq(tl))
                                        for pos, fp, tl in _annealing_positions(
                                                            str(p.seq),
                                                            tc,
                                                            limit) if pos<tcl))
                self.rev_primers.extend((Primer(p,
                                                pos + max(0, self.template.seq.ovhg),
                                                Seq(fp), Seq(tl))
                                         for pos, fp, tl in _annealing_positions(
                                                                         str(p.seq),
                                                                         tw,
                                                                         limit) if pos<twl))
            self.fwd_primers.sort(key = operator.attrgetter('position'))
            self.rev_primers.sort(key = operator.attrgetter('position'), reverse=True)

            for fp in self.fwd_primers:
                if fp.position-len(fp.footprint)>=0:
                    start = fp.position-len(fp.footprint)
                    end   = fp.position
                    self.template.features.append(SeqFeature(FeatureLocation(start, end),
                                                        type ="primer_bind",
                                                        strand = 1,
                                                        qualifiers = {"note":[fp.name],
                                                                      "ApEinfo_fwdcolor":["green"],
                                                                      "ApEinfo_revcolor":["red"]}))
                else:
                    start = len(self.template)-len(fp.footprint)+fp.position
                    end = start+len(fp.footprint)-len(self.template)
                    sf=SeqFeature(CompoundLocation([FeatureLocation(start,len(self.template)),
                                                    FeatureLocation(0, end)]),
                                                    type="primer_bind",
                                                    location_operator="join",
                                                    qualifiers = {"note":[fp.name]})

                    self.template.features.append(sf)

            for rp in self.rev_primers:
                if rp.position+len(rp.footprint)<=len(self.template):
                    start = rp.position
                    end   = rp.position + len(rp.footprint)
                    self.template.features.append(SeqFeature(FeatureLocation(start,end),
                                                        type ="primer_bind",
                                                        strand = -1,
                                                        qualifiers = {"note":[rp.name],
                                                                      "ApEinfo_fwdcolor":["green"],
                                                                      "ApEinfo_revcolor":["red"]}))
                else:
                    start = rp.position
                    end = rp.position+len(rp.footprint)-len(self.template)
                    self.template.features.append(SeqFeature(CompoundLocation([FeatureLocation(start,len(self.template)),
                                                                               FeatureLocation(0,end)]),
                                                        type ="primer_bind",
                                                        location_operator= "join",
                                                        strand = -1,
                                                        qualifiers = {"note":[rp.name]}))
            self.forward_primers = self.fwd_primers
            self.reverse_primers = self.rev_primers

        if os.environ["pydna_cache"] == "compare":
            self._compare(cached)

        if refresh or os.environ["pydna_cache"] == "refresh":
            self._save()

        elif cached and os.environ["pydna_cache"] not in ("nocache","refresh"):
            for key, value in cached.__dict__.items():
                setattr(self, key, value )
            cache.close()

    def _compare(self, cached):
        if str(self) != str(cached):
            module_logger.warning('amplify error')

    def _save(self):
        cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "amplify"), protocol=cPickle.HIGHEST_PROTOCOL, writeback=False)
        cache[self.key] = self
        cache.close()

    @property
    def products(self):

        if self._products:
            return self._products

        self._products = []

        for fp in self.fwd_primers:
            for rp in self.rev_primers:

                if self.template.circular and fp.position>rp.position:
                    tmpl = self.template.shifted(fp.position-len(fp.footprint))
                    tmpl = tmpl._multiply_circular(2)
                    tmpl = tmpl[:len(self.template) - (fp.position - rp.position) + len(rp.footprint) + len(fp.footprint)]
                    #print len(self.template) - (fp.position - rp.position) + len(rp.footprint) + len(fp.footprint)
                    #print len(tmpl)
                elif self.template.circular:
                    tmpl = self.template._multiply_circular(3)
                    tmpl = tmpl[fp.position-len(fp.footprint)+len(self.template):rp.position+len(rp.footprint)+len(self.template)]
                else:
                    tmpl = self.template[fp.position-len(fp.footprint):rp.position+len(rp.footprint)]

                prd = ( Dseqrecord(fp.tail) + tmpl + Dseqrecord(rp.tail).reverse_complement())

                prd.add_feature( 0, len(fp), label=fp.id)
                prd.add_feature( len(prd)-len(rp),len(prd),label=rp.id, strand=-1)

                #prd.seq = fp.seq+tmpl.seq[len(fp.footprint):len(tmpl)-len(rp.footprint)]+rp.seq.reverse_complement()

                #features = tmpl[fp.position-len(fp.footprint):rp.position+len(rp.footprint)].features
                #print   fp.position-len(fp.footprint), rp.position+len(rp.footprint), features
                #print "<<<<<<<<<<<<<<",features
                #prd.features = [f._shift(len(fp.tail)) for f in features]
                # description = Genbank LOCUS max 16 chars

                prd.name = "{0}bp_PCR_prod".format(len(prd))[:16]
                prd.id = "{0}bp {1}".format( str(len(prd))[:14], prd.seguid() )
                prd.description="Product_{0}_{1}".format( fp.description,
                                                          rp.description)

                self._products.append( Amplicon(prd,
                                                template=tmpl,
                                                forward_primer=fp,
                                                reverse_primer=rp,
                                                saltc=50,
                                                forward_primer_concentration=1000,
                                                reverse_primer_concentration=1000))
                assert " " not in str(prd.seq.watson)
                assert " " not in str(prd.seq.crick)

        return self._products



    def report(self):
        '''This method is an alias of str(Annealobj).
        Returns a short string representation.
        '''
        return self.__str__()

    def __repr__(self):
        ''' returns a short string representation '''
        return "Reaction(products = {})".format(len(self.fwd_primers*len(self.rev_primers)))

    def __str__(self):
        '''returns a short report describing if or where primer
       anneal on the template.
       '''
        mystring = "Template {name} {size} nt {top}:\n".format(name=self.template.name,
                                                               size=len(self.template),
                                                               top={True:"circular",
                                                                    False:"linear"}[self.template.circular]
                                                                    )
        if self.fwd_primers:
            for p in self.fwd_primers:
                mystring += "Primer {name} anneals forward at position {pos}\n".format(name=p.name, pos=p.position)
        else:
            mystring += "No forward primers anneal...\n"
        mystring +="\n"
        if self.rev_primers:
            for p in self.rev_primers:
                mystring += "Primer {name} anneals reverse at position {pos}\n".format(name=p.name, pos=p.position)
        else:
             mystring += "No reverse primers anneal...\n"
        return mystring.strip()






def pcr(*args,  **kwargs):
    '''pcr is a convenience function for Anneal to simplify its usage,
    especially from the command line. If more than one PCR product is
    formed, an exception is raised.

    args is any iterable of sequences or an iterable of iterables of sequences.
    args will be greedily flattened.

    Parameters
    ----------

    args : iterable containing sequence objects
        Several arguments are also accepted.

    limit : int = 13, optional
        limit length of the annealing part of the primers.

    Notes
    -----

    sequences in args could be of type:

    * string
    * Seq
    * SeqRecord
    * Dseqrecord

    The last sequence will be interpreted as the template
    all preceeding sequences as primers.

    This is a powerful function, use with care!

    Returns
    -------

    product : Dseqrecord
        a Dseqrecord object representing the PCR product.
        The direction of the PCR product will be the same as
        for the template sequence.

    Examples
    --------

    >>> import pydna
    >>> template = pydna.Dseqrecord("tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
    >>> from Bio.SeqRecord import SeqRecord
    >>> p1 = pydna.read(">p1\\ntacactcaccgtctatcattatc", ds = False)
    >>> p2 = pydna.read(">p2\\ncgactgtatcatctgatagcac", ds = False).reverse_complement()
    >>> pydna.pcr(p1, p2, template)
    Amplicon(51)
    >>> pydna.pcr([p1, p2], template)
    Amplicon(51)
    >>> pydna.pcr((p1,p2,), template)
    Amplicon(51)
    >>>

    '''

    from Bio.SeqRecord import SeqRecord
    # flatten args
    output = []
    stack = []
    stack.extend(reversed(args))
    while stack:
        top = stack.pop()
        if hasattr(top, "__iter__") and not isinstance(top, SeqRecord):
            stack.extend(reversed(top))
        else:
            output.append(top)
    new=[]

    for s in output:
        if isinstance(s, Seq):
            s = SeqRecord(s)
        elif isinstance(s, SeqRecord):
            pass
        elif hasattr(s, "watson"):
            s=s.watson
        elif isinstance(s, basestring):
            s = SeqRecord(Seq(s))
        else:
            raise TypeError("the record property needs to be a string, a Seq object or a SeqRecord object")
        new.append(s)

    anneal_primers = Anneal(  new[:-1],
                              new[-1],
                              **kwargs)

    if anneal_primers:
        if len(anneal_primers.products) == 1:
            return anneal_primers.products.pop()
        elif len(anneal_primers.products) == 0:
            raise Exception("No PCR products! {}".format(anneal_primers.report()))
        else:
            raise Exception("PCR not specific! {}".format(anneal_primers.report()))
    else:
        raise Exception(anneal_primers.report())
    return


def nopcr(*args,  **kwargs):
    '''pcr is a convenience function for Anneal to simplify its usage,
    especially from the command line. If more one or more PCR products are
    formed, an exception is raised.

    args is any iterable of sequences or an iterable of iterables of sequences.
    args will be greedily flattened.

    Parameters
    ----------

    args : iterable containing sequence objects
        Several arguments are also accepted.

    limit : int = 13, optional
        limit length of the annealing part of the primers.

    Notes
    -----

    sequences in args could be of type:

    * string
    * Seq
    * SeqRecord
    * Dseqrecord

    The last sequence will be interpreted as the template
    all preceeding sequences as primers.

    This is a powerful function, use with care!

    Returns
    -------

    product : Dseqrecord
        a Dseqrecord object representing the PCR product.
        The direction of the PCR product will be the same as
        for the template sequence.

    Examples
    --------

    >>> import pydna
    >>> template = pydna.Dseqrecord("tacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Seq import Seq
    >>> p1 = SeqRecord(Seq("tacactcaccgtctatcattatc"))
    >>> p2 = SeqRecord(Seq("gtgctatcagatgatacagtG")) # This primer does not anneal
    >>> pydna.nopcr(p1, p2, template)
    True
    '''

    from Bio.SeqRecord import SeqRecord
    # flatten args
    output = []
    stack = []
    stack.extend(reversed(args))
    while stack:
        top = stack.pop()
        if hasattr(top, "__iter__") and not isinstance(top, SeqRecord):
            stack.extend(reversed(top))
        else:
            output.append(top)
    new=[]

    for s in output:
        if isinstance(s, Seq):
            s = SeqRecord(s)
        elif isinstance(s, SeqRecord):
            pass
        elif hasattr(s, "watson"):
            s=s.watson
        elif isinstance(s, basestring):
            s = SeqRecord(Seq(s))
        else:
            raise TypeError("the record property needs to be a string, a Seq object or a SeqRecord object")
        new.append(s)

    anneal_primers = Anneal(  new[:-1],
                              new[-1],
                              **kwargs)

    if anneal_primers:
        if len(anneal_primers.products) != 0:
            raise Exception("PCR products formed! {}".format(anneal_primers.report()))
    else:
        raise Exception(anneal_primers.report())
    return True



def basictm(primer, *args, **kwargs):
    '''Returns the melting temperature (Tm) of the primer using
    the basic formula.

    | Tm = (wA+xT)*2 + (yG+zC)*4 assumed 50mM monovalent cations
    |
    | w = number of A in primer
    | x = number of T in primer
    | y = number of G in primer
    | z = number of C in primer

    Parameters
    ----------
    primer : string
        Primer sequence 5'-3'

    Returns
    -------
    tm : int

    Examples
    --------
    >>> from pydna.amplify import basictm
    >>> basictm("ggatcc")
    20
    >>>

    '''
    primer = str(primer).lower()
    return (primer.count("a") + primer.count("t"))*2 + (primer.count("g") + primer.count("c"))*4

# http://www.promega.com/techserv/tools/biomath/calc11.htm#melt_results

def tmstaluc98(primer, dnac=50, saltc=50, **kwargs):
    '''Returns the melting temperature (Tm) of the primer using
    the nearest neighbour algorithm. Formula and thermodynamic data
    is taken from SantaLucia 1998 [1]_. This implementation gives the same
    answer as the one provided by Biopython (See Examples).

    Thermodynamic data used:

    =====  ====  ====
    pair   dH    dS
    =====  ====  ====
    AA/TT  7.9   22.2
    AT/TA  7.2   20.4
    TA/AT  7.2   21.3
    CA/GT  8.5   22.7
    GT/CA  8.4   22.4
    CT/GA  7.8   21.0
    GA/CT  8.2   22.2
    CG/GC  10.6  27.2
    GC/CG  9.8   24.4
    GG/CC  8.0   19.9
    =====  ====  ====

    Parameters
    ----------
    primer : string
        Primer sequence 5'-3' in UPPERCASE

    Returns
    -------
    tm : float
        tm of the primer


    Examples
    --------

    >>> from pydna.amplify import tmstaluc98
    >>> from Bio.SeqUtils.MeltingTemp import Tm_staluc
    >>> tmstaluc98("ACGTCATCGACACTATCATCGAC")
    54.55597724052518
    >>>



    '''

    nntermsl={  "AA": (7.9  , 22.2),
                "TT": (7.9  , 22.2),
                "AT": (7.2  , 20.4),
                "TA": (7.2  , 21.3),
                "CA": (8.5  , 22.7),
                "TG": (8.5  , 22.7),
                "GT": (8.4  , 22.4),
                "AC": (8.4  , 22.4),
                "CT": (7.8  , 21.0),
                "AG": (7.8  , 21.0),
                "GA": (8.2  , 22.2),
                "TC": (8.2  , 22.2),
                "CG": (10.6 , 27.2),
                "GC": (9.8  , 24.4),
                "GG": (8    , 19.9),
                "CC": (8    , 19.9),
                "A" : (0    , 0   ),
                "C" : (0    , 0   ),
                "G" : (0    , 0   ),
                "T" : (0    , 0   )  }

    helixinit = {   "G": (-0.1 ,2.8),
                    "C": (-0.1 ,2.8),
                    "A": (-2.3, -4.1),
                    "T": (-2.3, -4.1) }

    dH, dS = helixinit[primer[0]]
    H ,  S = helixinit[primer[-1]]
    dH = dH+H
    dS = dS+S

    for p in range(len(primer)):
        dn = primer[p:p+2]
        H,S = nntermsl[dn]
        dH+=H
        dS+=S
    R = 1.987 # universal gas constant in Cal/degrees C*Mol
    k = (dnac/4.0)*1e-9
    dS = dS-0.368*(len(primer)-1)*math.log(float(saltc)/1e3)
    tm = ((1000* (-dH))/(-dS+(R * (math.log(k)))))-273.15
    return tm

def tmbreslauer86(primer, dnac=500.0, saltc=50,thermodynamics=False):
    '''Returns the melting temperature (Tm) of the primer using
    the nearest neighbour algorithm. Formula and thermodynamic data
    is taken from Breslauer 1986. These data are no longer widely used.


    Breslauer 1986, table 2 [#]_

    =====  ===== ====   ===
    pair   dH    dS     dG
    =====  ===== ====   ===
    AA/TT  9.1   24.0   1.9
    AT/TA  8.6   23.9   1.5
    TA/AT  6.0   16.9   0.9
    CA/GT  5.8   12.9   1.9
    GT/CA  6.5   17.3   1.3
    CT/GA  7.8   20.8   1.6
    GA/CT  5.6   13.5   1.6
    CG/GC  11.9  27.8   3.6
    GC/CG  11.1  26.7   3.1
    GG/CC  11.0  26.6   3.1
    =====  ===== ====   ===

    Parameters
    ----------
    primer : string
        Primer sequence 5'-3' in UPPERCASE

    Returns
    -------
    tm : float


    References
    ----------
    .. [#] K.J. Breslauer et al., “Predicting DNA Duplex Stability from the Base Sequence,” Proceedings of the National Academy of Sciences 83, no. 11 (1986): 3746.


    Examples
    ---------

    >>> from pydna.amplify import tmbreslauer86
    >>> tmbreslauer86("ACGTCATCGACACTATCATCGAC")
    64.28863985851899


    '''

    nntermbr={  "AA": (9.1   ,24.0   ,1.9),
                "TT": (9.1   ,24.0   ,1.9),
                "AT": (8.6   ,23.9   ,1.5),
                "TA": (6.0   ,16.9   ,0.9),
                "CA": (5.8   ,12.9   ,1.9),
                "TG": (5.8   ,12.9   ,1.9),
                "GT": (6.5   ,17.3   ,1.3),
                "AC": (6.5   ,17.3   ,1.3),
                "CT": (7.8   ,20.8   ,1.6),
                "AG": (7.8   ,20.8   ,1.6),
                "GA": (5.6   ,13.5   ,1.6),
                "TC": (5.6   ,13.5   ,1.6),
                "CG": (11.9  ,27.8   ,3.6),
                "GC": (11.1  ,26.7   ,3.1),
                "GG": (11.0  ,26.6   ,3.1),
                "CC": (11.0  ,26.6   ,3.1),
                "A" : (0     , 0     ,0),
                "C" : (0     , 0     ,0),
                "G" : (0     , 0     ,0),
                "T" : (0     , 0     ,0),     }
    dH=3.4
    dS=12.4
    dG=0
    for p in range(len(primer)):
        dn = primer[p:p+2]
        H,S,G = nntermbr[dn]
        dG+=G
        dH+=H
        dS+=S

    R = 1.9872          # universal gas constant in Cal/degrees C*Mol
    k = dnac*1E-9/2.0
    dH = dH - 5
    dS = dS-0.368*(len(primer)-1)*math.log(float(saltc)/1E3) # SantaLucia salt correction formula
    tm = 1000 * -dH /(-dS + R * math.log(k) )  - 273.15 # degrees Celsius

    if thermodynamics:
        return tm,dH,dS
    else:
        return tm


def tmbresluc(primer, primerc=500.0, saltc=50, thermodynamics=False):
    '''Returns the tm for a primer using a formula adapted to polymerases
    with a DNA binding domain.

    Parameters
    ----------

    primer : string
        primer sequence 5'-3'

    primerc : float
       concentration (nM)

    saltc : float, optional
       Monovalent cation concentration (mM)

    thermodynamics : bool, optional
        prints details of the thermodynamic data to stdout. For
        debugging only.

    Returns
    -------
    tm : float
        the tm of the primer

    tm,dH,dS : tuple
        tm and dH and dS used for the calculation

    '''

    import thermodynamic_data

    saltc = float(saltc)/1000
    pri  = primerc/10E7
    dS = -12.4
    dH = -3400

    STR = primer.lower();

    for i in range(len(STR)-1):
        n1=ord(STR[i])
        n2=ord(STR[i+1])
        dH += thermodynamic_data.dHBr[n1 - 97][n2 - 97]
        dS += thermodynamic_data.dSBr[n1 - 97][n2 - 97]

    tm = (dH / (1.9872 * math.log(pri / 1600) + dS) + (16.6 * math.log(saltc)) / math.log(10)) - 273.15

    if thermodynamics:
        return tm,dH,dS
    else:
        return tm

if __name__=="__main__":
    import doctest
    doctest.testmod()
