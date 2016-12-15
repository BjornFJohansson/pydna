#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# doctest: +NORMALIZE_WHITESPACE
# doctest: +SKIP
'''This module provides functions for PCR. Primers with 5' tails as well as inverse PCR on
circular templates are handled correctly.

'''

import pickle
import shelve

import itertools

import re
import copy
import operator
import os

import logging
module_logger = logging.getLogger("pydna."+__name__)

from Bio.Seq                        import Seq
from Bio.SeqUtils.CheckSum          import seguid
from Bio.SeqFeature                 import SeqFeature
from Bio.SeqFeature                 import CompoundLocation
from Bio.SeqFeature                 import FeatureLocation
from .utils                         import rc
from .dseqrecord                    import Dseqrecord
from .primer                        import Primer
from .amplicon                      import Amplicon


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
            footprint = rc(template[match_start:match_start+limit]+"".join([b for a,b in itertools.takewhile(lambda x: x[0].lower()==x[1].lower(), list(zip(tail, tm)))]))
            results.append((match_start, footprint, primer[: len(primer) - len(footprint) ]))
        return results
    return []

class Anneal(object):
    '''

    Parameters
    ----------
    primers : iterable containing SeqRecord objects
        Primer sequences 5'-3'.

    template : Dseqrecord object
        The template sequence 5'-3'.

    limit : int, optional
        limit length of the annealing part of the primers.

    fprimerc : float, optional
        Concentration of forward primer in nM, set to 1000.0 nM by default

    rprimerc : float, optional
        Concentration of reverse primer in nM, set to 1000.0 nM by default

    saltc  : float, optional
        Salt concentration (monovalet cations) :mod:`tmbresluc` set to 50.0 mM by default

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
    >>> print(ann.report())
    Template name? 51 nt linear:
    Primer p1 anneals forward at position 23
    <BLANKLINE>
    Primer p2 anneals reverse at position 29
    >>> ann.products
    [Amplicon(51)]
    >>> amplicon_list = ann.products
    >>> amplicon = amplicon_list.pop()
    >>> amplicon
    Amplicon(51)
    >>> print(amplicon.figure())
    5tacactcaccgtctatcattatc...cgactgtatcatctgatagcac3
                               |||||||||||||||||||||| tm 55.9 (dbd) 60.5
                              3gctgacatagtagactatcgtg5
    5tacactcaccgtctatcattatc3
     ||||||||||||||||||||||| tm 54.6 (dbd) 58.8
    3atgtgagtggcagatagtaatag...gctgacatagtagactatcgtg5
    >>> amplicon.annotations['date'] = '02-FEB-2013'   # Set the date for this example to pass the doctest
    >>> print(amplicon)
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
    >>> print(amplicon.program())
    <BLANKLINE>
    Taq (rate 30 nt/s) 35 cycles             |51bp
    95.0°C    |95.0°C                 |      |Tm formula: Biopython Tm_NN
    |_________|_____          72.0°C  |72.0°C|SaltC 50mM
    | 03min00s|30s  \         ________|______|Primer1C 1.0µM
    |         |      \ 45.4°C/ 0min 2s| 5min |Primer2C 1.0µM
    |         |       \_____/         |      |GC 39%
    |         |         30s           |      |4-12°C

    >>>

    '''
    def __init__( self,
                  primers,
                  template,
                  limit=13,
                  primerc=1000.0, # nM
                  saltc=50):

        refresh = False
        cached  = None

        self.primerc=primerc
        self.saltc = saltc

        primers = [p for p in primers if p.seq]

        key = str(template.seguid()) + "|".join(sorted([seguid(p.seq) for p in primers]))+str(limit)

        if os.environ["pydna_cache"] in ("compare", "cached"):
            cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "amplify"), protocol=pickle.HIGHEST_PROTOCOL, writeback=False)
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
                                                position = tcl-pos - min(self.template.seq.ovhg, 0),
                                                footprint = Seq(fp),
                                                tail = Seq(tl))
                                        for pos, fp, tl in _annealing_positions(
                                                            str(p.seq),
                                                            tc,
                                                            limit) if pos<tcl))
                self.rev_primers.extend((Primer(p,
                                                position = pos + max(0, self.template.seq.ovhg),
                                                footprint = Seq(fp),
                                                tail = Seq(tl))
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
            for key, value in list(cached.__dict__.items()):
                setattr(self, key, value )
            cache.close()

    def _compare(self, cached):
        if str(self) != str(cached):
            module_logger.warning('amplify error')

    def _save(self):
        cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "amplify"), protocol=pickle.HIGHEST_PROTOCOL, writeback=False)
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

                elif self.template.circular:
                    tmpl = self.template._multiply_circular(3)
                    tmpl = tmpl[fp.position-len(fp.footprint)+len(self.template):rp.position+len(rp.footprint)+len(self.template)]
                else:
                    tmpl = self.template[fp.position-len(fp.footprint):rp.position+len(rp.footprint)]

                prd = ( Dseqrecord(fp.tail) + tmpl + Dseqrecord(rp.tail).reverse_complement())

                prd.add_feature( 0, len(fp), label=fp.id)
                prd.add_feature( len(prd)-len(rp),len(prd),label=rp.id, strand=-1)

                prd.name = "{0}bp_PCR_prod".format(len(prd))[:16]
                prd.id = "{0}bp {1}".format( str(len(prd))[:14], prd.seguid() )
                prd.description="Product_{0}_{1}".format( fp.description,
                                                          rp.description)

                self._products.append( Amplicon(prd,
                                                template=tmpl,
                                                forward_primer=fp,
                                                reverse_primer=rp,
                                                saltc=self.saltc,
                                                fprimerc=self.primerc,
                                                rprimerc=self.primerc))
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
        elif isinstance(s, str):
            s = SeqRecord(Seq(s))
        else:
            raise TypeError("the record property needs to be a string, Seq, SeqRecord or Dseqrecord object")
        new.append(s)

    anneal_primers = Anneal(  new[:-1],
                              new[-1],
                              **kwargs)

    if anneal_primers:
        if len(anneal_primers.products) == 1:
            return anneal_primers.products[0]
        elif len(anneal_primers.products) == 0:
            raise Exception("No PCR products! {}".format(anneal_primers.report()))
        else:
            raise Exception("PCR not specific! {}".format(anneal_primers.report()))
    else:
        raise Exception(anneal_primers.report())
    return


def nopcr(*args,  **kwargs):
    '''no pcr is a convenience function for Anneal to simplify its usage,
    especially from the command line. If one or more PCR products are
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
        elif isinstance(s, str):
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


if __name__=="__main__":
    import os
    cache = os.getenv("pydna_cache")
    os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod()
    os.environ["pydna_cache"]=cache
