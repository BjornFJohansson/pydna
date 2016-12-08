#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''Provides a class for downloading sequences from genbank.
'''
import pickle
import shelve
import re
import os
import requests
import textwrap
from Bio import Entrez
from pydna.dsdna import read
from pydna.dsdna import Dseqrecord

class GenbankRecord(Dseqrecord):

    def __init__(self, record, *args, item = "", start=None, stop=None, strand=1,**kwargs):
        super().__init__(record, *args, **kwargs)
        self.item = item
        self.start = start
        self.stop = stop
        self.strand = strand

    def __repr__(self):
        '''returns a short string representation of the object'''
        return "Genbank({})({}{})".format(self.id, {True:"-", False:"o"}[self.linear],len(self))
        
    def _repr_pretty_(self, p, cycle):
        '''returns a short string representation of the object'''
        p.text("Genbank({})({}{})".format(self.id, {True:"-", False:"o"}[self.linear],len(self)))
            
    def _repr_html_(self):
        if not self.item:
            return self.__repr__()
        linktext = self.item
        if self.start != None and self.stop != None:
            linktext += " {}-{}".format(self.start, self.stop)
        return "<a href='https://www.ncbi.nlm.nih.gov/nuccore/{item}?from={start}&to={stop}&strand={strand}' target='_blank'>{linktext}</a>".format(item=self.item, start=self.start or "", stop=self.stop or "", strand=self.strand, linktext=linktext)

class Genbank(object):
    '''Class to facilitate download from genbank.

    Parameters
    ----------
    users_email : string
        Has to be a valid email address. You should always tell
        Genbanks who you are, so that they can contact you.
    tool : string, optional
        Default is "pydna". This is to tell Genbank which tool you are
        using.

    Examples
    --------

    >>> import pydna
    >>> gb=pydna.Genbank("bjornjobb@gmail.com")
    >>> rec = gb.nucleotide("L09137") # <- pUC19 from genbank
    >>> print(len(rec))
    2686
    '''

    def __init__(self, users_email, tool="pydna"):

        if not re.match("[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}",users_email,re.IGNORECASE):
            raise ValueError
        if email == "someone@example.com":
            raise ValueError("you have to set your email address in order to download from Genbank")
        self.email=users_email

    def __repr__(self):
        return "GenbankConnection({})".format(self.email)

    def nucleotide(self, item, start=None, stop=None, strand="watson" ):
        '''Download a genbank nuclotide record.

        Item is a string containing one genbank acession number.
        for a nucleotide file. Genbank nucleotide accession numbers have this format:

        | A12345   = 1 letter  + 5 numerals
        | AB123456 = 2 letters + 6 numerals
        
        The accession number is sometimes followed by a pint and version number
        
        | BK006936.2
        
        Item can also contain optional interval information in the following formats:
        
        | BK006936.2 REGION: complement(613900..615202)
        | NM_005546 REGION: 1..100
        | NM_005546 REGION: complement(1..100)
        | 21614549:1-100
        | 21614549:c100-1
        
        Start and stop are the sequence intervals to be downloaded. This 
        is useful for large genbank records.
        If strand is "c", "C", "crick", "Crick", "antisense","Antisense",
        "2" or 2, the antisense (Crick) strand is returned, otherwise
        the sense (Watson) strand is returned.

        Result is returned as a GenbankRecord object.

        References
        ----------

        .. [#]   http://www.dsimb.inserm.fr/~fuchs/M2BI/AnalSeq/Annexes/Sequences/Accession_Numbers.htm
        .. [#]   http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
        '''
        cached  = False
        refresh = False

        if os.environ["pydna_cache"] in ("compare", "cached"):
            cache = shelve.open(os.path.join(os.environ["pydna_data_dir"],"genbank"), protocol=pickle.HIGHEST_PROTOCOL, writeback=False)
            key = item+str(start)+str(stop)+str(strand)
            try:
                cached, item, start, stop = cache[key]
            except:
                if os.environ["pydna_cache"] == "compare":
                    raise Exception("no result for this key!")
                else:
                    refresh = True

        if refresh or os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):

            matches =((1, re.search("(REGION:\s(?P<start>\d+)\.\.(?P<stop>\d+))", item)),
                      (2, re.search("(REGION: complement\((?P<start>\d+)\.\.(?P<stop>\d+)\))",item)),
                      (1, re.search(":(?P<start>\d+)-(?P<stop>\d+)",item)),
                      (2, re.search(":c(?P<start>\d+)-(?P<stop>\d+)",item)),
                      (0, None))

            for strand, match in matches:
                if match:
                    start = match.group("start")
                    stop  = match.group("stop")
                    item = item[:match.start()]
                    break

            if str(strand).lower() in ("c","crick", "antisense", "2"):
                strand = 2
            else:
                strand = 1

            Entrez.email = self.email

            text = Entrez.efetch( db        ="nucleotide",
                                  id        = item,
                                  rettype   = "gbwithparts",
                                  seq_start = start,
                                  seq_stop  = stop,
                                  strand    = strand,
                                  retmode   = "text" ).read()
            dsr = read(text)
            result = GenbankRecord(dsr, item = item, start=start, stop=stop, strand=strand)
    
        if os.environ["pydna_cache"] == "compare":
            if result!=cached:
                module_logger.warning('download error')

        if refresh or os.environ["pydna_cache"] == "refresh":
            cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "genbank"), protocol=pickle.HIGHEST_PROTOCOL, writeback=False)
            cache[key] = result, item, start, stop
        elif cached and os.environ["pydna_cache"] not in ("nocache", "refresh"):
            result = cached
            cache.close()

        return result

def download_text(url):
    cached  = False
    refresh = False
    cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "web"), protocol=pickle.HIGHEST_PROTOCOL, writeback=False)
    key = str(url)

    if os.environ["pydna_cache"] in ("compare", "cached"):
        try:
            cached = cache[key]
        except KeyError:
            if os.environ["pydna_cache"] == "compare":
                raise Exception("no result for this key!")
            else:
                refresh = True

    if refresh or os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):
        result = requests.get(url).text
    if os.environ["pydna_cache"] == "compare":
        if result!=cached:
            module_logger.warning('download error')

    if refresh or os.environ["pydna_cache"] == "refresh":
        cache = shelve.open(os.path.join(os.environ["pydna_data_dir"],"genbank"), protocol=pickle.HIGHEST_PROTOCOL, writeback=False)
        cache[key] = result

    elif cached and os.environ["pydna_cache"] not in ("nocache", "refresh"):
        result = cached

    cache.close()
    
    result = textwrap.dedent(result).strip()
    result = result.replace( '\r\n', '\n')
    result = result.replace( '\r',   '\n')
    return result

email = os.getenv("pydna_email")

def genbank(accession):
    gb = Genbank(email)
    return gb.nucleotide(accession)

if __name__=="__main__":
    import doctest
    doctest.testmod()