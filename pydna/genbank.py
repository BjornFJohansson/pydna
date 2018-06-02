#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''Provides a class for downloading sequences from genbank.
'''
import re as _re
import os as _os
import logging as _logging
_module_logger = _logging.getLogger("pydna."+__name__)

from Bio import Entrez as _Entrez
from pydna.readers import read as _read
from pydna.genbankrecord import GenbankRecord as _GenbankRecord
from pydna.utils  import memorize as _memorize

#TODO http://httpbin.org/ use for testing?

class Genbank(object):
    '''Class to facilitate download from genbank.

    Parameters
    ----------
    users_email : string
        Has to be a valid email address. You should always tell
        Genbanks who you are, so that they can contact you.

    Examples
    --------

    >>> from pydna.genbank import Genbank
    >>> gb=Genbank("bjornjobb@gmail.com")
    >>> rec = gb.nucleotide("AJ515731")                 # <- entry from genbank
    >>> print(len(rec))
    19
    '''

    def __init__(self, users_email:str, *args, tool="pydna", **kwargs):

        if not _re.match("[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}",users_email,_re.IGNORECASE):
            raise ValueError("email address {} is not valid.".format(users_email))
        if users_email == "someone@example.com":
            raise ValueError("you have to set your email address in order to download from Genbank")
        self.email=users_email
        self.tool = tool

    def __repr__(self):
        return "GenbankConnection({})".format(self.email)

    @_memorize("Genbank_nucleotide")
    def nucleotide(self, item:str, seq_start=None, seq_stop=None, strand=1):
        '''Download a genbank nuclotide record.

        Item is a string containing one genbank acession number.
        for a nucleotide file. Genbank nucleotide accession numbers have this format:

        | A12345   = 1 letter  + 5 numerals
        | AB123456 = 2 letters + 6 numerals

        The accession number is sometimes followed by a point and version number

        | BK006936.2

        Item can also contain optional interval information in the following formats:

        | BK006936.2 REGION: complement(613900..615202)
        | NM_005546 REGION: 1..100
        | NM_005546 REGION: complement(1..100)
        | 21614549:1-100
        | 21614549:c100-1

        Start and stop are the sequence intervals to be downloaded. This
        is useful for large genbank records.
        If strand is 2. "c", "C", "crick", "Crick", "antisense","Antisense",
        "2" or 2, the antisense (Crick) strand is returned, otherwise
        the sense (Watson) strand is returned.

        Result is returned as a GenbankRecord object.

        References
        ----------

        .. [#]   http://www.dsimb.inserm.fr/~fuchs/M2BI/AnalSeq/Annexes/Sequences/Accession_Numbers.htm
        .. [#]   http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
        '''
        matches =((1, _re.search("(REGION:\s(?P<start>\d+)\.\.(?P<stop>\d+))", item)),
                  (2, _re.search("(REGION: complement\((?P<start>\d+)\.\.(?P<stop>\d+)\))",item)),
                  (1, _re.search(":(?P<start>\d+)-(?P<stop>\d+)",item)),
                  (2, _re.search(":c(?P<start>\d+)-(?P<stop>\d+)",item)))

        for strand_, match in matches:
            if match:
                seq_start = match.group("start")
                seq_stop  = match.group("stop")
                item = item[:match.start()]
                strand=strand_
                break
        

        if not strand in [1,2]:
            try:
                strand = {"c":2, "crick":2, "antisense":2, "2":2}[strand.lower()]
            except (KeyError,AttributeError):
                strand = 1


        _module_logger.info("#### Genbank download ####")
        _module_logger.info("item  %s", item)
        _module_logger.info("start %s", seq_start)
        _module_logger.info("stop  %s", seq_stop)

        _module_logger.info("strand  %s", str(strand))

        _Entrez.email = self.email
        _Entrez.tool  = self.tool
        
        _module_logger.info("Entrez.email  %s", self.email)
        text = _Entrez.efetch(db        ="nuccore",
                              id        = item,
                              rettype   = "gbwithparts",
                              seq_start = seq_start,
                              seq_stop  = seq_stop,
                              strand    = strand,
                              retmode   = "text" ).read()

        _module_logger.info("text[:160]  %s", text[:160])

        return _GenbankRecord(_read(text), item = item, start=seq_start, stop=seq_stop, strand=strand)


def genbank(accession:str="CS570233.1"):
    '''Download a genbank nuclotide record.'''
    email = _os.getenv("pydna_email")
    gb = Genbank(email)
    return gb.nucleotide(accession)


if __name__=="__main__":
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"]=""   
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"]=cached
    pass

    
    
    