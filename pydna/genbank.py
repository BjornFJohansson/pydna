#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''Provides a class for downloading sequences from genbank.
'''
import pickle as _pickle
import shelve as _shelve
import re as _re
import os as _os
import logging as _logging
_module_logger = _logging.getLogger("pydna."+__name__)

from Bio import Entrez as _Entrez
from .readers import read as _read
from .genbankrecord import GenbankRecord as _GenbankRecord

class Genbank(object):
    '''Class to facilitate download from genbank.

    Parameters
    ----------
    users_email : string
        Has to be a valid email address. You should always tell
        Genbanks who you are, so that they can contact you.

    Examples
    --------

    >>> import pydna                                                           
    >>> gb=pydna.Genbank("bjornjobb@gmail.com")                               
    >>> rec = gb.nucleotide("AJ515731")                 # <- entry from genbank                  
    >>> print(len(rec))                                                        
    19
    '''

    def __init__(self, users_email):

        if not _re.match("[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}",users_email,_re.IGNORECASE):
            raise ValueError
        if users_email == "someone@example.com":
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

        if _os.environ["pydna_cache"] in ("compare", "cached"):
            cache = _shelve.open(_os.path.join(_os.environ["pydna_data_dir"], "genbank"),  
                                protocol=_pickle.HIGHEST_PROTOCOL, 
                                writeback=False)
            key = item+str(start)+str(stop)+str(strand)
            try:
                cached, item, start, stop = cache[key]
            except:
                if _os.environ["pydna_cache"] == "compare":
                    raise Exception("no result for this key!")
                else:
                    refresh = True

        if refresh or _os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):

            matches =((1, _re.search("(REGION:\s(?P<start>\d+)\.\.(?P<stop>\d+))", item)),
                      (2, _re.search("(REGION: complement\((?P<start>\d+)\.\.(?P<stop>\d+)\))",item)),
                      (1, _re.search(":(?P<start>\d+)-(?P<stop>\d+)",item)),
                      (2, _re.search(":c(?P<start>\d+)-(?P<stop>\d+)",item)),
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

            _Entrez.email = self.email

            text = _Entrez.efetch( db        ="nucleotide",
                                  id        = item,
                                  rettype   = "gbwithparts",
                                  seq_start = start,
                                  seq_stop  = stop,
                                  strand    = strand,
                                  retmode   = "text" ).read()
            dsr = _read(text)
            result = _GenbankRecord(dsr, item = item, start=start, stop=stop, strand=strand)
    
        if _os.environ["pydna_cache"] == "compare":
            if result!=cached:
                _module_logger.warning('download error')

        if refresh or _os.environ["pydna_cache"] == "refresh":
            cache = _shelve.open(_os.path.join(_os.environ["pydna_data_dir"], "genbank"), protocol=_pickle.HIGHEST_PROTOCOL, writeback=False)
            cache[key] = result, item, start, stop
        elif cached and _os.environ["pydna_cache"] not in ("nocache", "refresh"):
            result = cached
            cache.close()

        return result

        
email = _os.getenv("pydna_email")

def genbank(accession):
    gb = Genbank(email)
    return gb.nucleotide(accession)
