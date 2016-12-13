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
import os
import re
import io
import textwrap
import glob

from Bio                    import SeqIO
from Bio.Alphabet.IUPAC     import IUPACAmbiguousDNA
from Bio.GenBank            import RecordParser
from .genbankfile            import GenbankFile
from .dseqrecord             import Dseqrecord
from .primer                 import Primer  

def parse2(data, ds = True):
    '''experimental'''

    pattern =  r"(?:>.+\n^(?:^[^>]+?)(?=\n\n|>|LOCUS|ID))|(?:(?:LOCUS|ID)(?:(?:.|\n)+?)^//)"

    def extract_seqs(raw):
        raw = textwrap.dedent(raw).strip()
        raw = raw.replace( '\r\n', '\n')
        raw = raw.replace( '\r',   '\n')
        return re.findall(pattern, textwrap.dedent(raw+ "\n\n"),re.MULTILINE)

    files=[]
    rawseqs=[]

    if not hasattr(data, '__iter__'):
        data = (data,)
    for item in data:
        for pth in glob.glob(item):
            if os.path.isfile(pth):
                files.append(os.path.abspath(pth))
            else:
                for dirpath,_,filenames in os.walk(pth):
                    for f in filenames:
                        files.append( os.path.abspath(os.path.join(dirpath, f)))
            for file_ in files:
                with open(file_,'r') as f:
                    rawseqs.extend(extract_seqs(f.read()))
            files=[]
        else:
            rawseqs.extend(extract_seqs(item))
    sequences = []
    while rawseqs:
        circular = False
        rawseq = rawseqs.pop(0)
        handle = io.StringIO(rawseq)
        try:
            parsed = SeqIO.read(handle, "embl", alphabet=IUPACAmbiguousDNA())
            #original_format = "embl"
            if "circular" in rawseq.splitlines()[0]:
                circular = True
        except ValueError:
            handle.seek(0)
            try:
                parsed = SeqIO.read(handle, "genbank", alphabet=IUPACAmbiguousDNA())
                #original_format = "genbank"
                handle.seek(0)
                parser = RecordParser()
                residue_type = parser.parse(handle).residue_type
                if "circular" in residue_type:
                    circular = True
            except ValueError:
                handle.seek(0)
                try:
                    parsed = SeqIO.read(handle, "fasta", alphabet=IUPACAmbiguousDNA())
                    if "circular" in rawseq.splitlines()[0]:
                        circular = True
                except ValueError:
                    continue

        if ds:
            sequences.append(Dseqrecord(parsed, circular = circular))
        else:
            sequences.append(parsed)
        handle.close()
    return sequences

def parse(data, ds = True):
    '''This function returns *all* DNA sequences found in data. If no
    sequences are found, an empty list is returned. This is a greedy
    function, use carefully.

    Parameters
    ----------
    data : string or iterable
        The data parameter is a string containing:

        1. an absolute path to a local file.
           The file will be read in text
           mode and parsed for EMBL, FASTA
           and Genbank sequences.

        2. a string containing one or more
           sequences in EMBL, GENBANK,
           or FASTA format. Mixed formats
           are allowed.

        3. data can be a list or other iterable where the elements are 1 or 2

    ds : bool
        If True double stranded :class:`Dseqrecord` objects are returned.
        If False single stranded :class:`Bio.SeqRecord` [#]_ objects are returned.

    Returns
    -------
    list
        contains Dseqrecord or SeqRecord objects

    References
    ----------

    .. [#] http://biopython.org/wiki/SeqRecord

    See Also
    --------
    read

    '''
    
    def embl_gb_fasta(raw, ds, path=None):
        
        pattern =  r"(?:>.+\n^(?:^[^>]+?)(?=\n\n|>|LOCUS|ID))|(?:(?:LOCUS|ID)(?:(?:.|\n)+?)^//)"
        
        result_list = []
        
        rawseqs = re.findall(pattern, textwrap.dedent(raw + "\n\n"), flags=re.MULTILINE)
        
        for rawseq in rawseqs:
            handle = io.StringIO(rawseq)
            circular=False
            try:
                parsed = SeqIO.read(handle, "embl", alphabet=IUPACAmbiguousDNA())
                if "circular" in rawseq.splitlines()[0]:
                    circular = True
            except ValueError:
                handle.seek(0)
                try:
                    parsed = SeqIO.read(handle, "genbank", alphabet=IUPACAmbiguousDNA())
                    handle.seek(0)
                    parser = RecordParser()
                    residue_type = parser.parse(handle).residue_type
                    if "circular" in residue_type:
                        circular = True
                except ValueError:
                    handle.seek(0)
                    try:
                        parsed = SeqIO.read(handle, "fasta", alphabet=IUPACAmbiguousDNA())
                        if "circular" in rawseq.splitlines()[0]:
                            circular = True
                    except ValueError:
                        parsed = ""
            handle.close()
            
            if ds and path:
                result_list.append( GenbankFile(parsed, circular=circular, path=path) )
            elif ds:
                result_list.append( Dseqrecord(parsed, circular =circular) )
            else:
                result_list.append( parsed )
                
        return result_list

    # a string is an iterable datatype but on Python2.x it doesn't have an __iter__ method.
    if not hasattr(data, '__iter__') or isinstance(data, (str, bytes)):
        data = (data,)
        
    sequences = []
    
    for item in data:
        try:
            # item is a path to a utf-8 encoded text file?
            with open(item, 'r', encoding="utf-8") as f:    
                raw = f.read()  
        except IOError:
            # item was not a path, add sequences parsed from item
            raw=item
            path=None
        else:   
            # item was a readable text file, seqences are parsed from the file
            path = item
        finally:
            sequences.extend( embl_gb_fasta(raw, ds, path) )
    return sequences
    
def parse_primers(data):
    return [Primer(x) for x in parse(data, ds=False)]

    
if __name__=="__main__":
    import os
    cache = os.getenv("pydna_cache")
    os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod()
    os.environ["pydna_cache"]=cache
