#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
'''This module provides miscellaneous functions.'''

import shelve      as _shelve
import os          as _os
import re          as _re
import logging     as _logging
import base64      as _base64
import pickle      as _pickle
import hashlib     as _hashlib
import keyword     as _keyword
import collections as _collections
import itertools   as _itertools
_module_logger = _logging.getLogger("pydna."+__name__)

from Bio.SeqUtils.CheckSum  import seguid   as _base64_seguid
from pydna._pretty       import pretty_str       as _pretty_str
from Bio.Seq             import _maketrans
#from Bio.Seq             import reverse_complement as _reverse_complement
from Bio.Data.IUPACData  import ambiguous_dna_complement as _ambiguous_dna_complement
#from Bio.Seq             import reverse_complement as _rc


_ambiguous_dna_complement.update({"U":"A"})
_complement_table = _maketrans(_ambiguous_dna_complement)


def rc(sequence:str):
    '''returns the reverse complement of sequence (string)
    accepts mixed DNA/RNA
    '''
    return sequence.translate(_complement_table)[::-1]

def complement(sequence:str):
    '''returns the complement of sequence (string)
    accepts mixed DNA/RNA
    '''
    return sequence.translate(_complement_table)

def memorize(filename):
    """Decorator for caching fucntions and classes, see pydna.download and pydna.Assembly for use"""
    def decorator(f):
        def wrappee( *args, **kwargs):
            _module_logger.info( "#### memorizer ####" )
            _module_logger.info( "cache filename                   = %s", filename )
            _module_logger.info( "os.environ['pydna_cached_funcs'] = %s", _os.environ["pydna_cached_funcs"] )
            if filename not in _os.environ["pydna_cached_funcs"]:
                _module_logger.info("cache filename not among cached functions, made it new!")
                return f(*args, **kwargs)               
            key = _base64.urlsafe_b64encode(_hashlib.sha1(_pickle.dumps((args, kwargs))).digest()).decode("ascii")
            _module_logger.info( "key = %s", key )
            cache = _shelve.open(_os.path.join(_os.environ["pydna_data_dir"], identifier_from_string(filename)), writeback=False)
            try:
                result = cache[key]
            except KeyError:
                _module_logger.info("no result for key %s in shelve %s", key, identifier_from_string(filename))
                result = f(*args, **kwargs)
                _module_logger.info("made it new!")
                cache[key] = result
                _module_logger.info("saved result under key %s", key)
            else:
                _module_logger.info( "found %s in cache", key)
            cache.close()
            return result
        return wrappee
    return decorator        


def eq(*args,**kwargs):
    '''Compares two or more DNA sequences for equality i.e. they
    represent the same DNA molecule. Comparisons are case insensitive.

    Parameters
    ----------
    args : iterable
        iterable containing sequences
        args can be strings, Biopython Seq or SeqRecord, Dseqrecord
        or dsDNA objects.
    circular : bool, optional
        Consider all molecules circular or linear
    linear : bool, optional
        Consider all molecules circular or linear

    Returns
    -------
    eq : bool
        Returns True or False

    Notes
    -----

    Compares two or more DNA sequences for equality i.e. if they
    represent the same DNA molecule.

    Two linear sequences are considiered equal if either:

    * They have the same sequence (case insensitive)
    * One sequence is the reverse complement of the other (case insensitive)

    Two circular sequences are considered equal if they are circular permutations:

    1. They have the same lengt, AND
    2. One sequence or can be found in the concatenation of the other sequence with it    , OR
    3. The reverse complement can be found in the concatenation of the other sequence with itself.

    The topology for the comparison can be set using one of the keywords
    linear or circular to True or False.

    If circular or linear is not set, it will be deduced from the topology of
    each sequence for sequences that have a linear or circular attribute
    (like Dseq and Dseqrecord).

    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.utils import eq
    >>> eq("aaa","AAA")
    True
    >>> eq("aaa","AAA","TTT")
    True
    >>> eq("aaa","AAA","TTT","tTt")
    True
    >>> eq("aaa","AAA","TTT","tTt", linear=True)
    True
    >>> eq("Taaa","aTaa", linear = True)
    False
    >>> eq("Taaa","aTaa", circular = True)
    True
    >>> a=Dseqrecord("Taaa")
    >>> b=Dseqrecord("aTaa")
    >>> eq(a,b)
    False
    >>> eq(a,b,circular=True)
    True
    >>> a=a.looped()
    >>> b=b.looped()
    >>> eq(a,b)
    True
    >>> eq(a,b,circular=False)
    False
    >>> eq(a,b,linear=True)
    False
    >>> eq(a,b,linear=False)
    True
    >>> eq("ggatcc","GGATCC")
    True
    >>> eq("ggatcca","GGATCCa")
    True
    >>> eq("ggatcca","tGGATCC")
    True


    '''

    args = flatten(args) # flatten
    
    topology = None

    if "linear" in kwargs:
        if kwargs["linear"]==True:
            topology = "linear"
        if kwargs["linear"]==False:
            topology = "circular"
    elif "circular" in kwargs:
        if kwargs["circular"]==True:
            topology = "circular"
        if kwargs["circular"]==False:
            topology = "linear"
    else:
        # topology keyword not set, look for topology associated to each sequence
        # otherwise raise exception
        topology = set([arg.circular if hasattr(arg, "circular") else None for arg in args])

        if len(topology)!=1:
            raise ValueError("sequences have different topologies")
        topology = topology.pop()
        if topology in (False, None):
            topology = "linear"
        elif topology==True:
            topology = "circular"

    #args_string_list    = [str(arg.seq).lower() if hasattr(arg,"seq") else str(arg).lower() for arg in args]

    args = [arg.seq if hasattr(arg, "seq") else arg for arg in args]
    args_string_list    = [arg.watson.lower() if hasattr(arg, "watson") else str(arg).lower() for arg in args]

    length = set((len(s) for s in args_string_list))

    if len(length)!=1:
        return False
    same = True

    if topology == "circular":
        # force circular comparison of all given sequences
        for s1, s2 in _itertools.combinations(args_string_list, 2):
            if not ( s1 in s2+s2 or rc(s1) in s2+s2):
                same = False
    elif topology == "linear":
        # force linear comparison of all given sequences
        for s1,s2 in _itertools.combinations(args_string_list, 2):
            if not ( s1==s2 or s1==rc(s2) ):
                same = False
    return same


def SmallestRotation(s):
    prev,rep = None,0
    ds=2*s
    lends=len(ds)
    old = 0
    k = 0
    w=""
    while k < lends:
        i,j = k,k+1
        while j < lends and ds[i] <= ds[j]:
            i = (ds[i] == ds[j]) and i+1 or k
            j += 1
        while k < i+1:
            k += j-i
            prev=w
            w=ds[old:k]
            old = k
            if w == prev:
                rep += 1
            else:
                prev,rep = w,1
            if len(w)*rep == len(s):
                return w*rep
            
            
#try:
#    import pyximport
#except ImportError:
#    pass
#else:
#    pyximport.install()
#    from pydna._smallest import SmallestRotation
    

def identifier_from_string(s:str) -> str:
    '''This function returns a string that is a valid python identifier based on the argument s or an empty string'''
    s=s.strip()
    s = _re.sub(r"\s+",r"_",s)
    s.replace("-", "_")
    s = _re.sub('[^0-9a-zA-Z_]', '', s)
    if s and not s[0].isidentifier() or _keyword.iskeyword(s):
        s="_{s}".format(s=s)
    assert s=="" or s.isidentifier()
    return s


def seguid(seq: str) -> _pretty_str:
    '''Returns the url safe SEGUID checksum for the sequence. This is the SEGUID
    checksum with the '+' and '/' characters of standard Base64 encoding are respectively
    replaced by '-' and '_'.
    '''
    return _pretty_str( _base64_seguid( seq.upper() ).replace("+","-").replace("/","_") )


def lseguid(seq: str) -> _pretty_str:
    '''Returns the url safe lSEGUID checksum for the sequence (seq). This is the SEGUID
    checksum with the '+' and '/' characters of standard Base64 encoding are respectively
    replaced by '-' and '_'.
    '''
    return seguid( min(seq.upper(), str(rc(seq)).upper() )).replace("+","-").replace("/","_")


def cseguid(seq: str) -> _pretty_str:
    '''Returns the url safe cSEGUID for the sequence. The cSEGUID is the SEGUID checksum
    calculated for the lexicographically minimal string rotation of a DNA sequence.
    Only defined for circular sequences.
    '''
    return seguid( min( SmallestRotation(seq.upper()), SmallestRotation(str(rc(seq)).upper())))


def flatten(*args): # flatten
    """Flattens an iterable of iterables down to str, bytes, bytearray or any of the pydna or Biopython seq objects"""
    output = []
    args=list(args)
    while args:
        top = args.pop()
        if isinstance(top, _collections.Iterable) and not isinstance(top, (str, bytes, bytearray)) and not hasattr(top, "features"):
            args.extend(top)
        else:
            output.append(top)
    return output[::-1]


if __name__=="__main__":
    import os as _os
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"]=""
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"]=cached

