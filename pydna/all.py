#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''The pretty_str class is same as str but has 
   a _repr_pretty_ method for
   for nicer string output 
   in the IPython shell'''
   
__all__=["Anneal", 
         "pcr", 
         "Assembly", 
         "genbank",
         "Genbank",
         "download_text",  
         "Dseqrecord", 
         "read", 
         "read_primer",  
         "parse", 
         "parse_primers", 
         "ape", 
         "primer_design", 
         "assembly_fragments",  
         "eq", 
         "shift_origin", 
         "pairwise", 
         "gbtext_clean", 
         "list_primers"]


from pydna.amplify    import Anneal
from pydna.amplify    import pcr
from pydna.assembly   import Assembly
from pydna.genbank    import genbank
from pydna.genbank    import Genbank
from pydna.download   import download_text
from pydna.dseqrecord import Dseqrecord
from pydna.readers    import read
from pydna.readers    import read_primer
from pydna.parsers    import parse
from pydna.parsers    import parse_primers
from pydna.editor     import ape
from pydna.design     import primer_design
from pydna.design     import assembly_fragments
from pydna.utils      import eq
from pydna.utils      import shift_origin
from pydna.utils              import pairwise
from pydna.genbankfixer       import gbtext_clean
from pydna.myprimers   import list_primers

if __name__=="__main__":
    import os as _os
    cache = _os.getenv("pydna_cache", "nocache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
