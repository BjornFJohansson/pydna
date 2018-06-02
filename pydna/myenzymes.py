#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''This module establish a RestrictionBatch based on enzymes found in a text file specified in the python.ini file
   or by the enviromnment variable 

'''
import os as _os
from Bio.Restriction import AllEnzymes as _AllEnzymes
from Bio.Restriction import RestrictionBatch as _RestrictionBatch
import logging   as _logging
import traceback as _traceback
_module_logger = _logging.getLogger("pydna."+__name__)


_text=""   # 894

try:
    with open( _os.environ["pydna_enzymes"], encoding="utf-8") as _f:
        _text = _f.read()
except FileNotFoundError:
    _module_logger.warning("%s not found.", _os.environ["pydna_enzymes"])
except IsADirectoryError:
    _module_logger.warning("%s is a directory.", _os.environ["pydna_enzymes"])
except IOError:
    _module_logger.warning("%s found, but could not be read.", _os.environ["pydna_enzymes"])
except Exception as e:
   _module_logger.warning(_traceback.format_exc())
    
myenzymes = _RestrictionBatch([e for e in _AllEnzymes if str(e).lower() in _text.lower()])

if __name__=="__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
