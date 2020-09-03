#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""This module establish a RestrictionBatch based on enzymes found in a text file specified in the enzymes entry 
in the python.ini file or by the environment variable pydna_enzymes.

The text file will be searched for all enzymes in the biopython 
AllEnzymes batch which is located in the Bio.Restriction package.

The pydna.myenzymes.myenzymes contains a new restriction batch with the enzymes contained
within the file specified.
"""

import os as _os
from Bio.Restriction import AllEnzymes as _AllEnzymes
from Bio.Restriction import RestrictionBatch as _RestrictionBatch
import logging as _logging
import traceback as _traceback

_module_logger = _logging.getLogger("pydna." + __name__)


_text = ""


try:
    with open(_os.environ["pydna_enzymes"], encoding="utf-8") as _f:
        _text = _f.read()
except FileNotFoundError:
    _module_logger.warning("%s not found.", _os.environ["pydna_enzymes"])
except IsADirectoryError:
    _module_logger.warning("%s is a directory.", _os.environ["pydna_enzymes"])
except IOError:
    _module_logger.warning(
        "%s found, but could not be read.", _os.environ["pydna_enzymes"]
    )
except Exception as e:
    _module_logger.warning(_traceback.format_exc())

myenzymes = _RestrictionBatch(
    [e for e in _AllEnzymes if str(e).lower() in _text.lower()]
)

if __name__ == "__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"] = cache
