#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""This module ...
"""

from glob import glob as _glob
import os as _os
import logging as _logging
import pathlib as _pathlib

_module_logger = _logging.getLogger("pydna." + __name__)

seqdir = _pathlib.Path(_os.environ["pydna_data_dir"])

mysequences = {
    u.name: u
    for u in (
        _pathlib.Path(p) for p in _glob(str(seqdir.joinpath("**/*.gb")), recursive=True)
    )
}

if __name__ == "__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"] = cache
