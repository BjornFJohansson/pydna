#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""docstring."""

from pydna.fakeseq import FakeSeq

PennStateLadder = [FakeSeq(int(n)) for n in (
                   10000, 7750, 5000, 4000,
                   3000, 2000, 1500, 1000,
                   750, 500)]

GeneRuler_1kb_ = [FakeSeq(int(n)) for n in (
                  10000, 8000, 6000, 5000, 4000, 3500, 3000,
                  2500, 2000, 1500, 1000, 750, 500, 250,)]

GeneRuler_1kb_plus = [FakeSeq(int(n)) for n in (
                      20000, 10000, 7000, 5000, 4000, 3000, 2000,
                      1500, 1000, 700, 500, 400, 300, 200, 75,)]

GeneRuler_Mix = [FakeSeq(int(n)) for n in (
                 10000, 8000, 6000, 5000, 4000, 3500, 3000,
                 2500, 2000, 1500, 1200, 1000, 900, 800, 700,
                 600, 500, 400, 300, 200, 100,)]


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
