#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""docstring."""

from pydna.fakeseq import FakeSeq as _FakeSeq

PennStateLadder = [_FakeSeq(int(n)) for n in (
                   10000,
                   7750,
                   5000,
                   4000,
                   3000,
                   2000,
                   1500,
                   1000,
                   750,
                   500)]


GeneRuler_1kb  = [_FakeSeq(int(n)) for n in (
                  10000,
                  8000,
                  6000,
                  5000,
                  4000,
                  3500,
                  3000,
                  2500,
                  2000,
                  1500,
                  1000,
                  750,
                  500,
                  250,)]


GeneRuler_1kb_plus = [_FakeSeq(ln, n=n*1e-15, rf=rf) for ln, n, rf in (
                      (20000, 1.538,  0.000),  # ( bp length, fmol, Rf )
                      (10000, 3.077,  0.040),
                      (7000,  4.396,  0.096),
                      (5000,  23.077, 0.154),
                      (4000,  7.692,  0.201),
                      (3000,  10.256, 0.261),
                      (2000,  15.385, 0.362),
                      (1500,  82.051, 0.443),
                      (1000,  38.462, 0.562),
                      (700,   54.945, 0.667),
                      (500,  230.769, 0.755),
                      (400,   96.154, 0.808),
                      (300,  128.205, 0.866),
                      (200,  192.308, 0.931),
                      (75,   512.821, 1.000),)]





# GeneRuler_Mix = [_FakeSeq(int(ln),) for ln,n in (
#                  10000,
#                  8000,
#                  6000,
#                  5000,
#                  4000,
#                  3500,
#                  3000,
#                  2500,
#                  2000,
#                  1500,
#                  1200,
#                  1000,
#                  900,
#                  800,
#                  700,
#                  600,
#                  500,
#                  400,
#                  300,
#                  200,
#                  100,)]


# The average weight of a single DNA base pair (bp) is 650 daltons. This can also be written as 650 g/mol (= molar mass)


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
