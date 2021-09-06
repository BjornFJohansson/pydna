#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""docstring."""


class FakeSeq(object):
    """docstring."""

    def __init__(self,
                 length,
                 n=50e-15,  # 50 fmol = 0.05 pmol
                 rf=0):
        self._length = int(length)
        self.n = n
        self.rf = rf

    def m(self):
        """Mass of the DNA molecule in grams."""
        # M(Da) * n (mol) = g
        return self.M() * self.n

    def M(self):
        """M grams/mol."""
        return (308.9 * self._length + 79.0) * 2

    def __len__(self):
        """docstring."""
        return self._length

    def __lt__(self, other):
        """docstring."""
        return self._length < len(other)

    def __repr__(self):
        """docstring."""
        return f"FakeSeq({self._length:.1e})"

    def __str__(self):
        """docstring."""
        return self.__repr__()


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
