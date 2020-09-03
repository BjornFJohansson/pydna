#!/usr/bin/env python
# -*- coding: utf-8 -*-
import timeit

s = """\
a="aaa"
"""
timeit.timeit(stmt=s, number=10000000)

setup = """\
from pydna.dseqrecord import Dseqrecord
from pydna.assembly import _Fragment
"""

s = """\
a=_Fragment(Dseqrecord("aaa"))
"""
timeit.timeit(stmt=s, setup=setup, number=10000000)

abcde = " 234"


print(abcde)
