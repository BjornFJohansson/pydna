#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import nose
import pydna

from pydna.utils import seguid

def test_download():
    a = pydna.Genbank("bjornjobb@gmail.com")
    v = os.environ["pydna_cache"]
    os.environ["pydna_cache"] = "nocache"
    gb = pydna.Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("L09137")
    assert len(result) == 2686
    assert seguid(result.seq) == "71B4PwSgBZ3htFjJXwHPxtUIPYE"
    os.environ["pydna_cache"] = v

if __name__ == '__main__':
    nose.runmodule()

