#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import nose
import pydna

def test_download():
    a = pydna.Genbank("bjornjobb@gmail.com")
    v = os.environ["pydna_cache"]
    os.environ["pydna_cache"] = "nocache"
    assert a.test()
    os.environ["pydna_cache"] = v


if __name__ == '__main__':
    nose.runmodule()

