#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
_pretty tests
'''
import nose, sys

import pydna

def test_pretty():

    x= '/label="2micron 2µ"'

    import io
    from Bio import SeqIO
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

    raw = open("tests/pydna_read_test.txt", 'r').read()

    handle = io.StringIO(raw)

    sr = SeqIO.read(handle, "genbank", alphabet=IUPACAmbiguousDNA())

    s = sr.format("gb").strip()

    assert s[559:578] == x

    y = pydna._pretty.pretty_string( s[:55]+"circular"+s[63:] )[559:578]
    
    assert x==y

    assert pydna.read("tests/pydna_read_test.txt").format("gb")[559:578] == x

    print(x)
    print(y)


if __name__ == '__main__':
    nose.runmodule(argv=[sys.argv[0], '--nocapture'])
