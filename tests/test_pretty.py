#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
_pretty tests
'''
import sys
import pytest
import pydna

def test_pretty():

    print(sys.getdefaultencoding())
    import locale
    print(locale.getpreferredencoding())

    import io
    from Bio import SeqIO
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

    raw = open("pydna_read_test.txt", 'r', encoding="UTF8").read()

    handle = io.StringIO(raw)

    sr = SeqIO.read(handle, "genbank", alphabet=IUPACAmbiguousDNA())

    s = sr.format("gb").strip()

    assert s[559:578] == '/label="2micron 2µ"'

    y = pydna._pretty.pretty_string( s[:55]+"circular"+s[63:] )[559:578]
    
    assert '/label="2micron 2µ"'==y

    assert pydna.read("pydna_read_test.txt").format("gb")[559:578] == '/label="2micron 2µ"'

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
