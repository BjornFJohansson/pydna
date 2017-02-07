#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from pydna._pretty import pretty_str

def test_pretty():

    import io
    from Bio import SeqIO
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

    with open("pydna_read_test.txt", 'r', encoding="utf-8") as f:
        raw = f.read()

    handle = io.StringIO(raw)

    sr = SeqIO.read(handle, "genbank", alphabet=IUPACAmbiguousDNA())

    s = sr.format("gb").strip()

    pretty_label = pretty_str( s[:55]+"circular"+s[63:] )[559:578]

    fe=sr.features[0]
    label_from_sr = fe.qualifiers["label"][0]                                   
    
    assert raw[295:312] == '/label=2micron 2µ'

    assert label_from_sr == "2micron 2µ"

    assert s[559:578] == '/label="2micron 2µ"'

    assert pretty_label == '/label="2micron 2µ"'

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
