#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from pydna.dseq import Dseq

def test_dseqs():
    ''' test Dseq dsdata and todata properties'''
    a = Dseq("ccGGATCC", "aaggatcc", -2)
    print()
    print(repr(a))
    assert a.todata == "ccGGATCCtt"
    assert a.dsdata == "GGATCC"
    print()
    
    a2 = Dseq("ccGGATCCaa", "ggatcc", -2)
    print(repr(a2))
    assert a2.todata  == "ccGGATCCaa"
    assert a2.dsdata == "GGATCC"
    print()
    a3 = Dseq("ccGGATCC", "ggatcc", -2)
    print(repr(a3))
    assert a3.todata == "ccGGATCC"
    assert a3.dsdata == "GGATCC"
    print()
    
    
    
    b = Dseq("GGATCC", "aaggatcccc", 2)
    print(repr(b))
    assert b.todata == "ggGGATCCtt"
    assert b.dsdata == "GGATCC"
    print()
    b2 = Dseq("GGATCCaa", "ggatcccc", 2)
    print(repr(b2))
    assert b2.todata == "ggGGATCCaa"
    assert b2.dsdata == "GGATCC"
    print()
    b3 = Dseq("GGATCC", "ggatcccc", 2)
    print(repr(b3))
    assert b3.todata == "ggGGATCC"
    assert b3.dsdata == "GGATCC"
    print()
    
    
    c = Dseq("GGATCCaaa", "ggatcc", 0)
    print(repr(c))
    assert c.todata == "GGATCCaaa"
    assert c.dsdata ==  "GGATCC"
    print()
    d = Dseq("GGATCC", "aaaggatcc", 0)
    print(repr(d))
    assert d.todata == "GGATCCttt"
    assert d.dsdata == "GGATCC"



if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])

