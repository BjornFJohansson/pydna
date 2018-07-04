#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

def test_myprimers(monkeypatch):
    monkeypatch.setenv("pydna_primers", "primers_linux_line_endings.txt")    
    from pydna import myprimers
    from pydna.parsers import parse_primers 
    from importlib import reload
    reload(myprimers)
    newlist = parse_primers("primers_linux_line_endings.txt")[::-1]    
    primerdict = myprimers.dict_primers
    
    assert len(primerdict) == 4
    primer_list = myprimers.list_primers
    assert primer_list==newlist
    
    with pytest.raises(NotImplementedError):
        myprimers.append_primer_list("slask")       

if __name__ == '__main__':
    pytest.main([__file__, "-v", "-s", "--cov=pydna","--cov-report=html"])
