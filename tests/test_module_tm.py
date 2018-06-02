#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

def test_tms():
    from pydna import tm
    args= []
    kwargs = {"key":"value"}
    
    primer = "AGTCTAGTCTGTGTAGTTTCGACTAGTCTATCG"
    primer = "agatcgactatctatcttatgcactatgtctat"
    
    assert 55.08264292462792 == tm.tmstaluc98(primer,*args, dnac=50, saltc=50, **kwargs)
    
    assert 63.86257335124958 == tm.tmbreslauer86(primer, *args, dnac=500.0, saltc=50, thermodynamics=False, **kwargs)
    assert (63.86257335124958, 229.4, 650.4777432533716) == tm.tmbreslauer86(primer, *args, dnac=500.0, saltc=50, thermodynamics=True, **kwargs)
    
    assert 63.598585735078075 == tm.tmbresluc(primer, *args, primerc=500.0, saltc=50, thermodynamics=False, **kwargs)
    assert (63.598585735078075, -234400, -615.1999999999998) == tm.tmbresluc(primer, *args, primerc=500.0, saltc=50, thermodynamics=True, **kwargs)    
    assert 88 == tm.basictm(primer, *args, **kwargs)
    
    with pytest.raises(NotImplementedError):
        tm.Q5()
    
if __name__ == '__main__':
    pytest.main([__file__, "-vv", "-s", "--cov=pydna","--cov-report=html"])