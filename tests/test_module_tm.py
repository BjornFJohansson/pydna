#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_tms():
    from pydna import tm

    # args = []
    # kwargs = {"key": "value"}

    primer = "AGTCTAGTCTGTGTAGTTTCGACTAGTCTATCG"
    # primer = "agatcgactatctatcttatgcactatgtctat"

    # assert 55.08264292462792 == tm.tmstaluc98(primer,*args, dnac=50, saltc=50, **kwargs)
    # assert 63.86257335124958 == tm.tmbreslauer86(primer, *args, dnac=500.0, saltc=50, thermodynamics=False, **kwargs)
    # assert (63.86257335124958, 229.4, 650.4777432533716) == tm.tmbreslauer86(primer, *args, dnac=500.0, saltc=50, thermodynamics=True, **kwargs)
    # assert 63.598585735078075 == tm.tmbresluc(primer, *args, primerc=500.0, saltc=50, thermodynamics=False, **kwargs)
    # assert (63.598585735078075, -234400, -615.1999999999998) == tm.tmbresluc(primer, *args, primerc=500.0, saltc=50, thermodynamics=True, **kwargs)
    # assert 88 == tm.basictm(primer, *args, **kwargs)

    assert tm.tm_default(primer) == pytest.approx(67.78918110181166)
    assert tm.tm_dbd(primer) == pytest.approx(62.74633103079093)
    assert tm.tm_product(primer * 20) == pytest.approx(76.27411419319003)
    assert tm.ta_default(primer, primer, primer * 20) == pytest.approx(
        58.82863426577652
    )
    assert tm.ta_dbd(primer, primer, primer * 20) == pytest.approx(65.74633103079093)
    assert tm.tmbresluc(primer) == pytest.approx(63.38496307044147)

    with pytest.raises(NotImplementedError):
        tm.Q5(primer)


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
