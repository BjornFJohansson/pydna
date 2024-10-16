#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_user():
    from pydna.user import USER
    from pydna.dseq import Dseq

    a = Dseq.from_representation(
        """\
    AACGAuGTCGACTTAGATCTCACAGGCTTTTTTCAAGaCGGCCTTGAATTCAGTCATTTGGATCCGGCCGATCTTT
    TTGCTACAGCTGAATCTAGAGTGTCCGAAAAAAGTTCTGCCGGAACTTAAGTCAGTAAACCTAGGCCGGCuAGAAA
    """
    )

    result_a = Dseq.from_representation(
        """\
          GTCGACTTAGATCTCACAGGCTTTTTTCAAGaCGGCCTTGAATTCAGTCATTTGGATCCGGCCGATCTTT
    TTGCTACAGCTGAATCTAGAGTGTCCGAAAAAAGTTCTGCCGGAACTTAAGTCAGTAAACCTAGGCCGGC
    """
    )

    print(result_a)

    us = USER(size=4)
    print(us.products(a)[0])
    assert us.products(a)[0] == result_a


def test_many_user_sites():
    pass


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
