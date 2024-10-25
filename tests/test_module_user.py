#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_user_enzyme():
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

    us = USER(size=4)
    # print(a.cut(us)[1].__repr__())
    # print(result_a)
    assert a.cut(us)[1] == result_a


def test_user_enzyme_short_motif():
    """
    A user site without enough upstream bases should not be recongnized.
    """
    from pydna.user import USER
    from pydna.dseq import Dseq

    a = Dseq.from_representation(
        """\
    AACGAuGTCGACTTAGATCTCACAGGCTTTTTTCAAGaCGGCCTTGAATTCAGTCATTTGGATCCGGCCGAT
    TTGCTACAGCTGAATCTAGAGTGTCCGAAAAAAGTTCTGCCGGAACTTAAGTCAGTAAACCTAGGCCGGCuA
    """
    )

    result_a = Dseq.from_representation(
        """\
          GTCGACTTAGATCTCACAGGCTTTTTTCAAGaCGGCCTTGAATTCAGTCATTTGGATCCGGCCGAT
    TTGCTACAGCTGAATCTAGAGTGTCCGAAAAAAGTTCTGCCGGAACTTAAGTCAGTAAACCTAGGCCGGCuA
    """
    )

    us = USER(size=4)
    # print(a.cut(us)[1].__repr__())
    # print(a.cut(us)[1].reverse_complement())
    # print(result_a)
    assert a.cut(us)[1] == result_a


def test_many_user_sites():
    pass


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
