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

    user = USER()
    assert a.cut(user)[1] == result_a


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
