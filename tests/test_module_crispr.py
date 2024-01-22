#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_crispr():
    from pydna.crispr import cas9, protospacer
    from pydna.dseqrecord import Dseqrecord
    from pydna.dseq import Dseq

    a = Dseq.from_representation(
        """\
    GTTACTTTACCCGACGT
    CAATGAAATGGGCTGCA
    """
    )

    b = Dseq.from_representation(
        """\
    CCCaGG
    GGGtCC
    """
    )

    containing_sgRNA = Dseqrecord("GTTACTTTACCCGACGTCCCgttttagagctagaaatagcaagttaaaataagg")

    target = Dseqrecord("GTTACTTTACCCGACGTCCCaGG")

    assert [f.seq for f in target.cut([cas9(ps) for ps in protospacer(containing_sgRNA)])] == [a, b]
    assert [f.seq for f in target.cut([cas9(ps) for ps in protospacer(containing_sgRNA.rc())])] == [a, b]
    assert [f.seq for f in target.rc().cut([cas9(ps) for ps in protospacer(containing_sgRNA)])] == [b.rc(), a.rc()]
    assert [f.seq for f in target.rc().cut([cas9(ps) for ps in protospacer(containing_sgRNA.rc())])] == [
        b.rc(),
        a.rc(),
    ]


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
