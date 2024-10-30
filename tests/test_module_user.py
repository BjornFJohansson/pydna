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


def test_user_enzyme2():
    from pydna.user import USER
    from pydna.dseq import Dseq

    a = Dseq.from_representation(
        """\
    AACGAuGTCGACTTAGATCTCACAGGCTTTTTTCAAGaCGGCCTTGAATTCAGTCATTTGGATCCGGCCGATCTTT
    TTGCTACAGCTGAATCTAGAGTGTCCGAAAAAAGTTCuGCCGGAACTTAAGTCAGTAAACCTAGGCCGGCTAGAAA
    """
    )

    result_a = Dseq.from_representation(
        """\
          GTCGACTTAGATCTCACAGGCTTTTTTCAAGaCGGCCTTGAATTCAGTCATTTGGATCCGGCCGATCTTT
    TTGCTACAGCTGAATCTAGAGTGTCCGAAAAAAGTTC
    """
    )

    us = USER(size=4)
    # print(a.cut(us)[1].__repr__())
    # print(result_a)
    assert a.cut(us)[1] == result_a


def test_user_enzyme_impossible():
    from pydna.user import USER
    from pydna.dseq import Dseq

    a = Dseq.from_representation(
        """\
    AACGATGTCGACTTAGATCTCACAGGCTTTTTTCAAGACGGCCTTGAATTCAGTCAuTTGGATCCGGCCGATCTTT
    TTGCTACAGCTGAAuCTAGAGTGTCCGAAAAAAGTTCTGCCGGAACTTAAGTCAGTAAACCTAGGCCGGCTAGAAA
    """
    )

    # These are the three expected results in vitro
    # result_a_1 is the most likely, since the shortest leftover sequences are more likely to detach
    result_a_1 = Dseq.from_representation(
        """\
    AACGATGTCGACTTAGATCTCACAGGCTTTTTTCAAGACGGCCTTGAATTCAGTCA
                   CTAGAGTGTCCGAAAAAAGTTCTGCCGGAACTTAAGTCAGTAAACCTAGGCCGGCTAGAAA
    """
    )
    result_a_2 = Dseq.from_representation(
        """\
                                              TTGGATCCGGCCGATCTTT
    CTAGAGTGTCCGAAAAAAGTTCTGCCGGAACTTAAGTCAGTAAACCTAGGCCGGCTAGAAA
    """
    )

    result_a_3 = Dseq.from_representation(
        """\
    AACGATGTCGACTTAGATCTCACAGGCTTTTTTCAAGACGGCCTTGAATTCAGTCA
    TTGCTACAGCTGAA
    """
    )

    us = USER(size=4)
    for t in [result_a_1, result_a_2, result_a_3]:
        assert t in a.cut(us)


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


def test_nickase_upstream():
    from pydna.user import NtBbvCI
    from pydna.dseq import Dseq

    a = Dseq.from_representation(
        """\
        GCTGAGGCTTAATTAAACCATCAGC
        CGACTCCGAATTAATTTGGTAGTCG
        """
    )

    result_a = Dseq.from_representation(
        """\
        GCTGAGGCTTAATTAAACCATCAGC
        CGACT
        """
    )

    nickase = NtBbvCI()
    assert a.cut(nickase)[0] == result_a


def test_nickase_downstream():
    from pydna.user import NtBbvCI
    from pydna.dseq import Dseq

    a = Dseq.from_representation(
        """\
        GCTGAGTGCTTAATTAAACCTCAGC
        CGACTCACGAATTAATTTGGAGTCG
        """
    )

    result_a = Dseq.from_representation(
        """\
                            TCAGC
        CGACTCACGAATTAATTTGGAGTCG
        """
    )

    nickase = NtBbvCI()
    assert a.cut(nickase)[1] == result_a


def test_nickase_with_crick_site_before_watson_site():
    from pydna.user import NtBbvCI
    from pydna.dseq import Dseq

    a = Dseq.from_representation(
        """\
        GCTGAGGCTTAATTAAACCTCAGC
        CGACTCCGAATTAATTTGGAGTCG
        """
    )

    result_a_1 = Dseq.from_representation(
        """\
        GCTGAGGCTTAATTAAA
               GAATTAATTTGGAGTCG
        """
    )

    result_a_2 = Dseq.from_representation(
        """\
        GCTGAGGCTTAATTAAACC
        CGACT
        """
    )

    result_a_3 = Dseq.from_representation(
        """\
                    TCAGC
        GAATTAATTTGGAGTCG
        """
    )

    for t in [result_a_1, result_a_2, result_a_3]:
        assert t in a.cut(NtBbvCI())

    # raise NotImplementedError("No code covering this behavior yet")


def test_nickase_with_watson_site_before_crick_site():
    from pydna.user import NtBbvCI
    from pydna.dseq import Dseq

    a = Dseq.from_representation(
        """\
        AACCTCAGCACTATCTTAGCTGAGGCTTAATTA
        TTGGAGTCGTGATAGAATCGACTCCGAATTAAT
        """
    )

    result_a_1 = Dseq("AACC", crick="", ovhg=0)

    result_a_2 = Dseq.from_representation(
        """\
            TCAGCACTATCTTAGCTGAGGCTTAATTA
        TTGGAGTCGTGATAGAATCGACT
        """
    )

    result_a_3 = Dseq(watson="", crick="TAATTAAGCC", ovhg=10)

    for t in [result_a_1, result_a_2, result_a_3]:
        assert t in a.cut(NtBbvCI())


def test_nickase_with_pac_before():
    from pydna.user import NtBbvCI
    from pydna.dseq import Dseq
    from Bio.Restriction import PacI

    a = Dseq.from_representation(
        """\
        GCTGAGGCTTAATTAAACCCAGC
        CGACTCCGAATTAATTTGGGTCG
        """
    )
    res_pac_1 = Dseq.from_representation(
        """\
        GCTGAGGCTTAAT
        CGACTCCGAAT
        """
    )

    res_pac_2 = Dseq.from_representation(
        """\
          TAAACCCAGC
        TAATTTGGGTCG
        """
    )

    # This should work regardless, comes from another module
    pac_cuts = a.cut(PacI)
    assert pac_cuts[0] == res_pac_1
    assert pac_cuts[1] == res_pac_2

    res_pac_1_nickase_1 = Dseq.from_representation(
        """\
        GCTGAGGCTTAAT
        CGACT
        """
    )

    nickase = NtBbvCI()
    assert res_pac_1_nickase_1 in res_pac_1.cut(nickase)

    res_pac_2_nickase_1 = Dseq.from_representation(
        """\
                CAGC
        TAATTTGGGTCG
        """
    )

    assert res_pac_2_nickase_1 in res_pac_2.cut(nickase)


def test_nickase_with_pac_together():
    from pydna.user import NtBbvCI
    from pydna.dseq import Dseq
    from Bio.Restriction import PacI

    a = Dseq.from_representation(
        """\
        GCTGAGGCTTAATTAAACCTCAGC
        CGACTCCGAATTAATTTGGAGTCG
        """
    )

    a.cut(NtBbvCI() + PacI())  # Not implemented yet, because NtBbvCI is not a RestrictionType


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
