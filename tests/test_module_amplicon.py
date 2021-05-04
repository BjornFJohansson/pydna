import pytest


def test_amplicon():

    from pydna.amplify import Anneal
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer

    template = Dseqrecord("AAAtacactcaccgtctatcattatctactatcgactgtatcatctgatagcacTTT")

    p1 = Primer("CCCtacactcaccgtctatcattatc")
    p2 = Primer("GGGgtgctatcagatgatacagtcg")

    ann = Anneal((p1, p2), template)

    prod = ann.products[0]

    assert repr(prod) == "Amplicon(57)"

    assert prod._repr_html_() == "Amplicon(57)"

    from unittest.mock import MagicMock

    pp = MagicMock()

    prod._repr_pretty_(pp, None)

    # assert pp.text.assert_called_with('Amplicon(57)')

    fig = """\
       5tacactcaccgtctatcattatc...cgactgtatcatctgatagcac3
                                  ||||||||||||||||||||||
                                 3gctgacatagtagactatcgtgGGG5
    5CCCtacactcaccgtctatcattatc3
        |||||||||||||||||||||||
       3atgtgagtggcagatagtaatag...gctgacatagtagactatcgtg5"""

    import textwrap

    assert prod.figure() == textwrap.dedent(fig)

    # assert prod.program() == prod.taq_program()

    # assert prod.pfu_sso7d_program() == prod.dbd_program()

    from pydna.amplicon import Amplicon

    from Bio.Seq import Seq

    from pydna.seqrecord import SeqRecord

    arg = SeqRecord(Seq("aaa"))

    x = Amplicon.from_SeqRecord(arg)


def test_amplicon_dbd():

    from pydna.amplify import Anneal
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    from textwrap import dedent

    template = Dseqrecord(
        "GCGTCCAGCGGCTGCCCGAGGCGCCAAGTGCCCGGGCCGAGCCCGCATCTGAGGCCGCCGCGGGC"
    )

    p1 = Primer("GCGTCCAGCGGCTGCCCGAGG")
    p2 = Primer("GCCCGCGGCGGCCTCAGATGCGG")

    ann = Anneal((p1, p2), template)

    prod = ann.products[0]

    assert repr(prod) == "Amplicon(65)"

    fig = r"""
              Pfu-Sso7d (rate 15s/kb)
              Two-step|    30 cycles |      |65bp
              98.0°C  |98.0C         |      |Tm formula: Pydna tmbresluc
              _____ __|_____         |      |SaltC 50mM
              00min30s|10s  \        |      |Primer1C 1.0µM
                      |      \ 72.0°C|72.0°C|Primer2C 1.0µM
                      |       \______|______|GC 81%
                      |       0min 0s|10min |4-12°C
              """[
        1:
    ]
    fig = r"""
            |98°C|98°C      |    |tmf:71.6
            |____|____      |    |tmr:75.3
            |30s |10s \ 72°C|72°C|15s/kb
            |    |     \____|____|GC 81%
            |    |      0: 0|5min|65bp
            """[
        1:
    ]

    fig = dedent(fig)
    assert str(prod.dbd_program()) == fig


def test_amplicon_dbd_low_gc():

    from pydna.amplify import Anneal
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    from textwrap import dedent

    template = Dseqrecord("AAAATATTTTTATACATAATACAATTGTATATTCTTAAATAAAAAATACGTCATC")

    p1 = Primer("AAAATATTTTTATACAT")
    p2 = Primer("GATGACGTATTTTTTAT")

    ann = Anneal((p1, p2), template)

    prod = ann.products[0]

    assert repr(prod) == "Amplicon(55)"

    fig = r"""
            Pfu-Sso7d (rate 15s/kb)                 |55bp
            Three-step|          30 cycles   |      |Tm formula: Pydna tmbresluc
            98.0°C    |98.0°C                |      |SaltC 50mM
            __________|_____          72.0°C |72.0°C|Primer1C 1.0µM
            00min30s  |10s  \ 39.0°C ________|______|Primer2C 1.0µM
                      |      \______/ 0min 0s|10min |GC 14%
                      |        10s           |      |4-12°C
            """[
        1:
    ]

    fig = r"""
              |98°C|98°C               |    |tmf:32.6
              |____|_____          72°C|72°C|tmr:39.6
              |30s |10s  \ 35.6°C _____|____|15s/kb
              |    |      \______/ 0: 0|5min|GC 14%
              |    |       10s         |    |55bp
              """[
        1:
    ]
    fig = dedent(fig)

    assert str(prod.dbd_program()) == fig


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
