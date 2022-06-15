#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_cas9():

    from pydna.dseq import Dseq

    s = Dseq("gattcatgcatgtagcttacgtagtct")

    RNA = "catgcatgtagcttacgtag"

    assert slice(0, 21, 1), slice(21, 27, 1) == s.cas9(RNA)


def test_initialization():

    import pytest
    from pydna.dseq import Dseq

    obj = Dseq(b"aaa")
    assert obj.crick == "ttt"
    assert obj.watson == "aaa"

    obj = Dseq("a", "t", 0)
    assert obj * 3 == Dseq("aaa", "ttt", 0)
    assert not obj == 123
    assert obj * 0 == Dseq("")

    with pytest.raises(TypeError):
        obj * 2.3

    assert obj.lseguid() == "bc1M4j2I4u6VaLpUbAB8Y9kTHBs"

    assert obj == Dseq("a", "t", circular=False, linear=True)

    with pytest.raises(ValueError):
        Dseq("a", ovhg=0)

    with pytest.raises(ValueError):
        Dseq("ttt", "tt")

    with pytest.raises(ValueError):
        Dseq("ttt", "aa")

    obj2 = Dseq("gata")

    assert obj2.linear == True
    assert obj2.circular == False

    l = Dseq("gt")
    c = l.looped()

    assert l.linear
    assert not l.circular
    assert c.circular
    assert not c.linear

    assert Dseq("gt", linear=None, circular=None) == l
    assert Dseq("gt", linear=None, circular=False) == l
    assert Dseq("gt", linear=None, circular=True) == c
    assert Dseq("gt", linear=False, circular=None) == c
    assert Dseq("gt", linear=False, circular=False) == l
    assert Dseq("gt", linear=False, circular=True) == c
    assert Dseq("gt", linear=True, circular=None) == l
    assert Dseq("gt", linear=True, circular=False) == l
    assert Dseq("gt", linear=True, circular=True) == l

    assert Dseq.from_string("A") == Dseq("A") == Dseq("A", linear=True)
    assert (
        Dseq.from_string("A", linear=False, circular=True)
        == Dseq("A", circular=True)
        == Dseq("A", linear=False)
    )


    obj3 = Dseq.from_representation("""
                                   GGATCC
                                   CCTAGG
                                   """)
    assert obj3.ovhg == 0

    obj3 = Dseq.from_representation("""
                                   aGGATCC
                                    CCTAGGg
                                   """)
    assert obj3.ovhg == -1

    obj3 = Dseq.from_representation("""
                                    GGATCCg
                                   aCCTAGG
                                   """)
    assert obj3.ovhg == 1


def test_cut_around_and_religate():
    from pydna.dseq import Dseq
    from pydna.utils import eq
    from Bio.Restriction import KpnI, BamHI, Acc65I

    def cut_and_religate_Dseq(seq_string, enz, top):
        ds = Dseq(seq_string, linear=top)
        frags = list(ds.cut(enz))
        if not frags:
            return
        a = frags.pop(0)
        for f in frags:
            a += f
        if not top:
            a = a.looped()
        assert eq(a, ds)

    seqs = [
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", BamHI, False),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", Acc65I, False),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", KpnI, False),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", Acc65I, True),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", KpnI, True),
        ("aaaGGTACCcccGGTACCCgggGGTACCttt", BamHI, True),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [Acc65I, BamHI], False),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [KpnI, BamHI], False),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI, Acc65I], False),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI, KpnI], False),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [Acc65I, BamHI], True),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [KpnI, BamHI], True),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI, Acc65I], True),
        ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI, KpnI], True),
    ]

    for s in seqs:
        sek, enz, lin = s
        for i in range(len(sek)):
            zek = sek[i:] + sek[:i]
            cut_and_religate_Dseq(zek, enz, lin)


def test_Dseq_cutting_adding():
    from pydna.dseq import Dseq
    from Bio.Restriction import BamHI, PstI, EcoRI

    a = Dseq(
        "GGATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGGATCC",
        "CCTAGGagtagatgatagtagcatcgcatgactagataagacgacgagtagtagccatgagagatattaatatatatatacgcgcaCCTAGG"[
            ::-1
        ],
        linear=True,
        ovhg=0,
    )

    b = a.cut(BamHI)[1]

    assert (
        b.watson
        == "GATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtG"
    )
    assert (
        b.crick
        == "GATCCacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaG"
    )
    c = Dseq(
        "nCTGCAGtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGAATTCn",
        "nGACGTCagtagatgatagtagcatcgcatgactagataagacgacgagtagtagccatgagagatattaatatatatatacgcgcaCTTAAGn"[
            ::-1
        ],
        linear=True,
        ovhg=0,
    )

    f, d, l = c.cut((EcoRI, PstI))

    assert (
        d.watson
        == "GtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtG"
    )
    assert (
        d.crick
        == "AATTCacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaCTGCA"
    )

    e = Dseq(
        "nGAATTCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtCTGCAGn",
        "nCTTAAGagtagatgatagtagcatcgcatgactagataagacgacgagtagtagccatgagagatattaatatatatatacgcgcaGACGTCn"[
            ::-1
        ],
        linear=True,
        ovhg=0,
    )

    f = e.cut((EcoRI, PstI))[1]

    assert (
        f.watson
        == "AATTCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtCTGCA"
    )
    assert (
        f.crick
        == "GacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaG"
    )


def test_dseq():

    import textwrap
    from pydna.dseq import Dseq

    obj1 = Dseq("a", "t", circular=True)
    obj2 = Dseq("a", "t")

    with pytest.raises(TypeError):
        obj1 + obj2

    with pytest.raises(TypeError):
        obj2 + obj1

    with pytest.raises(TypeError):
        obj1 + ""

    with pytest.raises(AttributeError):
        obj2 + ""

    obj1 = Dseq("at", "t")
    obj2 = Dseq("a", "t")

    with pytest.raises(TypeError):
        obj1 + obj2

    obj = Dseq("aaa", "ttt", circular=True)
    assert obj[1:2] == Dseq("a", "t", 0)

    assert obj[:] == Dseq("aaa", "ttt", circular=False)

    obj = Dseq("atg", "cat", 0, circular=False)

    assert obj[1:2]._data == b"atg"[1:2]

    assert obj[2:1]._data == b"atg"[2:1]

    assert obj.reverse_complement() == obj.rc() == Dseq("cat", "atg", 0)

    obj = Dseq("atg", "cat", circular=True)

    assert obj.looped() == obj

    assert obj[:] == Dseq("atg", "cat", 0, circular=False)

    assert obj[1:2]._data == b"atg"[1:2]

    assert obj[2:1]._data == b"ga"

    obj = Dseq("G", "", 0)
    assert obj.five_prime_end() == ("5'", "g")
    obj = Dseq("", "C", 0)
    assert obj.five_prime_end() == ("3'", "c")

    obj = Dseq("ccGGATCC", "aaggatcc", -2)
    assert obj._data == b"ccGGATCCtt"
    assert str(obj.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-10)
    ccGGATCC
      cctaggaa
    """
    ).strip()
    assert repr(obj) == rpr

    assert obj[3] == Dseq("G", "c", 0)

    assert obj.fill_in() == Dseq("ccGGATCCtt", "aaggatccgg", 0)

    assert obj + Dseq("") == obj
    assert Dseq("") + obj == obj

    obj = Dseq("gatcAAAAAA", "gatcTTTTTT")
    assert obj.fill_in("gatc") == Dseq("gatcAAAAAAgatc", "gatcTTTTTTgatc")
    assert obj.fill_in("atc") == obj
    assert obj.fill_in("ac") == obj
    assert obj.fill_in("at") == obj

    obj = Dseq("AAAAAAgatc", "TTTTTTgatc")
    assert obj.fill_in("gatc") == obj
    assert obj.fill_in("atc") == obj
    assert obj.fill_in("ac") == obj
    assert obj.fill_in("at") == obj

    obj = Dseq("gatcAAAAAA", "gatcTTTTTT")
    assert obj.t4() == Dseq("gatcAAAAAAgatc", "gatcTTTTTTgatc")

    assert obj.t4("at") == obj
    assert obj.t4("atg") == Dseq("gatcAAAAAAgat", "gatcTTTTTTgat")
    assert obj.t4("atgc") == Dseq("gatcAAAAAAgatc", "gatcTTTTTTgatc")
    assert obj.mung() == Dseq("AAAAAA", "TTTTTT")

    obj = Dseq("AAAAAAgatc", "TTTTTTgatc")
    assert obj.t4() == obj.t4("at") == Dseq("AAAAAA")
    assert obj.t4("atc") == obj.t4("atg") == obj.t4("atcg") == Dseq("AAAAAA")

    assert Dseq("GGATCC", "GGATCC").t4() == Dseq("GGATCC", "GGATCC")
    assert Dseq("GGATCCa", "GGATCC").t4() == Dseq("GGATCC", "GGATCC")
    assert Dseq("aGGATCC", "GGATCC").t4() == Dseq("aGGATCC", "GGATCCt")
    assert Dseq("aGGATCCa", "GGATCC").t4() == Dseq("aGGATCC", "GGATCCt")
    assert Dseq("GGATCC", "aGGATCC").t4() == Dseq("GGATCCt", "aGGATCC")
    assert Dseq("GGATCC", "GGATCCa").t4() == Dseq("GGATCC", "GGATCC")
    assert Dseq("GGATCC", "aGGATCCa").t4() == Dseq("GGATCCt", "aGGATCC")

    assert Dseq("GGATCC", "ATCC").t4("g") == Dseq("gg", "", ovhg=0)
    assert Dseq("GGATCC", "GGATCC").t4("gat") == Dseq("ggat", "ggat", ovhg=-2)

    a2 = Dseq("ccGGATCCaa", "ggatcc", -2)
    assert a2._data == b"ccGGATCCaa"
    assert a2._data == b"ccGGATCCaa"
    assert str(a2.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-10)
    ccGGATCCaa
      cctagg
    """
    ).strip()
    assert repr(a2) == rpr

    a3 = Dseq("ccGGATCC", "ggatcc", -2)
    assert a3._data == b"ccGGATCC"
    assert a3._data == b"ccGGATCC"
    assert str(a3.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-8)
    ccGGATCC
      cctagg
    """
    ).strip()
    assert repr(a3) == rpr

    b = Dseq("GGATCC", "aaggatcccc", 2)
    assert b._data == b"ggGGATCCtt"
    assert b._data == b"ggGGATCCtt"
    assert str(b.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-10)
      GGATCC
    cccctaggaa
    """
    ).strip()
    assert repr(b) == rpr

    b2 = Dseq("GGATCCaa", "ggatcccc", 2)
    assert b2._data == b"ggGGATCCaa"
    assert b2._data == b"ggGGATCCaa"
    assert str(b2.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-10)
      GGATCCaa
    cccctagg
    """
    ).strip()
    assert repr(b2) == rpr

    assert b2.lseguid() == "yz2Dl2jN1FxbE1BNV1hvT2r4S3Y"
    assert b2.rc().lseguid() == "yz2Dl2jN1FxbE1BNV1hvT2r4S3Y"

    b3 = Dseq("GGATCC", "ggatcccc", 2)
    assert b3._data == b"ggGGATCC"
    assert b3._data == b"ggGGATCC"
    assert str(b3.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-8)
      GGATCC
    cccctagg
    """
    ).strip()
    assert repr(b3) == rpr

    c = Dseq("GGATCCaaa", "ggatcc", 0)
    assert c._data == b"GGATCCaaa"
    assert c._data == b"GGATCCaaa"
    assert str(c.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-9)
    GGATCCaaa
    cctagg
    """
    ).strip()
    assert repr(c) == rpr

    d = Dseq("GGATCC", "aaaggatcc", 0)
    assert d._data == b"GGATCCttt"
    assert d._data == b"GGATCCttt"
    assert str(d.mung()) == "GGATCC"
    rpr = textwrap.dedent(
        """
    Dseq(-9)
    GGATCC
    cctaggaaa
    """
    ).strip()
    assert repr(d) == rpr

    obj = Dseq("GGATCCaaa", "ggatcc", 0)
    from Bio.Restriction import BamHI

    frag1 = Dseq("G", "gatcc", 0)
    frag2 = Dseq("GATCCaaa", "g", -4)

    assert obj.cut(BamHI) == (frag1, frag2)

    assert frag1 + frag2 == obj

    assert obj.lseguid() == "XtNVVlmD0oknESP6ErbhOUlPp9M"

    assert frag1.lseguid() == "uTTkMqUMJ524xcTzAKlAPsq6q7g"
    assert frag2.lseguid() == "pXjshRKBvk2SnAZeuVgsMWXmkfE"

    assert frag1.rc().lseguid() == "uTTkMqUMJ524xcTzAKlAPsq6q7g"
    assert frag2.rc().lseguid() == "pXjshRKBvk2SnAZeuVgsMWXmkfE"

    obj = Dseq("tagcgtagctgtagtatgtgatctggtcta", "tagaccagatcacatactacagctacgcta")
    assert (
        repr(obj)
        == "Dseq(-30)\ntagcgtagctgtagtatgtgatctggtcta\natcgcatcgacatcatacactagaccagat"
    )

    obj2 = Dseq("tagcgtagctgtagtatgtgatctggtcta")

    obj3 = obj = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta", "tagaccagatcacatactacagctacgcta", 0
    )

    assert obj == obj2 == obj3

    assert obj.find("ggatcc") == -1

    assert obj.find("tgtagta") == 9

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaa", "ttagaccagatcacatactacagctacgcta")

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaa", "CCCttagaccagatcacatactacagctacgcta")

    assert repr(obj) == "Dseq(-34)\ntagc..ctaa   \natcg..gattCCC"

    obj = Dseq("tagcgtagctgtagtatgtgatctggtctaaCCC", "ttagaccagatcacatactacagctacgcta")

    assert repr(obj) == "Dseq(-34)\ntagc..ctaaCCC\natcg..gatt   "

    obj = Dseq("agcgtagctgtagtatgtgatctggtctaa", "ttagaccagatcacatactacagctacgcta")
    assert repr(obj) == "Dseq(-31)\n agcg..ctaa\natcgc..gatt"

    obj = Dseq("Atagcgtagctgtagtatgtgatctggtctaa", "ttagaccagatcacatactacagctacgcta")
    assert repr(obj) == "Dseq(-32)\nAtagc..ctaa\n atcg..gatt"

    obj = Dseq(
        "tagcgtagctgtagtatgtgatctggtctaa", "tatcgcatcgacatcatacactagaccagatt"[::-1]
    )

    assert repr(obj) == "Dseq(-32)\n tagc..ctaa\ntatcg..gatt"

    assert round(obj.mw(), 1) == 19535.6

    obj1 = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta",
        "tagaccagatcacatactacagctacgcta",
        circular=True,
        linear=False,
    )
    obj2 = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta",
        "tagaccagatcacatactacagctacgcta",
        circular=True,
    )
    obj3 = Dseq(
        "tagcgtagctgtagtatgtgatctggtcta", "tagaccagatcacatactacagctacgcta", linear=False
    )

    assert obj1 == obj2 == obj3

    assert obj1.find("ggatcc") == -1

    assert obj1.find("tgtagta") == 9

    assert (
        Dseq(
            "tagcgtagctgtagtatgtgatctggtcta", "tagaccagatcacatactacagctacgcta"
        ).looped()
        == obj1
    )

    from Bio.Restriction import BglII, BamHI

    obj = Dseq("ggatcc")

    assert BglII in obj.no_cutters()
    assert BamHI not in obj.no_cutters()

    assert BamHI in obj.unique_cutters()

    assert BamHI in obj.once_cutters()

    assert BamHI in (obj + obj).twice_cutters()
    assert BamHI not in obj.twice_cutters()

    assert BamHI in obj.n_cutters(1)
    assert BamHI in obj.cutters()

    from Bio.Restriction import RestrictionBatch

    rb = RestrictionBatch((BamHI, BglII))

    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("ggatccAGATCT")

    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("AGATCTggatcc")

    assert obj.cut(rb) == obj.cut(BamHI, BglII) == obj.cut(BglII, BamHI)

    obj = Dseq("ggatccAGATCT", circular=True)

    assert obj.cut(rb) == obj.cut(BamHI, BglII) != obj.cut(BglII, BamHI)

    obj = Dseq("AGATCTggatcc", circular=True)

    assert obj.cut(rb) == obj.cut(BglII, BamHI) != obj.cut(BamHI, BglII)


def test_Dseq_slicing():
    from pydna.dseq import Dseq
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec
    from Bio.Restriction import BamHI

    a = Dseq("ggatcc", "ggatcc", 0)

    assert a[:].watson == a.watson
    assert a[:].crick == a.crick
    assert a.ovhg == a[:].ovhg
    b, c = a.cut(BamHI)
    d = b[1:5]
    e = d.rc()
    # assert  d+e == Dseq("gatc","gatc",0)
    assert e + d == Dseq("gatc", "gatc", 0)


def test_Dseq_slicing2():
    from pydna.dseq import Dseq
    from Bio.Restriction import BamHI, EcoRI, KpnI

    a = Dseq("aaGGATCCnnnnnnnnnGAATTCccc", circular=True)

    assert a.cut(EcoRI, BamHI, KpnI,) == a.cut(
        BamHI,
        EcoRI,
        KpnI,
    )[::-1]


def test_Dseq___getitem__():
    from pydna.dseq import Dseq

    s = Dseq("GGATCC", circular=False)
    assert s[1:-1] == Dseq("GATC", circular=False)
    t = Dseq("GGATCC", circular=True)
    assert t[1:5] == Dseq("GATC")
    assert t[1:5].__dict__ == Dseq("GATC").__dict__
    assert s[1:5] == Dseq("GATC")
    assert s[1:5] == Dseq("GATC", circular=False)
    assert s[5:1:-1] == Dseq("CCTA")

    assert t[5:1] == Dseq("CG")

    assert s[9:1] == Dseq("")
    assert t[9:1] == Dseq("")


def test_cut_circular():

    from pydna.dseq import Dseq
    from Bio.Restriction import BsaI, KpnI, Acc65I, NotI

    test = "aaaaaaGGTACCggtctcaaaa"

    for i in range(len(test)):

        nt = test[i:] + test[:i]

        a = Dseq(nt, circular=True).cut(Acc65I)[0]  # G^GTACC

        assert a.watson.upper() == "GTACCGGTCTCAAAAAAAAAAG"
        assert a.crick.upper() == "GTACCTTTTTTTTTTGAGACCG"
        assert a.ovhg == -4  # CggtctcaaaaaaaaaaGGTAC
        b = Dseq(nt, circular=True).cut(KpnI)[0]  # GGTAC^C
        assert b.watson.upper() == "CGGTCTCAAAAAAAAAAGGTAC"
        assert b.crick.upper() == "CTTTTTTTTTTGAGACCGGTAC"
        assert b.ovhg == 4
        c = Dseq(nt, circular=True).cut(BsaI)[0]  # ggtctcnnn
        assert c.watson.upper() == "AAAAAAAAAGGTACCGGTCTCA"
        assert c.crick.upper() == "TTTTTGAGACCGGTACCTTTTT"
        assert c.ovhg == -4
        d = Dseq(nt, circular=True).cut(NotI)
        assert d == ()


def test_repr():
    from pydna.dseq import Dseq

    a = Dseq("gattcgtatgctgatcgtacgtactgaaaac")

    assert repr(a) == "Dseq(-31)\ngatt..aaac\nctaa..tttg"

    b = Dseq("gattcgtatgctgatcgtacgtactgaaaac", "gactagcatgcatgacttttc"[::-1])

    assert repr(b) == "Dseq(-31)\ngattcgtatgctga..aaac\n          gact..tttc"

    c = Dseq("gattcgtatgctgatcgtacgtactgaaaac", "actagcatgcatgacttttc"[::-1])

    assert repr(c) == "Dseq(-31)\ngatt..atgctgat..aaac\n          acta..tttc"

    d = Dseq("gattcgtatgctgatcgtacg", "gactagcatgc"[::-1])

    assert repr(d) == "Dseq(-21)\ngattcgtatgctgatcgtacg\n          gactagcatgc"

    e = Dseq("gactagcatgcatgacttttc", "gattcgtatgctgatcgtacgtactgaaaac"[::-1])

    assert repr(e) == "Dseq(-31)\n          gact..tttc\ngattcgtatgctga..aaac"

    f = Dseq("Ggactagcatgcatgacttttc", "gattcgtatgctgatcgtacgtactgaaaac"[::-1])

    assert repr(f) == "Dseq(-31)\n         Ggac..tttc\ngattcgtatgctg..aaac"

    g = Dseq("gattcgtatgctgatcgtacgtactgaaaac", "ctaagcatacgactagc"[::-1])

    assert repr(g) == "Dseq(-31)\ngatt..atcgtacg..aaac\nctaa..tagc          "

    h = Dseq("cgtatgctgatcgtacgtactgaaaac", "gcatacgactagc"[::-1])

    assert repr(h) == "Dseq(-27)\ncgtatgctgatcgtacgtactgaaaac\ngcatacgactagc"

    i = Dseq("cgtatgctgatcgtacgtactgaaaacagact", "gcatacgactagc"[::-1])

    assert repr(i) == "Dseq(-32)\ncgta..atcgtacg..gact\ngcat..tagc          "

    j = Dseq("gattcgtatgctgatcgtacgtactgaaaac", "acAAGGAGAGAtg", ovhg=11)

    assert repr(j) == "Dseq(-42)\n          gattcg..aaac\ngtAG..GGAAca          "

    k = Dseq("g", "gattcgtatgctgatcgtacgtactgaaaac", ovhg=0)

    assert repr(k) == "Dseq(-31)\ng          \ncaaaa..ttag"

    x = Dseq("gattcgtatgctgatcgtacgtactgaaaa")

    assert (
        repr(x)
        == "Dseq(-30)\ngattcgtatgctgatcgtacgtactgaaaa\nctaagcatacgactagcatgcatgactttt"
    )

    y = Dseq("gattcgtatgctgatcgtacgtactgaaaa", "gactagcatgcatgactttt"[::-1])

    assert (
        repr(y)
        == "Dseq(-30)\ngattcgtatgctgatcgtacgtactgaaaa\n          gactagcatgcatgactttt"
    )

    z = Dseq("gattcgtatgctgatcgtacgtactgaaaa", "actagcatgcatgactttt"[::-1])

    assert (
        repr(z)
        == "Dseq(-30)\ngattcgtatgctgatcgtacgtactgaaaa\n           actagcatgcatgactttt"
    )


def test_shifted():
    from pydna.dseq import Dseq

    a = Dseq("gatc", circular=True)

    assert a.shifted(1) == Dseq("atcg", circular=True)

    assert a.shifted(4) == a

    b = Dseq("gatc", circular=False)
    with pytest.raises(TypeError):
        b.shifted(1)


def test_misc():

    from pydna.dseq import Dseq

    x = Dseq("ctcgGCGGCCGCcagcggccg", circular=True)

    from Bio.Restriction import NotI

    a, b = x.cut(NotI)

    z = (a + b).looped()

    assert z.shifted(5) == x


def test_cut_missing_enzyme():

    from pydna.dseq import Dseq

    x = Dseq("ctcgGCGGCCGCcagcggccg")

    from Bio.Restriction import PstI

    assert x.cut(PstI) == ()

    x = Dseq("ctcgGCGGCCGCcagcggccg", circular=True)

    assert x.cut(PstI) == ()


def test_cut_with_no_enzymes():

    from pydna.dseq import Dseq

    x = Dseq("ctcgGCGGCCGCcagcggccg")

    assert x.cut([]) == ()

    x = Dseq("ctcgGCGGCCGCcagcggccg", circular=True)

    assert x.cut([]) == ()


def test_transcribe():

    from pydna.dseq import Dseq

    x = Dseq("ATGAAATAA")

    assert str(x.transcribe()) == 'AUGAAAUAA'

    assert str(x.reverse_complement().transcribe()) == 'UUAUUUCAU'


def test_translate():

    from pydna.dseq import Dseq

    x = Dseq("ATGAAATAA")

    assert str(x.translate()) == 'MK*'

    assert str(x.reverse_complement().translate()) == 'LFH'


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
