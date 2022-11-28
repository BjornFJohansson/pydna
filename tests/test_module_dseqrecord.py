#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from pydna import _PydnaWarning

def test_orfs():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("atgaaattttaa")

    assert s.orfs(2) == (s,)

def test_cas9():

    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("gattcatgcatgtagcttacgtagtct")

    RNA = "catgcatgtagcttacgtag"

    (f1, f2), (f3,) = s.cas9(RNA)

    assert f1.seq == Dseqrecord("gattcatgcatgtagcttacg").seq
    assert f2.seq == Dseqrecord("tagtct").seq
    assert f3.seq == Dseqrecord("gattcatgcatgtagcttacgtagtct").seq


def test_FadiBakoura():
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord

    dseq_record = Dseqrecord(Dseq("ACTCTTTCTCTCTCT", circular=True))
    dseq_record.features = [SeqFeature(FeatureLocation(start=2, end=4))]
    # below assert fails (not anymore)
    assert len(dseq_record[6:1].features) == 0


def test_IPython_missing(monkeypatch):
    import IPython
    from pydna import dseqrecord

    # assert dseqrecord._display is IPython.display.display
    assert dseqrecord._display_html == IPython.display.display_html
    import sys
    import copy

    fakesysmodules = copy.copy(sys.modules)
    fakesysmodules["IPython.display"] = None
    monkeypatch.delitem(sys.modules, "IPython.display")
    monkeypatch.setattr("sys.modules", fakesysmodules)
    from importlib import reload

    reload(dseqrecord)
    from pydna import dseqrecord

    assert dseqrecord._display_html("item") == "item"
    #assert dseqrecord._HTML("item") == "item"


def test_initialization():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    a = []

    a.append(Dseqrecord("attt"))
    a.append(Dseqrecord(Dseq("attt")))
    a.append(Dseqrecord(Seq("attt")))
    a.append(Dseqrecord(Srec(Seq("attt"))))
    a.append(Dseqrecord(Dseqrecord("attt")))
    a.append(Dseqrecord(Dseqrecord(Dseq("attt", circular=True)), circular=False))

    for b in a:
        assert type(b.seq) == Dseq
        assert str(b.seq.watson) == "attt"
        assert str(b.seq.crick) == "aaat"
        assert str(b.seq) == "attt"
        assert str(b.seq) == "attt"
        assert b.linear == b.seq.linear
        assert b.linear == True
        assert b.circular == False
        assert b.seq.linear == True
        assert b.seq.circular == False

    a = []
    a.append(Dseqrecord("attt", circular=True))
    a.append(Dseqrecord(Dseq("attt"), circular=True))
    a.append(Dseqrecord(Seq("attt"), circular=True))
    a.append(Dseqrecord(Srec(Seq("attt")), circular=True))
    a.append(Dseqrecord(Dseqrecord("attt"), circular=True))

    for b in a:
        assert type(b.seq) == Dseq
        assert str(b.seq.watson) == "attt"
        assert str(b.seq.crick) == "aaat"
        assert str(b.seq) == "attt"
        assert str(b.seq) == "attt"
        assert b.linear == b.seq.linear
        assert b.linear == False
        assert b.circular == True
        assert b.seq.linear == False
        assert b.seq.circular == True

    a = []
    a.append(Dseqrecord(Dseq("attt", circular=True), circular=True))
    a.append(Dseqrecord(Dseq("attt", circular=False), circular=True))
    a.append(Dseqrecord(Dseq("attt", circular=True), circular=False))
    a.append(Dseqrecord(Dseq("attt", circular=False), circular=False))

    circular = [True, True, False, False]
    linear = [False, False, True, True]

    for b, ci, li in zip(a, circular, linear):
        assert type(b.seq) == Dseq
        assert str(b.seq.watson) == "attt"
        assert str(b.seq.crick) == "aaat"
        assert str(b.seq) == "attt"
        assert str(b.seq) == "attt"
        assert b.linear == b.seq.linear
        assert b.linear == li
        assert b.circular == ci
        assert b.seq.linear == li
        assert b.seq.circular == ci

    a = []
    ds = Dseq("attt", "taaa")
    assert ds.linear == True
    assert ds.ovhg == -1
    assert str(ds.watson) == "attt"
    assert str(ds.crick) == "taaa"

    #   attt
    #    aaat

    a.append(Dseqrecord(ds, circular=False))
    assert ds.linear == True
    a.append(Dseqrecord(ds, linear=True))
    assert ds.linear == True

    a.append(Dseqrecord(ds, circular=True))
    assert ds.linear == True
    a.append(Dseqrecord(ds, linear=False))
    assert ds.linear == True

    circular = [False, False, True, True]
    linear = [True, True, False, False]
    crick = ["taaa", "taaa", "aaat", "aaat"]
    sek = ["attta", "attta", "attt", "attt"]
    for b, ci, li, s, cri in zip(a, circular, linear, sek, crick):

        assert type(b.seq) == Dseq
        assert str(b.seq.watson) == "attt"
        assert str(b.seq.crick) == cri
        assert str(b.seq) == s

        assert b.linear == b.seq.linear
        assert b.linear == li
        assert b.circular == ci
        assert b.seq.linear == li
        assert b.seq.circular == ci

    a = []
    ds = Dseq("attt", "caaa")
    assert ds.linear == True
    assert ds.ovhg == -1

    a.append(Dseqrecord(ds, circular=False))
    assert ds.linear == True
    a.append(Dseqrecord(ds, linear=True))
    assert ds.linear == True

    with pytest.raises(TypeError):
        Dseqrecord(ds, circular=True)

    assert ds.linear == True

    with pytest.raises(TypeError):
        Dseqrecord(ds, linear=False)

    assert ds.linear == True

    with pytest.raises(ValueError):
        b = Dseqrecord([])

    with pytest.raises(ValueError):
        b = Dseqrecord(("a",))

    with pytest.raises(ValueError):
        b = Dseqrecord(0)

    input = """
            LOCUS       New_DNA                    4 bp ds-DNA     linear       30-MAR-2013
            DEFINITION  .
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            COMMENT
            COMMENT     ApEinfo:methylated:1
            FEATURES             Location/Qualifiers
                 misc_feature    2..3
                                 /label=NewFeature
                                 /ApEinfo_fwdcolor=cyan
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
            ORIGIN
                    1 acgt
            //
            """
    a = read(input)

    assert a.features[0].extract(a).seq.watson == "CG"

    b = a + a

    for f in b.features:
        assert b.features[0].extract(a).seq.watson == "CG"
    feature = a.features[0]

    s = Dseq("agctt", "agcta")
    # print s.fig()
    # Dseq(-6)
    # agctt
    # atcga
    b = Dseqrecord(s)
    b.features.append(feature)
    cb = Dseqrecord(b, circular=True)
    assert (
        b.features[0].extract(b).seq.watson.lower()
        == cb.features[0].extract(b).seq.watson.lower()
    )
    assert (
        b.features[0].extract(b).seq.crick.lower()
        == cb.features[0].extract(b).seq.crick.lower()
    )
    s = Dseq("aagct", "aagct")
    # print s.fig()
    # Dseq(-6)
    # aagct
    # tcgaa
    b = Dseqrecord(s)
    with pytest.raises(TypeError):
        cb = Dseqrecord(b, circular=True)

    s = Dseq("agctt", "agcta")
    # print s.fig()
    # Dseq(-6)
    # agcta
    # ttcga

    b = Dseqrecord(s)
    b.features.append(feature)
    cb = Dseqrecord(b, circular=True)
    assert (
        b.features[0].extract(b).seq.watson.lower()
        == cb.features[0].extract(b).seq.watson.lower()
    )
    assert (
        b.features[0].extract(b).seq.crick.lower()
        == cb.features[0].extract(b).seq.crick.lower()
    )


def test_linear_circular():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    """ test Dseqrecord linear & circular property"""
    a = Dseqrecord("attt")
    a.stamp("useguid")
    assert a.stamp("useguid")
    a = Dseqrecord("attt", linear=True)
    assert a.linear == True
    assert a.circular == False
    assert a.rc().linear == True
    assert a.rc().circular == False
    assert a.seq.linear == True
    assert a.seq.circular == False

    a = Dseqrecord("attt", linear=False)
    assert a.linear == False
    assert a.circular == True
    assert a.rc().linear == False
    assert a.rc().circular == True
    assert a.seq.linear == False
    assert a.seq.circular == True

    a = Dseqrecord("attt", circular=True)
    assert a.linear == False
    assert a.circular == True
    assert a.rc().linear == False
    assert a.rc().circular == True
    assert a.seq.linear == False
    assert a.seq.circular == True

    a = Dseqrecord("attt", circular=False)
    assert a.linear == True
    assert a.circular == False
    assert a.rc().linear == True
    assert a.rc().circular == False
    assert a.seq.linear == True
    assert a.seq.circular == False

import pydna
def test_stamp():
    from pydna.dseqrecord import Dseqrecord


    lin = Dseqrecord("attt")
    assert lin.stamp("useguid") == "ot6JPLeAeMmfztW1736Kc6DAqlo"
    assert lin.stamp("useguid")[:35] == "uSEGUID ot6JPLeAeMmfztW1736Kc6DAqlo"

    crc = Dseqrecord("attt", circular=True)
    assert crc.stamp("cseguid") == "oopV-6158nHJqedi8lsshIfcqYA"
    crc.annotations["comment"] = "cSEGUID_oopV-6158nHJqedi8lsshIfcqYZ"
    with pytest.warns(_PydnaWarning):
        crc.stamp("cseguid")

    from pydna.dseq import Dseq

    blunt = Dseqrecord(Dseq("aa"))

    assert blunt.stamp("useguid") == "gBw0Jp907Tg_yX3jNgS4qQWttjU"

    staggered = Dseqrecord(Dseq("aa", "tta"))
    assert staggered.stamp("lseguid") == "qlQKjCD04EhPHtsdnEQ0DriOr10"

    staggered = Dseqrecord(Dseq("aa", "att"))
    assert staggered.stamp("lseguid") == "mvhPR1scFLQfgDn3olxegJmckgg"

    staggered = Dseqrecord(Dseq("aa", "atta"))
    assert staggered.stamp("lseguid") == "hPVQIzC3k29tDuAu_C5xEgLsFZs"


def test_revcomp():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    #      ----
    #     attcccgggg
    #     taagggcccc
    #
    #      ttcc
    #
    #
    #     ccccgggaat
    #     ggggccctta
    #          ----

    a = Dseqrecord("attcccgggg")

    a.add_feature(1, 5)

    assert str(a.features[0].extract(a).seq) == "ttcc"

    assert a.features[0].location.strand == 1

    rc = a.rc()

    assert str(rc.features[0].extract(rc).seq) == "ttcc"

    assert rc.features[0].location.strand == -1

    a = Dseqrecord("attcccgggg")

    a.add_feature(1, 5)

    a.features[0].location.strand = None

    rc = a.rc()

    assert rc.features[0].strand == None


def test_m():
    from pydna.dseqrecord import Dseqrecord
    s = Dseqrecord("A" * 5000)
    assert f'{s.m():.3e}' == '1.544e-07'


def test_extract_feature():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    s.add_feature(3, 9)
    s.extract_feature(0).seq == Dseqrecord("GGATCC").seq


def test_find():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    assert s.find(Dseqrecord("ggatcc")) == 3
    assert s.find(Dseqrecord("ggatTc")) == -1
    s = Dseqrecord("tttGGATCCaaa", circular=True)
    assert s.find(Dseqrecord("aaattt")) == 9


def test_find_aa():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    assert s.find_aa("FGSK") == slice(0, 12, None)
    assert s.find_aa("MMMM") == None


def test_str():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    s.annotations = {"date": "03-JAN-2018"}
    assert (
        str(s)
        == "Dseqrecord\ncircular: False\nsize: 12\nID: id\nName: name\nDescription: description\nNumber of features: 0\n/date=03-JAN-2018\nDseq(-12)\ntttGGATCCaaa\naaaCCTAGGttt"
    )
    s = s.looped()
    assert (
        str(s)
        == "Dseqrecord\ncircular: True\nsize: 12\nID: id\nName: name\nDescription: description\nNumber of features: 0\n/date=03-JAN-2018\nDseq(o12)\ntttGGATCCaaa\naaaCCTAGGttt"
    )


def test___contains__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    assert "ggatcc" in s


def test_cseguid():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa")
    with pytest.raises(TypeError):
        s.cseguid()
    s = Dseqrecord("tttGGATCCaaa", circular=True)
    assert s.cseguid() == "54zt6fixMmNwFTsUVILdrdIRi2U"


def test_lseguid():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("tttGGATCCaaa", circular=True)
    with pytest.raises(TypeError):
        s.lseguid()
    s = Dseqrecord("tttGGATCCaaa")
    assert s.lseguid() == "qr7ePrqTrFyQt0fsTsiltN20oEk"


def test_format():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=True)
    s.format("gb")
    s.format("genbank")
    s = Dseqrecord("GGATCC", circular=False)
    s.format("gb")
    s.format("genbank")

    s = Dseqrecord("GGATCC", circular=True)
    s.format("fasta")
    s = Dseqrecord("GGATCC", circular=False)
    s.format("fasta")


def test_write():
    from unittest.mock import patch
    from unittest.mock import mock_open
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=True)
    m = mock_open()
    with patch("pydna.dseqrecord.open", m):
        s.write(filename="AAA.gb")
    m.mock_calls
    m.assert_called_once_with("AAA.gb", "w", encoding='utf8')
    handle = m()
    handle.write.assert_called_once_with(Dseqrecord("GGATCC", circular=True).format())

    m = mock_open()
    with patch("pydna.dseqrecord.open", m):
        s.write()
    m.mock_calls
    m.assert_called_once_with("name.gb", "w", encoding='utf8')
    handle = m()
    handle.write.assert_called_once_with(Dseqrecord("GGATCC", circular=True).format())

    with pytest.raises(ValueError):
        s.write(filename=123)


def test_write_same_seq_to_existing_file(monkeypatch):

    import builtins
    from unittest.mock import patch
    from unittest.mock import mock_open
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read

    s = Dseqrecord("Ggatcc", circular=True)

    monkeypatch.setattr("pydna.dseqrecord._os.path.isfile", lambda x: True)
    m = mock_open(read_data=s.format())

    with patch("builtins.open", m) as d:
        s.write(filename="Ggatcc.gb")


def test_write_different_file_to_existing_file(monkeypatch):

    import builtins
    from unittest.mock import patch
    from unittest.mock import mock_open
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read

    s = Dseqrecord("Ggatcc", circular=True)
    d = Dseqrecord("Ggatcn", circular=True)

    monkeypatch.setattr("pydna.dseqrecord._os.path.isfile", lambda x: True)
    monkeypatch.setattr("pydna.dseqrecord._os.rename", lambda x, y: True)
    m = mock_open(read_data=d.format())

    with patch("builtins.open", m) as d:
        s.write(filename="Ggatcc.gb")


def test_write_different_file_to_stamped_existing_file(monkeypatch):

    import builtins
    from unittest.mock import patch
    from unittest.mock import mock_open
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read

    new = Dseqrecord("Ggatcc", circular=True)
    new.stamp("cSEGUID")
    old = Dseqrecord("Ggatcc", circular=True)
    old.stamp("cSEGUID")

    assert new.description[:35] == old.description[:35]

    monkeypatch.setattr("pydna.dseqrecord._os.path.isfile", lambda x: True)
    monkeypatch.setattr("pydna.dseqrecord._os.rename", lambda x, y: True)
    m = mock_open(read_data=old.format())

    with patch("builtins.open", m) as d:
        new.write(filename="Ggatcc.gb")

    new = Dseqrecord("Ggatcc", circular=True)

    with patch("builtins.open", m) as d:
        new.write(filename="Ggatcc.gb")

    new.description = "cSEGUID_NNNNNNNNNNNNNNNNNNNNNNNNNNN_2018-06-01T05:05:51.778066"

    with patch("builtins.open", m) as d:
        new.write(filename="Ggatcc.gb")

    new.description = "cSEGUID_N"

    m = mock_open(read_data=old.format())
    with patch("builtins.open", m) as d:
        new.write(filename="Ggatcc.gb")

    assert m.called
    # m.write().assert_called_once_with(new.format())
    assert m.call_count == 6
    assert m.mock_calls[0]
    assert m.mock_calls[4]


def test_write_different_file_to_stamped_existing_file2(monkeypatch):

    import builtins
    from unittest.mock import patch
    from unittest.mock import mock_open
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read

    new = Dseqrecord("Ggatcc", circular=True)
    new.stamp("cSEGUID")
    old = Dseqrecord("Ggatcc", circular=True)
    old.stamp("cSEGUID")

    assert new.description[:35] == old.description[:35]

    monkeypatch.setattr("pydna.dseqrecord._os.path.isfile", lambda x: True)
    monkeypatch.setattr("pydna.dseqrecord._os.rename", lambda x, y: True)
    m = mock_open(read_data=old.format())

    with patch("builtins.open", m) as d:
        new.write(filename="Ggatcc.gb")

    new = Dseqrecord("Ggatcc", circular=True)

    with patch("builtins.open", m) as d:
        new.write(filename="Ggatcc.gb")

    new.description = "cSEGUID_NNNNNNNNNNNNNNNNNNNNNNNNNNN_2018-06-01T05:05:51.778066"

    with patch("builtins.open", m) as d:
        new.write(filename="Ggatcc.gb")

    new.description = "cSEGUID_N"

    m = mock_open(read_data=old.format())
    with patch("builtins.open", m) as d:
        new.write(filename="Ggatcc.gb")

    assert m.called
    # m.write().assert_called_once_with(new.format())
    assert m.call_count == 6
    assert m.mock_calls[0]
    assert m.mock_calls[4]


# [call('Ggatcc.gb', 'r', encoding='utf-8'),
# call().__enter__(),
# call().read(),
# call().__exit__(None, None, None),
# call('Ggatcc.gb', 'w'),
# call().__enter__(),
# call().write('LOCUS       name                       6 bp    DNA     circular UNK 01-JUN-2018\nDEFINITION  cSEGUID_N\n            cSEGUID_6WVYnCK97MOPMOlbLHvMnd4XIEY_2018-06-01T06:05:08.951398.\nACCESSION   id\nVERSION     id\nKEYWORDS    .\nSOURCE      .\n  ORGANISM  .\n            .\nFEATURES             Location/Qualifiers\nORIGIN\n        1 ggatcc\n//'),
# call().__exit__(None, None, None)]


# monkeypatch.setitem(readers.__builtins__, 'open', open)
# def test_write_to_existing_file():
#    from unittest.mock import patch
#    from unittest.mock import mock_open
#    from pydna.dseqrecord  import Dseqrecord
#    from pydna.readers     import read
#
#    s = Dseqrecord("GGATCC", circular=True)
#    m = mock_open()
#    m.read_data=Dseqrecord("GGATCC", circular=True).format()
#    with patch('pydna.dseqrecord.open', m):
#        s.write(filename="AAA.gb")
#    m.mock_calls
#    m.assert_called_once_with("AAA.gb", 'w')
#    handle = m()
#    handle.write.assert_called_once_with(Dseqrecord("GGATCC", circular=True).format())
#
#    s = Dseqrecord("GGATCC", circular=True)
#    m = mock_open(read_data = s.format())
#    with patch('pydna.readers.open')as d:
#        x=read("AAA.gb")
#
#    #"ABC.gb"
#    m.mock_calls
#    m.assert_called_once_with("AAA.gb", 'w')
#    handle = m()
#    handle.write.assert_called_once_with(Dseqrecord("GGATCC", circular=True).format())


def test_cut_args():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC")
    from Bio.Restriction import BamHI, BglII, RestrictionBatch

    rb = RestrictionBatch((BamHI, BglII))
    assert s.cut(BamHI)[0].seq == s.cut(BamHI + BglII)[0].seq == s.cut(rb)[0].seq
    assert s.cut(BamHI)[1].seq == s.cut(BamHI + BglII)[1].seq == s.cut(rb)[1].seq


def test_cut_circular():

    from pydna.dseqrecord import Dseqrecord
    from Bio.Restriction import BsaI, KpnI, Acc65I, NotI

    test = "aaaaaaGGTACCggtctcaaaa"

    for i in range(len(test)):

        nt = test[i:] + test[:i]

        a = Dseqrecord(nt, circular=True).cut(Acc65I)[0]

        assert a.seq.watson.upper() == "GTACCGGTCTCAAAAAAAAAAG"
        assert a.seq.crick.upper() == "GTACCTTTTTTTTTTGAGACCG"
        assert a.seq.ovhg == -4
        b = Dseqrecord(nt, circular=True).cut(KpnI)[0]
        assert b.seq.watson.upper() == "CGGTCTCAAAAAAAAAAGGTAC"
        assert b.seq.crick.upper() == "CTTTTTTTTTTGAGACCGGTAC"
        assert b.seq.ovhg == 4
        c = Dseqrecord(nt, circular=True).cut(BsaI)[0]
        assert c.seq.watson.upper() == "AAAAAAAAAGGTACCGGTCTCA"
        assert c.seq.crick.upper() == "TTTTTGAGACCGGTACCTTTTT"
        assert c.seq.ovhg == -4
        d = Dseqrecord(nt, circular=True).cut(NotI)
        assert d == ()


def test_cut_add():

    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    from Bio.Seq import Seq
    from Bio.Restriction import BamHI, EcoRI, PstI, EcoRV, SmaI

    from Bio.SeqUtils.CheckSum import seguid

    a = Dseqrecord(
        "GGATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGGATCC"
    ).seq
    b = a.cut(BamHI)[1]
    c = Dseqrecord(
        "nCTGCAGtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGAATTCn"
    ).seq
    f, d, l = c.cut((EcoRI, PstI))

    pUC19 = read("pUC19.gb")

    assert pUC19.linear == False

    assert len(pUC19) == 2686
    assert len(pUC19.seq.watson) == 2686
    assert len(pUC19.seq.crick) == 2686

    assert pUC19.seq.circular == True
    assert pUC19.seq.linear == False

    pUC19_SmaI = pUC19.cut(SmaI)
    assert len(pUC19_SmaI) == 1
    pUC19_SmaI = pUC19_SmaI[0]

    assert pUC19_SmaI.linear
    assert len(pUC19_SmaI) == 2686
    assert pUC19_SmaI.linear

    pUC19_SmaI_a = pUC19_SmaI.seq + a

    assert pUC19_SmaI_a.linear
    assert pUC19_SmaI_a.circular == False

    pUC19_SmaI_a = pUC19_SmaI_a.looped()
    assert len(pUC19_SmaI_a) == 2778

    assert pUC19_SmaI_a.circular
    assert pUC19_SmaI_a.linear == False
    assert eq(pUC19_SmaI_a, read("pUC19-SmaI-a.gb"))

    """ sticky end cloning """

    pUC19_BamHI = pUC19.cut(BamHI)

    assert len(pUC19_BamHI) == 1

    pUC19_BamHI = pUC19_BamHI[0].seq

    assert len(pUC19_BamHI.watson) == len(pUC19_BamHI.crick) == 2686

    pUC19_BamHI_a = pUC19_BamHI + b

    assert len(pUC19_BamHI_a.watson) == len(pUC19_BamHI_a.crick) == 2772

    assert pUC19_BamHI_a.circular == False
    assert pUC19_BamHI_a.linear == True

    pUC19_BamHI_a = pUC19_BamHI_a.looped()

    assert pUC19_BamHI_a.circular == True
    assert pUC19_BamHI_a.linear == False

    assert eq(pUC19_BamHI_a, read("pUC19-BamHI-a.gb"))

    pUC19_BamHI_a_rc = pUC19_BamHI + b.rc()

    pUC19_BamHI_a_rc = pUC19_BamHI_a_rc.looped()

    assert pUC19_BamHI_a.circular == True
    assert pUC19_BamHI_a.linear == False
    assert eq(pUC19_BamHI_a_rc, read("pUC19-BamHI-a-rc.gb"))

    """ adding (ligating) dsDNA objects """
    with pytest.raises(TypeError) as excinfo:
        pUC19 + a
    assert "circular" in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        a + pUC19
    assert "circular" in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        a + b
    assert "compatible" in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        b + a
    assert "compatible" in str(excinfo.value)
    with pytest.raises(TypeError) as excinfo:
        d + d
    assert "compatible" in str(excinfo.value)

    """ directional cloning """

    pUC19_EcoRI_PstI = pUC19.cut(EcoRI, PstI)[1]

    with pytest.raises(TypeError) as excinfo:
        pUC19_EcoRI_PstI + d
    assert "compatible" in str(excinfo.value)

    pUC19_EcoRI_PstI_d = pUC19_EcoRI_PstI + d.rc()

    pUC19_EcoRI_PstI_d = pUC19_EcoRI_PstI_d.looped()

    assert eq(pUC19_EcoRI_PstI_d, read("pUC19-EcoRI_PstI-d-rc.gb"))
    assert eq(pUC19_EcoRI_PstI_d.rc(), read("pUC19-EcoRI_PstI-d-rc.gb"))

    from Bio.Restriction import Bsu36I, BstAPI

    pCAPs = read("pCAPs.gb")
    a, b = pCAPs.cut(Bsu36I, BstAPI)
    c = (a + b).looped()
    assert eq(c, pCAPs)


def test_Dseqrecord_cutting_adding_2():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    a = (
        Dseqrecord(
            Dseq(
                "AATTCACANGGTACCNGGTACCNGCGGATATC",
                "GTGTNCCATGGNCCATGGNCGCCTATAG"[::-1],
                -4,
            )
        ),
        Dseqrecord(
            Dseq(
                "CACANGGTACCNGGTACCNGCGGATATC", "GTGTNCCATGGNCCATGGNCGCCTATAG"[::-1], 0
            )
        ),
        Dseqrecord(
            Dseq(
                "CACANGGTACCNGGTACCNGCGGATATC",
                "AATTGTGTNCCATGGNCCATGGNCGCCTATAG"[::-1],
                4,
            )
        ),
    )

    from Bio.Restriction import KpnI, Acc65I, NlaIV

    enzymes = [Acc65I, NlaIV, KpnI]

    for enz in enzymes:
        for f in a:
            b, c, d = f.cut(enz)
            e = b + c + d
            assert str(e.seq).lower() == str(f.seq).lower()


def test_Dseqrecord_cutting_adding_3():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec
    from Bio.Restriction import KpnI, BamHI, Acc65I, NlaIV, EcoRI, EcoRV

    a = read(
        """

LOCUS       New_DNA                   10 bp ds-DNA     linear       02-APR-2013
DEFINITION
ACCESSION   New_DNA
VERSION     New_DNA
KEYWORDS    .
SOURCE
  ORGANISM  . .
COMMENT
COMMENT     ApEinfo:methylated:1
FEATURES             Location/Qualifiers
     misc_feature    1..1
                     /label=1
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    2..2
                     /label=2
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    3..3
                     /label=3
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    4..4
                     /label=4
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    5..5
                     /label=5
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    6..6
                     /label=6
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    7..7
                     /label=7
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    8..8
                     /label=8
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    9..9
                     /label=9
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    10..10
                     /label=10
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
ORIGIN
        1 ttGGTACCgg
//"""
    )

    b, c = a.cut(Acc65I)

    assert [f.qualifiers["label"] for f in b.features] == [
        ["1"],
        ["2"],
        ["3"],
        ["4"],
        ["5"],
        ["6"],
        ["7"],
    ]
    assert [f.qualifiers["label"] for f in c.features] == [
        ["4"],
        ["5"],
        ["6"],
        ["7"],
        ["8"],
        ["9"],
        ["10"],
    ]


def test_Dseqrecord_cutting_adding_4():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec
    from Bio.Restriction import KpnI, Acc65I, NlaIV, EcoRI, EcoRV

    a = read(
        """
LOCUS       New_DNA                   33 bp ds-DNA     linear       08-NOV-2012
DEFINITION  .
ACCESSION   .
VERSION     .
SOURCE      .
  ORGANISM  .
COMMENT
COMMENT     ApEinfo:methylated:1
FEATURES             Location/Qualifiers
     misc_feature    1..11
                     /label=Acc65I-1
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    12..18
                     /label=Acc65I-2
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    19..33
                     /label=Acc65I-3
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    1..15
                     /label=KpnI-1
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    16..22
                     /label=KpnI-2
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    23..33
                     /label=KpnI-3
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    1..13
                     /label=NlaIV-1
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    14..20
                     /label=NlaIV-2
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    21..33
                     /label=NlaIV-3
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
ORIGIN
        1 GAATTCacan ggtaccnGGT ACCngcgGAT ATC
//

    """
    )

    assert a.lseguid() == "di3hL8t2G4iQQsxlm_CtvnUMBz8"

    assert [x.qualifiers["label"][0] for x in a.features] == [
        "Acc65I-1",
        "Acc65I-2",
        "Acc65I-3",
        "KpnI-1",
        "KpnI-2",
        "KpnI-3",
        "NlaIV-1",
        "NlaIV-2",
        "NlaIV-3",
    ]

    b, c, d = a.cut(Acc65I)

    assert [x.qualifiers["label"][0] for x in b.features] == [
        "Acc65I-1",
        "KpnI-1",
        "NlaIV-1",
    ]
    assert [x.qualifiers["label"][0] for x in c.features] == [
        "Acc65I-2",
        "KpnI-2",
        "NlaIV-2",
    ]
    assert [x.qualifiers["label"][0] for x in d.features] == [
        "Acc65I-3",
        "KpnI-3",
        "NlaIV-3",
    ]
    e = b + c + d
    assert sorted([x.qualifiers["label"][0] for x in e.features]) == [
        x.qualifiers["label"][0] for x in a.features
    ]
    assert str(a.seq) == str(e.seq)

    b, c, d = a.cut(KpnI)
    assert [x.qualifiers["label"][0] for x in b.features] == [
        "Acc65I-1",
        "KpnI-1",
        "NlaIV-1",
    ]
    assert [x.qualifiers["label"][0] for x in c.features] == [
        "Acc65I-2",
        "KpnI-2",
        "NlaIV-2",
    ]  # ???
    assert [x.qualifiers["label"][0] for x in d.features] == [
        "Acc65I-3",
        "KpnI-3",
        "NlaIV-3",
    ]
    e = b + c + d
    assert sorted([x.qualifiers["label"][0] for x in e.features]) == [
        x.qualifiers["label"][0] for x in a.features
    ]

    b, c, d = a.cut(NlaIV)
    assert [x.qualifiers["label"][0] for x in b.features] == ["Acc65I-1", "NlaIV-1"]
    assert [x.qualifiers["label"][0] for x in c.features] == ["NlaIV-2"]
    assert [x.qualifiers["label"][0] for x in d.features] == ["KpnI-3", "NlaIV-3"]
    e = b + c + d
    assert str(a.seq) == str(e.seq)

    b, c = a.cut(EcoRI)
    e = b + c
    assert str(a.seq) == str(e.seq)

    b, c = a.cut(EcoRV)
    e = b + c
    assert str(a.seq) == str(e.seq)

    b, c, d = a.cut(EcoRI, EcoRV)
    e = b + c + d

    assert str(a.seq) == str(e.seq)

    b, c, d, f = a.cut(Acc65I, EcoRI)
    e = b + c + d + f
    assert str(a.seq) == str(e.seq)

    b, c, d, f = a.cut(EcoRI, Acc65I)
    e = b + c + d + f
    assert str(a.seq) == str(e.seq)


def test_features_on_slice():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.seqfeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation

    dseq_record = Dseqrecord(Dseq("ACTCTTTCTCTCTCT", circular=True))
    dseq_record.features = [SeqFeature(FeatureLocation(start=2, end=4))]
    assert len(dseq_record[6:1].features) == 0
    assert len(dseq_record[6:3].features) == 0
    assert len(dseq_record[6:4].features) == 1
    assert len(dseq_record[6:5].features) == 1

    dseq_record2 = Dseqrecord(Dseq("ACTCTTTCTCTCTCT", linear=True))
    assert dseq_record2[6:1].seq == Dseq("")


def test_features_change_ori():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec

    s = read(
        """
                    LOCUS       New_DNA                   13 bp ds-DNA     circular     08-JAN-2018
                    DEFINITION  .
                    SOURCE      .
                    COMMENT
                    COMMENT     ApEinfo:methylated:1
                    FEATURES             Location/Qualifiers
                         misc_feature    join(2..3,5..7,9..12)
                                         /locus_tag="hej"
                                         /label="hej"
                                         /ApEinfo_label="hej"
                                         /ApEinfo_fwdcolor="cyan"
                                         /ApEinfo_revcolor="green"
                                         /ApEinfo_graphicformat="arrow_data {{0 1 2 0 0 -1} {} 0}
                                         width 5 offset 0"
                    ORIGIN
                            1 atcaccgatt tta
                    //"""
    )

    for i in range(1, len(s)):
        b = s.shifted(i)
        assert str(b.features[0].extract(b).seq).lower() == "tcccgtttt"

    s = read(
        """
                    LOCUS       New_DNA                   13 bp ds-DNA     circular     08-JAN-2018
                    DEFINITION  .
                    SOURCE      .
                    COMMENT
                    COMMENT
                    COMMENT     ApEinfo:methylated:1
                    FEATURES             Location/Qualifiers
                         misc_feature    complement(join(2..5,7..9,11..12))
                                         /locus_tag="hej"
                                         /label="hej"
                                         /ApEinfo_label="hej"
                                         /ApEinfo_fwdcolor="cyan"
                                         /ApEinfo_revcolor="green"
                                         /ApEinfo_graphicformat="arrow_data {{0 1 2 0 0 -1} {} 0}
                                         width 5 offset 0"
                    ORIGIN
                            1 taaaatcggt gat
                    //"""
    )

    for i in range(1, len(s)):
        b = s.shifted(i)
        assert str(b.features[0].extract(b).seq).lower() == "tcccgtttt"

    s = read(
        """
            LOCUS       New_DNA                   21 bp ds-DNA     circular     03-APR-2013
            DEFINITION  a
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            COMMENT
            COMMENT     ApEinfo:methylated:1
            FEATURES             Location/Qualifiers
                 misc_feature    join(18..21,1..4)
                                 /label=bb
                                 /ApEinfo_fwdcolor=cyan
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
                 misc_feature    5..17
                                 /label=ins
                                 /ApEinfo_fwdcolor=#e03c2b
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
            ORIGIN
                    1 aaaGGTACCt ttGGATCCggg
            //
        """
    )

    assert str(s.features[0].extract(s).seq) == "CGGGAAAG"        # bb
    assert str(s.features[1].extract(s).seq) == "GTACCTTTGGATC"   # ins

    for i in range(1, len(s)):

        b = s.shifted(i)

        assert [
            str(f.extract(b).seq)
            for f in b.features
            if f.qualifiers["label"][0] == "ins"
        ][0] == "GTACCTTTGGATC"
        assert [
            str(f.extract(b).seq)
            for f in b.features
            if f.qualifiers["label"][0] == "bb"
        ][0] == "CGGGAAAG"

    from Bio.Restriction import Acc65I, BamHI

    bb1, ins1 = sorted(s.cut(Acc65I, BamHI), key=len, reverse=True)

    for i in range(1, len(s)):

        b = s.shifted(i)

        bb, ins = sorted(b.cut(Acc65I, BamHI), key=len, reverse=True)

        assert eq(bb1, bb)
        assert eq(ins1, ins)

        assert bb.features[0].extract(bb).seq.watson == "CGGGAAAG"
        assert bb.features[0].extract(bb).seq.crick == "CTTTCCCG"

        assert eq(bb.features[0].extract(bb), s.features[0].extract(s))

        assert ins.features[0].extract(ins).seq.watson == "GTACCTTTG"
        assert ins.features[0].extract(ins).seq.crick == "GATCCAAAG"

        assert str(ins.features[0].extract(ins).seq) == str(
            s.features[1].extract(s).seq
        )


def test_synced():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read

    bb_ins = Dseqrecord("tcgcgcgtttcgATGTAgtgatgacggtgaA", circular=True)

    for i in range(1, len(bb_ins)):
        cand = bb_ins.shifted(i)
        assert str(cand.synced("tcgcgcgtttcg", limit=11).seq).upper() == str(
            bb_ins.seq.upper()
        )
        assert str(cand.synced(Dseq("tcgcgcgtttcg"), limit=11).seq).upper() == str(
            bb_ins.seq.upper()
        )
        assert str(cand.synced(Seq("tcgcgcgtttcg"), limit=11).seq).upper() == str(
            bb_ins.seq.upper()
        )
        assert str(
            cand.synced(SeqRecord(Seq("tcgcgcgtttcg")), limit=11).seq
        ).upper() == str(bb_ins.seq.upper())

    pGUP1 = read("pGUP1_correct.gb")
    pGREG505 = read("pGREG505.gb")
    pGUP1_not_synced = read("pGUP1_not_synced.gb")
    assert (
        pGUP1_not_synced.synced(pGREG505).useguid()
        == "42wIByERn2kSe_Exn405RYwhffU"
        == pGUP1.useguid()
    )

    bb_ins = Dseqrecord("tcgcgcgtttcgAgtgatgacggtgaA", circular=True)

    for i in range(1, len(bb_ins)):
        cand = bb_ins.shifted(i)
        assert str(cand.synced("tcgcgcgtttcg", limit=11).seq).upper() == str(
            bb_ins.seq.upper()
        )

    a = Dseqrecord("gtatcgttctagctattggctagta", circular=True)
    b = Dseqrecord("ggttagtcagttatatcggcttatc", circular=True)

    with pytest.raises(TypeError):
        a.synced(b)

    a = Dseqrecord("gtatcgttctagctattggctagta", circular=False)
    with pytest.raises(TypeError):
        a.synced(b)

    a = Dseqrecord("gtatcgttctagctattggctagta", circular=True)
    b = Dseqrecord("acgatactactagccaatagctaga", circular=True)

    assert str(a.synced(b).seq) == "acgatactactagccaatagctaga"
    assert str(a.synced(a).seq) == "gtatcgttctagctattggctagta"
    assert str(a.rc().synced(a).seq) == "gtatcgttctagctattggctagta"


def test_map_pCR_MCT1_HA46():

    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read

    pCR_MCT1_HA46 = read("pCR_MCT1_HA46.gb")

    slc = pCR_MCT1_HA46.find_aa("VFFKE YPYDVPDYA IEG".replace(" ", ""))

    pCR_MCT1_HA46.map_target = slc

    with pytest.raises(ValueError):
        pCR_MCT1_HA46.map_trace_files("/subfolder/*.ab1")

    map_ = pCR_MCT1_HA46.map_trace_files("*.ab1")

    assert set(map_) == set(
        ["28-1rev_D04_026.ab1", "32-3rev_H04_018.ab1", "36-5rev_D05_041.ab1"]
    )

    assert set([x.fname for x in pCR_MCT1_HA46.matching_reads]) == set(
        ["28-1rev_D04_026.ab1", "32-3rev_H04_018.ab1", "36-5rev_D05_041.ab1"]
    )

    assert set([x.fname for x in pCR_MCT1_HA46.not_matching_reads]) == set(
        ["02-G1_B01_013.ab1"]
    )

    assert pCR_MCT1_HA46.find_aa("YPYDVPDYA".replace(" ", "")) == slice(
        1088, 1115, None
    )

    assert pCR_MCT1_HA46.find_aa("VFFKE YPYDVPDYA IEG".replace(" ", "")) == slice(
        1073, 1124, None
    )

    pCR_MCT1_HA46.add_feature(1073, 1124)
    mt = pCR_MCT1_HA46.features[-1]
    pCR_MCT1_HA46.map_target = mt

    map_ = pCR_MCT1_HA46.map_trace_files("*.ab1")
    assert set(map_) == set(
        ["28-1rev_D04_026.ab1", "32-3rev_H04_018.ab1", "36-5rev_D05_041.ab1"]
    )

    pCR_MCT1_HA46.map_target = None

    map_ = pCR_MCT1_HA46.map_trace_files("*.ab1")

    assert set(map_) == set(
        ["32-3rev_H04_018.ab1", "36-5rev_D05_041.ab1", "28-1rev_D04_026.ab1"]
    )


def test_map_short():
    from pydna.dseqrecord import Dseqrecord

    t = Dseqrecord("AAGTTAAAATAAGGCTAGTCCGTTAT")
    t.map_target = slice(0, 26)
    t.map_trace_files("*.ab1")
    assert t.map_trace_files("*.ab1") == ["02-G1_B01_013.ab1"]


def test_map_too_short():
    from pydna.dseqrecord import Dseqrecord

    t = Dseqrecord("AAGTTAAAATAAGGCTAGTCCGTT")
    # t.map_target = slice(0,23)
    t.map_trace_files("*.ab1")
    assert t.map_trace_files("*.ab1") == []


def test_map_no_match():
    from pydna.dseqrecord import Dseqrecord

    t = Dseqrecord(
        "AGGAGGGGCCCACACCGAGGAAGTAGACTGTTATACGTCGGCGATGGTGGTAGCTAACTATGTTGCCTGCCACTACAACAGTATCTAAGCCGTGTAAAGG"
    )  # random DNA!
    t.map_trace_files("*.ab1")
    assert t.map_trace_files("*.ab1") == []


def test_slicing2():
    from pydna.dseqrecord import Dseqrecord
    from Bio.Restriction import BamHI, EcoRI, KpnI

    a = Dseqrecord("aaGGATCCnnnnnnnnnGAATTCccc", circular=True)
    assert (
        a.cut(
            EcoRI,
            BamHI,
            KpnI,
        )[0].seq
        == a.cut(
            BamHI,
            EcoRI,
            KpnI,
        )[1].seq
    )


def test_rogerstager():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.utils import eq
    from Bio.Seq import Seq
    from Bio.Restriction import BsaI

    answ = []
    answ.append(Dseq("aaaaaaaaaaaaggtctca", "ttttttttccagagttttt"[::-1]))
    answ.append(Dseq("aaaaaaaaaggtctca", "tttttccagagttttt"[::-1]))

    tests = [Seq("aaaaaaggtctcaaaaaaa"), Seq("aaaaaaggtctcaaaa")]

    for s in tests:
        d = Dseqrecord(s).looped()
        for f in d.cut(BsaI):
            a = answ.pop(0)
            assert f.seq.watson == a.watson
            assert f.seq.crick == a.crick
            assert f.seq.ovhg == a.ovhg
            assert eq(f.seq, a)


def test___mul__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    assert s * 3 == Dseqrecord("GGATCCGGATCCGGATCC", circular=False)
    assert s * 0 == Dseqrecord("")
    s = Dseqrecord("GGATCC", circular=True)
    with pytest.raises(TypeError):
        s * 3
    with pytest.raises(TypeError):
        s * 3.1


def test___repr__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    assert repr(s) == "Dseqrecord(-6)"
    s = Dseqrecord("GGATCC", circular=True)
    assert repr(s) == "Dseqrecord(o6)"


def test__repr_pretty_():
    from unittest.mock import MagicMock
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    pp = MagicMock()
    s._repr_pretty_(pp, None)
    pp.text.assert_called_with("Dseqrecord(-6)")
    s = Dseqrecord("GGATCC", circular=True)
    s._repr_pretty_(pp, None)
    pp.text.assert_called_with("Dseqrecord(o6)")


def test___getitem__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    assert s[1:-1].seq == Dseqrecord("GATC", circular=False).seq
    t = Dseqrecord("GGATCC", circular=True)
    assert t[1:5].seq == Dseqrecord("GATC").seq
    #    assert t[1:5].__dict__ == Dseqrecord("GATC").__dict__
    #    assert str(t[1:5].__dict__) == str(Dseqrecord("GATC").__dict__)
    assert s[1:5].seq == Dseqrecord("GATC").seq
    assert s[1:5].seq == Dseqrecord("GATC", circular=False).seq
    assert s[5:1:-1].seq == Dseqrecord("CCTA", circular=False).seq

    assert t[1:1].seq == Dseqrecord("").seq
    assert t[5:1].seq == Dseqrecord("CG", circular=False).seq
    assert t[9:1].seq == Dseqrecord("").seq
    assert t[1:9].seq == Dseqrecord("").seq
    assert t[9:10].seq == Dseqrecord("").seq
    assert t[10:9].seq == Dseqrecord("").seq


def test___eq__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    t = Dseqrecord("GGATCC", circular=False)
    u = Dseqrecord("GGATCC", circular=True)
    v = Dseqrecord("GGATCC", circular=True)

    assert s == t
    assert u == v
    assert s != u

    assert s != 123


def test___ne__():
    pass


def test___hash__():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    t = Dseqrecord("GGATCC", circular=False)
    u = Dseqrecord("GGATCC", circular=True)
    v = Dseqrecord("GGATCC", circular=True)

    assert hash(s) == hash(t)
    assert hash(u) == hash(v)
    assert hash(s) != hash(u)


def test_linearize():
    from Bio.Restriction import BamHI, BglII
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord

    u = Dseqrecord("GGATCC", circular=True)
    frag = u.linearize(BamHI)
    assert frag.seq == Dseq("GATCCG", "GCCTAG"[::-1], -4)
    result = Dseqrecord(Dseq("GATCCG", "GCCTAG"[::-1], -4))
    assert frag.seq == result.seq

    t = Dseqrecord("GGATCC", circular=False)
    with pytest.raises(TypeError):
        t.linearize(BamHI)
    with pytest.raises(TypeError):
        u.linearize(BglII)
    u = Dseqrecord("GGATCCagatct", circular=True)
    u.linearize(BglII)
    u.linearize(BamHI)
    with pytest.raises(TypeError):
        u.linearize(BamHI, BglII)
    with pytest.raises(TypeError):
        u.linearize(BamHI + BglII)
    u = Dseqrecord("GGATCCGGATCC", circular=True)
    with pytest.raises(TypeError):
        u.linearize(BamHI)


def test_cutters():
    from pydna.dseqrecord import Dseqrecord

    obj = Dseqrecord("ggatcc")
    from Bio.Restriction import BglII, BamHI

    assert BglII in obj.no_cutters()
    assert BamHI not in obj.no_cutters()

    assert BamHI in obj.unique_cutters()

    assert BamHI in obj.once_cutters()

    assert BamHI in (obj + obj).twice_cutters()
    assert BamHI not in obj.twice_cutters()

    assert BamHI in obj.n_cutters(1)
    assert BamHI in obj.cutters()


def test_number_of_cuts():
    from Bio.Restriction import EcoRI, BamHI, BglII
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCCgaattc", circular=False)
    s.cut(BamHI)
    s.cut(EcoRI)
    s.cut(BglII)
    from Bio.Restriction import RestrictionBatch

    rb = RestrictionBatch((EcoRI, BamHI, BglII))
    assert s.number_of_cuts(BamHI) == 1
    assert s.number_of_cuts(EcoRI) == 1
    assert s.number_of_cuts(BglII) == 0
    assert s.number_of_cuts(rb) == 2
    assert s.number_of_cuts(EcoRI + BamHI) == 2


def test_reverse_complement():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCC", circular=False)
    assert s.reverse_complement() == s.rc()
    assert s.rc().seq == s.seq
    assert s.reverse_complement().seq == s.seq

    from pydna.readers import read

    s = read(
        """
                LOCUS       New_DNA                   13 bp ds-DNA     circular     12-NOV-2013
                DEFINITION  .
                ACCESSION
                VERSION
                SOURCE      .
                ORGANISM  .
                COMMENT
                COMMENT     ApEinfo:methylated:1
                FEATURES             Location/Qualifiers
                     misc_feature    join(9..10,12..13,1..1,3..6)
                                     /label=hej
                                     /ApEinfo_fwdcolor=cyan
                                     /ApEinfo_revcolor=green
                                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                     width 5 offset 0
                ORIGIN
                        1 gattttaatc acc
                //"""
    )


def test_shifted():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("GGATCCgaattc", circular=False)
    with pytest.raises(TypeError):
        s.shifted(1)


def test_looped():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord

    a = Dseqrecord("aAAa")
    a.add_feature()
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord("aAAa")
    a.add_feature(0, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord("aAAa")
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord("aAAa")
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord("aAAa")
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features






    a = Dseqrecord(Dseq("gAAa", "cTTt", ovhg=-1))
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features[0].extract(a).seq == b.features[0].extract(b).seq
    assert str(a.features[0].extract(a).seq) == str(b.features[0].extract(b).seq)
    assert a.features == b.features

    a = Dseqrecord(Dseq("caaa", "gttt", ovhg=-1))
    a.add_feature(0, 5)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("caaa", "gttt", ovhg=-1))
    a.add_feature(0, 5, strand=-1)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("aaac", "tttg", ovhg=1))
    a.add_feature(2, 4)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("aaaa", "tttt", ovhg=1))
    a.add_feature(0, 5)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("aaaa", "tttt", ovhg=1))
    a.add_feature(0, 5, strand=-1)
    b = a.looped()
    assert a.features == b.features

    a = Dseqrecord(Dseq("aaaa", "tttt", ovhg=-1))
    a.add_feature(0, 6)
    b = a.looped()
    assert a.features == b.features


def test_upper():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("Gc")
    s.annotations["sample"] = ["sample"]
    u = s.upper()
    assert s.upper() == u
    assert s.seq.upper() == u.seq.upper()
    del u.__dict__["_seq"]
    del s.__dict__["_seq"]
    assert u.__dict__ == s.__dict__


def test_lower():
    from pydna.dseqrecord import Dseqrecord

    s = Dseqrecord("Gc")
    s.annotations["sample"] = ["sample"]
    l = s.lower()
    assert s.lower() == l
    assert s.seq.upper() == l.seq.upper()
    del l.__dict__["_seq"]
    del s.__dict__["_seq"]
    assert l.__dict__ == s.__dict__


def test_map():
    from pydna.dseq import Dseq
    from pydna.dseqrecord import Dseqrecord
    from pydna.readers import read
    from pydna.utils import eq

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord as Srec
    from Bio.SeqIO import read as abiread

    traces = []

    import glob

    for name in glob.glob("*.ab1"):
        traces.append(abiread(name, "abi"))

    for t in traces:

        d = Dseqrecord(t.seq)

        if "ITVFFKEYPYDVPDYAIEGIFHAT" in d:

            tag = "tat cca tat gac gtt cca gac tat gca".replace(" ", "")
            trc = "ata ggt ata ctg caa ggt ctg ata cgt"[::-1].replace(" ", "")

            s = Dseqrecord(Dseq(tag, trc, 0))
            sl = s.find_aa("YPYDVPDYA")
            assert str(s[sl].seq.translate()) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s

            tag = "AAA tat cca tat gac gtt cca gac tat gca".replace(" ", "")
            trc = "    ata ggt ata ctg caa ggt ctg ata cgt"[::-1].replace(" ", "")

            s = Dseqrecord(Dseq(tag, trc, 0))
            sl = s.find_aa("YPYDVPDYA")
            assert str(s[sl].seq.translate()) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s

            tag = "    tat cca tat gac gtt cca gac tat gca".replace(" ", "")
            trc = "AAA ata ggt ata ctg caa ggt ctg ata cgt"[::-1].replace(" ", "")

            s = Dseqrecord(Dseq(tag, trc, 0))
            sl = s.find_aa("YPYDVPDYA")
            assert str(s[sl].seq.translate()) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s

            tag = "    tat cca tat gac gtt cca gac tat gca".replace(" ", "")
            trc = "AAA ata ggt ata ctg caa ggt ctg ata cgt"[::-1].replace(" ", "")

            s = Dseqrecord(Dseq(tag, trc, 0))
            sl = s.find_aa("YPYDVPDYA")
            assert str(s[sl].seq.translate()) == "YPYDVPDYA"

            tag = "tat cca tat gac gtt cca gac tat gca".replace(" ", "")
            trc = "ata ggt ata ctg caa ggt ctg ata cgt"[::-1].replace(" ", "")

            tag, trc = trc, tag

            s = Dseqrecord(Dseq(tag, trc, 0))
            sl = s.rc().find_aa("YPYDVPDYA")

            assert str(s.rc()[sl].seq.translate()) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s.rc()

            tag = "aaa tat cca tat gac gtt cca gac tat gca".replace(" ", "")
            trc = "ttt ata ggt ata ctg caa ggt ctg ata cgt"[::-1].replace(" ", "")

            s = Dseqrecord(Dseq(tag, trc, 0, circular=True))
            sl = s.find_aa("YPYDVPDYA")
            assert str(s[sl].seq.translate()) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s


def test_assemble_YEp24PGK_XK():
    """
    test YEp24PGK_XK
    """
    import pytest
    import sys

    from pydna.readers import read
    from pydna.utils import eq
    from pydna.amplify import pcr

    """ test YEp24PGK_XK"""
    p1 = read("primer1.txt", ds=False)
    p3 = read("primer3.txt", ds=False)
    XKS1 = read("XKS1_orf.txt")
    YEp24PGK = read("YEp24PGK.txt")

    PCR_prod = pcr(p1, p3, XKS1)

    from Bio.Restriction import BamHI

    stuffer1, insert, stuffer2 = PCR_prod.cut(BamHI)

    from Bio.Restriction import BglII

    YEp24PGK_BglII = YEp24PGK.cut(BglII)[0]

    YEp24PGK_XK = YEp24PGK_BglII + insert

    YEp24PGK_XK = YEp24PGK_XK.looped()

    YEp24PGK_XK = YEp24PGK_XK.synced(
        "gaattctgaaccagtcctaaaacgagtaaataggaccggcaattc"
    )  # YEp24PGK)

    assert YEp24PGK_XK.useguid() == "HRVpCEKWcFsKhw_W-25ednUfldI"
    assert YEp24PGK_XK.cseguid() == "t9fs_9UvEuD-Ankyy8XEr1hD5DQ"

    YEp24PGK_XK_correct = read("YEp24PGK_XK_manually_assembled.txt")
    assert YEp24PGK_XK_correct.cseguid() == "t9fs_9UvEuD-Ankyy8XEr1hD5DQ"
    assert eq(YEp24PGK_XK, YEp24PGK_XK_correct)


if __name__ == "__main__":
    args = [
    __file__,
    "--cov=pydna",
    "--cov-append",
    "--cov-report=html:../htmlcov",
    "--cov-report=xml",
    "--capture=no",
    "--durations=10",
    "--import-mode=importlib",
    "--nbval",
    "--current-env",
    "--doctest-modules",
    "--capture=no",
    "-vvv"]
    pytest.main(args)
