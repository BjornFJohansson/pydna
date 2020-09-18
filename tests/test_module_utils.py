#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from unittest import mock


def test_eq():

    from pydna.dseqrecord import Dseqrecord

    from pydna.utils import eq
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    assert eq("AAA", "TTT", linear=True)
    assert eq("AAA", "TTT", linear=False)

    assert eq("aAA", "TtT", linear=True)
    assert eq("AAa", "TtT", linear=False)

    assert eq("ATA", "AAT", circular=True)
    assert not eq("ATA", "AAT", circular=False)
    assert eq("AAA", "AAA", linear=True)
    assert eq("AAA", "AAA", linear=False)

    assert eq("ATA", Seq("AAT"), circular=True)
    assert not eq("ATA", Seq("AAT"), circular=False)
    assert eq("AAA", Seq("AAA"), linear=True)
    assert eq("AAA", Seq("AAA"), linear=False)

    assert eq("ATA", SeqRecord("AAT"), circular=True)
    assert not eq("ATA", SeqRecord("AAT"), circular=False)
    assert eq("AAA", SeqRecord("AAA"), linear=True)
    assert eq("AAA", SeqRecord("AAA"), linear=False)

    assert eq("ATA", Dseqrecord("AAT"), circular=True)
    assert not eq("ATA", Dseqrecord("AAT"), circular=False)
    assert eq("AAA", Dseqrecord("AAA"), linear=True)
    assert eq("AAA", Dseqrecord("AAA"), linear=False)

    assert eq(Seq("ATA"), SeqRecord("AAT"), circular=True)
    assert not eq(Seq("ATA"), SeqRecord("AAT"), circular=False)
    assert eq(Seq("AAA"), SeqRecord("AAA"), linear=True)
    assert eq(Seq("AAA"), SeqRecord("AAA"), linear=False)

    assert eq(Seq("ATA"), Dseqrecord("AAT"), circular=True)
    assert not eq(Seq("ATA"), Dseqrecord("AAT"), circular=False)
    assert eq(Seq("AAA"), Dseqrecord("AAA"), linear=True)
    assert eq(Seq("AAA"), Dseqrecord("AAA"), linear=False)

    assert eq(Dseqrecord("AAA", circular=False), Dseqrecord("AAA", circular=False))
    assert eq(Dseqrecord("AAA", circular=True), Dseqrecord("AAA", circular=True))
    assert not eq(Dseqrecord("ATA", circular=False), Dseqrecord("AAT", circular=False))
    assert eq(Dseqrecord("ATA", circular=True), Dseqrecord("AAT", circular=True))

    with pytest.raises(ValueError):
        eq(Dseqrecord("ATA", circular=True), Dseqrecord("ATA", circular=False))

    assert not eq(Dseqrecord("ATA", circular=True), Dseqrecord("ATAA", circular=True))

    assert eq(Dseqrecord("ATA"), Dseqrecord("ATA"), circular=True)
    assert not eq(Dseqrecord("ATA"), Dseqrecord("CCC"), circular=True)
    assert not eq(
        Dseqrecord("ATA"), Dseqrecord("ATA"), Dseqrecord("CCC"), circular=True
    )


# def test_shift_origin():
#    from pydna.readers import read
#    from pydna.dseqrecord import Dseqrecord
#
#
#    from pydna.utils import shift_origin, eq
#    from Bio.Seq import Seq
#    from Bio.SeqRecord import SeqRecord
#    pCAPs   = read("pCAPs.gb")
#    assert pCAPs.circular
#    pCAPs_b = shift_origin(pCAPs, 200)
#    assert len(pCAPs) == len(pCAPs_b)
#    assert pCAPs_b.circular
#    assert eq(pCAPs, pCAPs_b)
#    pCAPs_b_linear = pCAPs_b.tolinear()
#    assert eq(pCAPs, pCAPs_b_linear, circular=True)
#    pCAPs_c = pCAPs[200:]+pCAPs[:200]
#    assert eq(pCAPs, pCAPs_c, circular=True)
#    #with self.assertRaisesRegex(ValueError, "shift"):
#    #    pCAPs_b = shift_origin(pCAPs, 20000)


# def test_copy_features():
#    from pydna.readers import read
#    from pydna.utils import seguid, copy_features
#
#    a=read("pCAPs.gb")
#    b=read("pCAPs_fasta.txt")
#
#    for sh in [1,2,3,3127,3128,3129]:
#        newb = (b[sh:]+b[:sh]).looped()
#        copy_features(a, newb)
#        #print "a",[len(str(f.extract(a).seq.lower()) for f in a.features if len(f)>10]
#        #print "b",[len(str(f.extract(newb).seq).lower()) for f in newb.features]
#        x= sorted([str(f.extract(a).seq).lower() for f in a.features if len(f)>10],key=len)
#        y= sorted([str(f.extract(newb).seq).lower() for f in newb.features],key=len)
#        assert x==y
#
#    b=b.rc()
#
#    for sh in [1,2,3,3127,3128,3129]:
#        newb = b[sh:]+b[:sh]
#        copy_features(a, newb)
#
#
#        x = sorted([str(f.extract(a).seq).lower() for f in a.features if len(f)>10],key=len)
#        y = sorted([str(f.extract(newb).seq).lower() for f in newb.features],key=len)
#        assert x==y
#
#    seguid_bla = "riT98j2v4NxVS8sbw_Q8epCwQwo"
#    seguid_cre = "xLZ2xs2O8CUMmWh2OrhmNFp5ZLg"
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre, seguid_cre, seguid_bla, seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#
#    b=b.synced("attaacgagtgccgtaaacgacgatggttttacc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("ttaacgagtgccgtaaacgacgatggttttacc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("taacgagtgccgtaaacgacgatggttttacc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("gttaccaatgcttaatcagtgaggcacctatctcagc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("ttaccaatgcttaatcagtgaggcacctatctcagc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]
#
#    b=read("pCAPs_fasta.txt").looped()
#    b=b.synced("taccaatgcttaatcagtgaggcacctatctcagc")
#
#    copy_features(a, b)
#    assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,]


def test_cseguid():
    from pydna.utils import cseguid

    x = "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"
    assert (
        cseguid(x)
        == cseguid(x.upper())
        == cseguid(x.lower())
        == "naaZmDzyMa58OsNXROe5SvjC7WU"
    )


def test_lseguid():
    from pydna.utils import lseguid

    x = "tcgcgcgtttcggtgatgacggtgAAAAcctctgacacatgcagctcccggattgtactgagagtgc"
    assert (
        lseguid(x)
        == lseguid(x.upper())
        == lseguid(x.lower())
        == "bHrqalTJ793oAigMQ5_qCttJRTk"
    )


def test_seq31():
    pass


def test_parse_text_table():
    pass


def test_join_list_to_table():
    pass


def test_expandtolist():
    pass


def test_randomRNA():
    pass


def test_randomDNA():
    pass


def test_randomORF():
    pass


def test_randomprot():
    pass


def test_smallest_rotation():
    from pydna.utils import SmallestRotation as sr

    assert sr("tttaaa") == "aaattt"


def test_memorize(monkeypatch):
    import pytest
    from unittest import mock

    from pydna.utils import memorize as _memorize

    @_memorize("mf")
    def mf(*args, **kwargs):
        return args, kwargs

    import base64 as _base64
    import pickle as _pickle
    import hashlib as _hashlib

    args = (1,)
    kwargs = {"kw": 1}

    dump = _pickle.dumps((args, kwargs))

    hash_ = _hashlib.sha1(dump).digest()

    bkey = _base64.urlsafe_b64encode(hash_)

    key = bkey.decode("ascii")

    assert key == "6pHTTwgXP8xcXoEMEzdKSzN6EeM=" or "ux_W9TiWkWBAkQD_FgZTO-pXuYk="

    class Fakedict(dict):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def close(self):
            pass

    cache = Fakedict()
    cache[key] = "saved!"
    mockshelve_open = mock.MagicMock()
    mockshelve_open.return_value = cache

    monkeypatch.setenv("pydna_cached_funcs", "mf")
    monkeypatch.setattr("pydna.utils._shelve.open", mockshelve_open)

    monkeypatch.setenv("pydna_cached_funcs", "mf")

    assert mf(1, kw=1) == "saved!"

    cache[key] = ((1,), {"kw": 1})

    assert mf(1, kw=1) == ((1,), {"kw": 1})

    assert mf(2, kw=2) == ((2,), {"kw": 2})

    monkeypatch.setenv("pydna_cached_funcs", "")

    assert mf(1, kw=1) == ((1,), {"kw": 1})


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s", "--cov=pydna", "--cov-report=html"])
