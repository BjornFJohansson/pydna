#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

def test_pydna_gbtext_clean():

    pytest.importorskip("pyparsing")

    from pydna.readers import read
    from pydna.genbankfixer import gbtext_clean

    files = [
        ("sequence.gb", "j2yAlBCZ-txSTCkakAmykAielRI"),
        ("NCBI_example.gb", "j2yAlBCZ-txSTCkakAmykAielRI"),
        ("YEplac181.txt", "lbnQtxi5LyDONRswRdG88-l8NF0"),
        ("pGADT7-Rec.gb", "rhnYE78wGKdqAZWiyVJzQ7HXXys"),
        ("P30350%20(2013-10-11%2013_49_14).dna.txt", "_aEPoGLctHcOZdQdZIh-KyBt5WY"),
        ("ApE_example.gb", "c47i2ifiNZVvvnLQbX5anTVVoPE"),
        ("VectorNTI_example.gb", "bDPbx5P4yigGWh1zK7FiG_SF8qQ"),
        ("hej.txt", "lbnQtxi5LyDONRswRdG88-l8NF0"),
        ("fakeGenBankFile.gb", "ATrCXrjheFhltm8HhLJuFNtWXGw"),
    ]

    for file_, seg in files:
        with open("broken_genbank_files/" + file_, "r") as f:
            infile = f.read()
        if file_ == "hej.txt":
            from Bio import BiopythonParserWarning

            with pytest.warns(BiopythonParserWarning):
                assert read(gbtext_clean(infile).gbtext).useguid() == seg
        else:
            assert read(gbtext_clean(infile).gbtext).useguid() == seg


def test_wrapstring():

    pytest.importorskip("pyparsing")

    from pydna.genbankfixer import wrapstring

    assert (
        wrapstring(
            "0123456789",
            0,
            5,
        )
        == "01234\n56789\n\n"
    )
    assert wrapstring("0123456789", 1, 5, padfirst=False) == "0123\n 4567\n 89\n"


if __name__ == "__main__":
    pytest.main([__file__, "-x", "-vv", "-s"])
