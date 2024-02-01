#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_pydna_gbtext_clean():
    pytest.importorskip("pyparsing")

    from pydna.readers import read
    from pydna.genbankfixer import gbtext_clean

    files = [
        ("sequence.gb", "ldseguid-iHVblovQ8vMlc6r4I0WtJ7U86FE"),
        ("NCBI_example.gb", "ldseguid-iHVblovQ8vMlc6r4I0WtJ7U86FE"),
        ("YEplac181.txt", "ldseguid-Ko6R38rKGk-BM3cW-Yd3oSP1YII"),
        ("pGADT7-Rec.gb", "cdseguid-uR971QAEetpd7Mptme2A5yybbek"),
        ("P30350%20(2013-10-11%2013_49_14).dna.txt", "cdseguid-C9vd9PyLFUV4LF8llPM3DgM8sA4"),
        ("ApE_example.gb", "ldseguid-5v2YtLDJRhA0_O_Ut6HUXr_4EIc"),
        ("VectorNTI_example.gb", "cdseguid-_nYOPeDWBC7OaB2Arnt5x5fAFZM"),
        ("hej.txt", "ldseguid-Ko6R38rKGk-BM3cW-Yd3oSP1YII"),
        ("fakeGenBankFile.gb", "cdseguid-mIbMOKcoiRlnUs7ZaYwBxQnFego"),
    ]

    for file_, seg in files:
        with open("broken_genbank_files/" + file_, "r") as f:
            infile = f.read()
        if file_ == "hej.txt":
            from Bio import BiopythonParserWarning

            with pytest.warns(BiopythonParserWarning):
                assert read(gbtext_clean(infile).gbtext).seguid() == seg
        else:
            assert read(gbtext_clean(infile).gbtext).seguid() == seg


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
    pytest.main([__file__, "-x", "-vvv", "-s"])
