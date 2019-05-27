#!/usr/bin/env python
# -*- coding: utf-8 -*-
# https://www.ncbi.nlm.nih.gov/nuccore/5

import pytest
from unittest import mock

def test_set_email():
    from pydna.genbank import Genbank
    with pytest.raises(ValueError):
        gb = Genbank("invalidemailaddress")
    with pytest.raises(ValueError):
        gb = Genbank("someone@example.com")

def test_repr():
    from pydna.genbank import Genbank
    gb = Genbank("someoneelse@example.com")
    assert repr(gb) == 'GenbankConnection(someoneelse@example.com)'

@mock.patch('Bio.Entrez._urlopen')
def test_biopython_Entrez_efetch_whole(urlopenMock,monkeypatch):
    from Bio import Entrez
    urlopenMock.return_value = open("X60065.gb", "rb")
    Entrez._urlopen = urlopenMock
    Entrez.email = "bjornjobb@gmail.com"
    Entrez.tool  = "pydna"
    handle = Entrez.efetch(db="nuccore",
                           id="X60065.1",
                           rettype="gb",
                           retmode="text")
    from Bio import SeqIO
    result = SeqIO.read(handle, "genbank")
    canned = SeqIO.read("X60065.gb", "genbank")
    assert str(result.seq) == str(canned.seq)

@mock.patch('Bio.Entrez._urlopen')
def test_pydna_Genbank_fresh(urlopenMock,monkeypatch):
    from Bio import Entrez
    from pydna.genbank import Genbank
    urlopenMock.return_value = open("X60065.gb", "rb")
    Entrez._urlopen = urlopenMock
    monkeypatch.setenv("pydna_cached_funcs", "")
    gb = Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("X60065.1")
    from Bio import SeqIO
    canned = SeqIO.read("X60065.gb", "genbank")
    assert str(result.seq) == str(canned.seq)

@mock.patch('Bio.Entrez._urlopen')
def test_pydna_Genbank_from_cache(urlopenMock,monkeypatch):
    from Bio import Entrez
    from pydna.genbank import Genbank
    urlopenMock.return_value = open("X60065.gb", "rb")
    Entrez._urlopen = urlopenMock
    monkeypatch.setenv("pydna_cached_funcs", "pydna.genbank.Genbank.nucleotide")
    gb = Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("X60065.1")
    from Bio import SeqIO
    canned = SeqIO.read("X60065.gb", "genbank")
    assert str(result.seq) == str(canned.seq)
    urlopenMock.return_value = None
    result = gb.nucleotide("X60065.1")
    assert str(result.seq) == str(canned.seq)

def test_genbank_function_set_email(monkeypatch):
    from pydna.genbank import Genbank
    mock_Gb = mock.MagicMock()
    monkeypatch.setenv("pydna_cached_funcs", "")
    monkeypatch.setenv("pydna_email", "someoneelse@example.com")
    monkeypatch.setattr("pydna.genbank.Genbank", mock_Gb)
    from pydna.genbank import genbank
    s=genbank("X60065")
    mock_Gb.assert_called_with("someoneelse@example.com")


def test_pydna_Genbank_fresh_part(monkeypatch):
    monkeypatch.setenv("pydna_cached_funcs", "")
    import pytest
    from unittest import mock
    mock_efetch = mock.MagicMock(name="mock_efetch1")
    mock_efetch().read.side_effect = open("X60065-100-110.gb", "r").read
    monkeypatch.setenv("pydna_email", "someoneelse@example.com")
    monkeypatch.setattr("pydna.genbank._Entrez.efetch", mock_efetch)
    from pydna.genbank import Genbank
    gb = Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("X60065.1", seq_start=1,seq_stop=10)
    assert(str(result.seq).lower() == "ctgaaacggac")


def test_pydna_Genbank_fresh_partII(monkeypatch):
    monkeypatch.setenv("pydna_cached_funcs", "")
    import pytest
    from unittest import mock
    mock_efetch = mock.MagicMock(name="mock_efetch1")
    mock_efetch().read.side_effect = open("X60065-100-110.gb", "r").read
    monkeypatch.setenv("pydna_email", "someoneelse@example.com")
    monkeypatch.setattr("pydna.genbank._Entrez.efetch", mock_efetch)
    from pydna.genbank import Genbank
    gb = Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("X60065.1 REGION: 100..110")
    assert(str(result.seq).lower() == "ctgaaacggac")

@mock.patch('Bio.Entrez._urlopen')
def test_pydna_Genbank_fresh_circular(urlopenMock,monkeypatch):
    from Bio import Entrez
    from pydna.genbank import Genbank
    urlopenMock.return_value = open("pUC19.gb", "rb")
    Entrez._urlopen = urlopenMock
    monkeypatch.setenv("pydna_cached_funcs", "")
    gb = Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("L09137.2")
    from Bio import SeqIO
    canned = SeqIO.read("pUC19.gb", "genbank")
    assert str(result.seq) == str(canned.seq)
    assert result.circular
    assert not result.linear

@mock.patch('Bio.Entrez._urlopen')
def test_pydna_Genbank_set_strand(urlopenMock,monkeypatch):
    from Bio import Entrez
    from pydna.genbank import Genbank
    urlopenMock.return_value = open("X60065.gb", "rb")
    Entrez._urlopen = urlopenMock
    monkeypatch.setenv("pydna_cached_funcs", "")
    gb = Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("X60065.1", strand=1)
    from Bio import SeqIO
    canned = SeqIO.read("X60065.gb", "genbank")
    assert str(result.seq) == str(canned.seq)

@mock.patch('Bio.Entrez._urlopen')
def test_pydna_Genbank_set_strand_not_valid(urlopenMock,monkeypatch):
    from Bio import Entrez
    from pydna.genbank import Genbank
    urlopenMock.return_value = open("X60065.gb", "rb")
    Entrez._urlopen = urlopenMock
    monkeypatch.setenv("pydna_cached_funcs", "")
    gb = Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("X60065.1", strand="notvalid")
    from Bio import SeqIO
    canned = SeqIO.read("X60065.gb", "genbank")
    assert str(result.seq) == str(canned.seq)

@mock.patch('Bio.Entrez._urlopen')
def test_pydna_Genbank_set_strand_antisense(urlopenMock,monkeypatch):
    from Bio import Entrez
    from pydna.genbank import Genbank
    urlopenMock.return_value = open("X60065.gb", "rb")
    Entrez._urlopen = urlopenMock
    monkeypatch.setenv("pydna_cached_funcs", "")
    gb = Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("X60065.1", strand="antisense")
    from Bio import SeqIO
    canned = SeqIO.read("X60065.gb", "genbank")
    assert str(result.seq) == str(canned.seq)

if __name__ == '__main__':
    pytest.main([__file__, "-vv", "-s", "--cov=pydna","--cov-report=html"])



