#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest
import pydna

def test_efetch_download_text():
    # see https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=E05006&strand=1&rettype=gb&retmode=gbwithparts"
    cachevar = os.environ["pydna_cache"]
    gbdata = pydna.download_text(url)
    os.environ["pydna_cache"] = cachevar
    with open("E05006.gb") as f:
        localdata = f.read().strip()
    assert localdata==gbdata

def test_biopython_download():
    from Bio import Entrez
    from Bio import SeqIO    
    Entrez.email = "bjornjobb@gmail.com"    
    handle = Entrez.efetch(db="nuccore",
                           id="E05006",
                           rettype="gb",
                           retmode="text")    
    result = SeqIO.read(handle, "genbank")    
    assert str(result.seq) == "ATATGGGTACCGATCGTACGGACCA"

def test_pydna_download_fresh():
    cachevar = os.environ["pydna_cache"]
    os.environ["pydna_cache"] = "nocache"
    gb = pydna.Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("E05006")
    assert len(result) == 25
    assert str(result.seq) == "ATATGGGTACCGATCGTACGGACCA"
    os.environ["pydna_cache"] = cachevar
    
def test_pydna_download_cache():
    cachevar = os.environ["pydna_cache"]
    os.environ["pydna_cache"] = "cached"
    gb = pydna.Genbank("bjornjobb@gmail.com")
    result = gb.nucleotide("E05006")
    assert len(result) == 25
    assert str(result.seq) == "ATATGGGTACCGATCGTACGGACCA"
    result = gb.nucleotide("E05006")
    assert len(result) == 25
    assert str(result.seq) == "ATATGGGTACCGATCGTACGGACCA"    
    os.environ["pydna_cache"] = cachevar 
    

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
    

    
