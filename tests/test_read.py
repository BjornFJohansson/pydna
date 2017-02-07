#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import sys
from pydna.readers import read
from pydna.parsers import parse

from Bio.SeqIO import read as BPread
from Bio.SeqIO import parse as BPparse

def pydna_read_test():
    print("sys.getdefaultencoding()",sys.getdefaultencoding())
    import locale
    print("locale.getpreferredencoding()", locale.getpreferredencoding())
    assert read("pydna_read_test.txt").format("gb")[559:578] == '/label="2micron 2µ"'

def test_parse_and_read_with_biopython_and_pydna():

    q = BPread("read1.gb", "gb")
    w = BPread("read2.gb", "gb")
    e = BPread("read3.fasta", "fasta")
    r = BPread("read4.fasta", "fasta")

    #a, b = BPparse("pth1.txt", "gb")
    with open("pth1.txt", "r", encoding="utf-8") as f:
        a, b = BPparse(f, "gb")
        
    print("|"+a.features[13].qualifiers['label'][0]+"|")
    print("|"+a.format("gb")[3314:3324]+"|")
    
    assert a.features[13].qualifiers['label'][0] == '2micron 2µ'
    assert a.format("gb")[3314:3324] == '2micron 2µ'

    x, y = parse("pth1.txt")

    assert "".join(a.format("gb").splitlines()[1:]) == "".join(x.format("gb").splitlines()[1:])
    assert "".join(b.format("gb").strip().splitlines()[4:]) == "".join(y.format("gb").splitlines()[4:])

def test_read_from_string():

    input_ ='''
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
            '''
    a = read(input_)
    
    assert str(a.seq)=="ACGT"

    input_ ='''>hej
               acgt'''
               
    assert str(a.seq)=="ACGT"

    input_ ='''
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
            '''
    a = read(input_)
    assert str(a.seq)=="ACGT"

    input_ ='''>hej
                acgt'''
    assert str(a.seq)=="ACGT"

    input_ ='''>hej öööh!
                acgt'''
    assert str(a.seq)=="ACGT"

    input_ ='''
                LOCUS       New_DNA                    4 bp ds-DNA     linear       30-MAR-2013
                DEFINITION  öööh!
                ACCESSION
                VERSION
                SOURCE      .
                  ORGANISM  .
                COMMENT
                COMMENT     ApEinfo:methylated:1
                FEATURES             Location/Qualifiers
                     misc_feature    2..3
                                     /label=öööh!
                                     /ApEinfo_fwdcolor=cyan
                                     /ApEinfo_revcolor=green
                                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                     width 5 offset 0
                ORIGIN
                        1 acgt
                //
            '''
    a = read(input_)
    assert str(a.seq)=="ACGT"

def test_read_from_unicode():
    with open("pth1.txt", "r", encoding="utf-8") as f: 
        text = f.read()
    assert type(text) == str
    x,y = parse( text )
    assert x.format()[3314:3324] == '2micron 2µ'

def test_read_from_file():
    a = read("read1.gb")
    b = read("read2.gb")
    c = read("read3.fasta")
    d = read("read4.fasta")
    x,y = parse( "pth1.txt" )

    a.format("gb")
    b.format("gb")
    c.format("gb")
    d.format("gb")
    x.format("gb")
    y.format("gb")
    assert x.format()[3314:3324] == '2micron 2µ'
    assert x.features[13].qualifiers['label'][0] == u'2micron 2µ'
    assert str(a.seq).lower()==str(b.seq).lower()==str(c.seq).lower()==str(d.seq).lower()

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])

