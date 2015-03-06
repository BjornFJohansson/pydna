#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nose
from pydna import parse, read

from Bio.SeqIO import read as BPread
from Bio.SeqIO import parse as BPparse

def test_bp():

    q = BPread("read1.gb", "gb")
    w = BPread("read2.gb", "gb")
    e = BPread("read3.fasta", "fasta")    
    r = BPread("read4.fasta", "fasta") 
    
    q.format("gb")
    w.format("gb")
    
    v, b = BPparse("pth1.txt", "gb")
    
    print type(v), type(b)
    print type(v.format("gb")), type(b.format("gb"))
    

def test_parse_from_file():
    
    v, b = parse("pth1.txt")
    
    print type(v), type(b)
    
    print type(v), type(b)
    print type(v.format("gb")), type(b.format("gb"))


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
    
    input_ =u'''
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
    
    input_ =u'''>hej
                acgt'''
    assert str(a.seq)=="ACGT"

    input_ =u'''>hej öööh!
                acgt'''
    assert str(a.seq)=="ACGT"                       
    
    input_ =u'''
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
    

def test_read_from_file():
    a = read("./read1.gb") 
    b = read("./read2.gb")
    c = read("./read3.fasta")  
    d = read("./read4.fasta")

    a.format("gb")
    b.format("gb")
    c.format("gb")
    d.format("gb")
    
    assert str(a.seq).lower()==str(b.seq).lower()==str(c.seq).lower()==str(d.seq).lower()     
    
if __name__ == '__main__':
    nose.runmodule()

