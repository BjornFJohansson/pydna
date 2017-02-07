#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test parse
'''

import pytest
import sys
import os

from pydna.assembly   import Assembly
from pydna.dseqrecord import Dseqrecord
from pydna.parsers    import parse

from pydna.design import primer_design
from pydna.design import assembly_fragments

frags = parse('''
>49
ccaaggacacaatcgagctccgatccgtactgtcgagaaacttgtatcc
>50
ctctaactagtatggatagccgtgtcttcactgtgctgcggctacccatc
>51
gtagtgaaacatacacgttgctcgggttcaccccggtccgttctgagtcga
''')

from Bio.Restriction import BamHI
bam = Dseqrecord(BamHI.site)




def test_primer_design_all_pcr_products():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments(x, 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]
    assert result.seq == (frags[0]+frags[1]+frags[2]).seq
    
def test_primer_design_first_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([frags[0],x[1],x[2]], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]
    assert result.seq == (frags[0]+frags[1]+frags[2]).seq
    
def test_primer_design_second_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0],frags[1],x[2]], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]
    assert result.seq == (frags[0]+frags[1]+frags[2]).seq
    
def test_primer_design_third_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([x[0],x[1],frags[2]], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]
    assert result.seq == (frags[0]+frags[1]+frags[2]).seq
    
def test_primer_design_linker_first():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([ bam, x[0], x[1], x[2] ], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]
    assert result.seq == (bam + frags[0]+frags[1]+frags[2]).seq
    
def test_primer_design_linker_second():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([ x[0], bam, x[1], x[2] ], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == (frags[0] + bam + frags[1] + frags[2]).seq

def test_primer_design_linker_third():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([ x[0], x[1], bam, x[2] ], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == (frags[0] + frags[1] + bam + frags[2]).seq

def test_primer_design_linker_last():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([ x[0], x[1], x[2], bam], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == (frags[0] + frags[1] + frags[2] + bam).seq

def test_primer_design_linker_second_after_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([ frags[0], bam, x[1], x[2] ], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == (frags[0] + bam + frags[1] + frags[2]).seq
                    
def test_primer_design_linker_second_before_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([ x[0], bam, frags[1], x[2] ], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == (frags[0] + bam + frags[1] + frags[2]).seq
    
def test_primer_design_linker_third_after_Dseqrecord():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments([ x[0], frags[1], bam, x[2] ], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == (frags[0] + frags[1] + bam + frags[2]).seq
                          
def test_primer_design_two_fragments():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments( [ x[0], x[1] ], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == ( frags[0] + frags[1] ).seq

def test_primer_design_four_fragments():
    x = [primer_design(f) for f in frags]
    fourth = Dseqrecord("TAAAAATAAAATTGTTGACAGCAGAAGTGATATAGAAATTTGTTAATTATTA")
    y = assembly_fragments( x + [fourth], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == ( frags[0] + frags[1] + frags[2] + fourth).seq
                          
def test_primer_design_two_fragments_linker_in_between():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments( [ x[0], bam ,x[1] ], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == ( frags[0] + bam +  frags[1] ).seq

def test_primer_design_two_fragments_flanking_linkers():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments( [ bam, x[0], x[1], bam], 20)
    z = Assembly(y, limit=20)
    result = z.linear_products[0]   
    assert result.seq == ( bam + frags[0] + frags[1] + bam).seq
                          
def test_primer_design_one_fragment_flanking_linkers():
    x = [primer_design(f) for f in frags]
    y = assembly_fragments( [ bam, x[0], bam], 20)
    assert y[0].seq == ( bam + frags[0] + bam).seq
                          
                                  
                          

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])