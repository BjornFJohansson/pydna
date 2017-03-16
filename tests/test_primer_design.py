#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test parse
'''

import pytest

from pydna.amplify       import pcr
from pydna.assembly      import Assembly
from pydna.dseqrecord    import Dseqrecord
from pydna.design        import primer_design
from pydna.design        import assembly_fragments


def test_primer_Design():
    ''' test_primer_design'''

    a=Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    b=Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc")
    c=Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt")

    frags = assembly_fragments( [primer_design(r) for r in (a,b,c)] )

    asm=Assembly(frags)

    assert asm.linear_products[0].seguid() == "1eNv3d_1PqDPP8qJZIVoA45Www8"

    frags = assembly_fragments( [primer_design(r) for r in (a,b,c,a)] )
    
    a2 = pcr(frags[-1].forward_primer, frags[0].reverse_primer, a)

    asm=Assembly( (a2, frags[1], frags[2]) )

    assert asm.circular_products[0].cseguid() == "V3Mi8zilejgyoH833UbjJOtDMbc"

def test_primer_Design_with_linker():
    ''' test_primer_design'''

    b  = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    l  = Dseqrecord("AAATTTCCCGGG")
    c  = Dseqrecord("tctgatctactatctgactgtactgatctattgacactgtgatcattctagtgtattactc")

    frags = assembly_fragments( (primer_design(b),l,primer_design(c)) )

    asm1 = Assembly(frags)

    assert asm1.linear_products[0].seguid(),(b+l+c).seguid() == 'l95igKB8iKAKrvvqE9CYksyNx40'

    frags = assembly_fragments( (primer_design(b),l,primer_design(c), primer_design(b)) )
    
    b2 = pcr(frags[-1].forward_primer, frags[0].reverse_primer, b)
    
    asm2 = Assembly( (b2, frags[1], frags[2]) )

    assert (b+l+c).looped().cseguid() == asm2.circular_products[0].cseguid()

    assert (b+l+c).looped().cseguid() == 'jdHXfQI5k4Sk2ESiZYfKv4oP2FI'

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
