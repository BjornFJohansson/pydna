#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test parse
'''

import pytest

from pydna.amplify       import pcr
from pydna.assembly      import Assembly
from pydna.dseqrecord    import Dseqrecord
from pydna.design import assembly_primers
from pydna.design import cloning_primers
from pydna.parsers       import parse


def test_primer_Design():
    ''' test_primer_design'''

    a=Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    b=Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc")
    c=Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt")

    primer_pairs = assembly_primers([a,b,c])

    frags=[]

    for (f,r),t in zip(primer_pairs,[a,b,c]):
        frags.append(pcr(f,r,t))

    asm=Assembly(frags)

    assert asm.linear_products[0].seguid() == "1eNv3d_1PqDPP8qJZIVoA45Www8"

    frags=[]

    primer_pairs = assembly_primers([a,b,c], circular=True)

    for (f,r),t in zip(primer_pairs,[a,b,c]):
        frags.append(pcr(f,r,t))

    #print frags

    asm=Assembly(frags)

    assert asm.circular_products[0].cseguid() == "V3Mi8zilejgyoH833UbjJOtDMbc"



def test_primer_Design_with_linker():
    ''' test_primer_design'''

    b  = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    l  = Dseqrecord("AAATTTCCCGGG")
    c  = Dseqrecord("tctgatctactatctgactgtactgatctattgacactgtgatcattctagtgtattactc")

    ((bf,br),(cf,cr)) = assembly_primers((b,l,c))

    nb = pcr((bf,br),b)
    nc = pcr((cf,cr),c)

    asm1 = Assembly((nb,nc))

    assert asm1.linear_products[0].seguid(),(b+l+c).seguid() == 'l95igKB8iKAKrvvqE9CYksyNx40'


    b  = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    l  = Dseqrecord("AAATTTCCCGGG")
    c  = Dseqrecord("tctgatctactatctgactgtactgatctattgacactgtgatcattctagtgtattactc")

    ((bf,br),(cf,cr)) = assembly_primers((b,l,c), circular = True)

    nb = pcr((bf,br),b)
    nc = pcr((cf,cr),c)

    asm = Assembly((nb,nc))

    #print (b+l+c).looped().seq

    assert (b+l+c).looped().cseguid() == asm.circular_products[0].cseguid()
    #print (b+l+c).looped().cseguid() == 'jdHXfQI5k4Sk2ESiZYfKv4oP2FI'

    assert (b+l+c).looped().cseguid() == 'jdHXfQI5k4Sk2ESiZYfKv4oP2FI'

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
