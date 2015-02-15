#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test parse
'''

import unittest
from pydna import pcr, Assembly, Dseqrecord, assembly_primers

class test_primer_design(unittest.TestCase):

    def test_primer_Design(self):
        ''' test_primer_design'''



        a=Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
        b=Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc")
        c=Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt")

        primer_pairs = assembly_primers([a,b,c])

        frags=[]

        for (f,r),t in zip(primer_pairs,[a,b,c]):
            frags.append(pcr(f,r,t))

        asm=Assembly(frags)

        self.assertEqual(asm.linear_products[0].seguid(), "1eNv3d/1PqDPP8qJZIVoA45Www8")

        frags=[]

        primer_pairs = assembly_primers([a,b,c], circular=True)

        for (f,r),t in zip(primer_pairs,[a,b,c]):
            frags.append(pcr(f,r,t))

        print frags

        asm=Assembly(frags)

        self.assertEqual(asm.circular_products[0].cseguid(), "V3Mi8zilejgyoH833UbjJOtDMbc")



    def test_primer_Design_with_linker(self):
        ''' test_primer_design'''


        from pydna import Dseqrecord, Assembly, pcr, assembly_primers

        b  = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
        l  = Dseqrecord("AAATTTCCCGGG")
        c  = Dseqrecord("tctgatctactatctgactgtactgatctattgacactgtgatcattctagtgtattactc")

        ((bf,br),(cf,cr)) = assembly_primers((b,l,c))

        nb = pcr((bf,br),b)
        nc = pcr((cf,cr),c)

        asm1 = Assembly((nb,nc))

        self.assertEqual(asm1.linear_products[0].seguid(),(b+l+c).seguid(),'l95igKB8iKAKrvvqE9CYksyNx40')


        b  = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
        l  = Dseqrecord("AAATTTCCCGGG")
        c  = Dseqrecord("tctgatctactatctgactgtactgatctattgacactgtgatcattctagtgtattactc")

        ((bf,br),(cf,cr)) = assembly_primers((b,l,c), circular = True)

        nb = pcr((bf,br),b)
        nc = pcr((cf,cr),c)

        asm = Assembly((nb,nc))

        self.assertEqual((b+l+c).looped().cseguid(), asm.circular_products[0].cseguid(), 'qMEHxKkTsWIXkbqGA5O35631eMU')


if __name__ == '__main__':
    unittest.main()









