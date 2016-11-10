#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import nose
import pydna

def test_pydna_download():

    files = [   ("sequence.gb", 'ZoC2nuFNWNTkogGMG19xAziKIn8'),
                ("NCBI_example.gb", 'ZoC2nuFNWNTkogGMG19xAziKIn8'),
                ("hej.txt", 'PYhRds7e1hBcIA4xbVICtYquFBY'),
                ("YEplac181.txt", '22AjVUSbrYu6RbTiG96Ruq6s08I'),
                ("pGADT7-Rec.gb", 'st4DzsZyh-_ZYDPSt0farDJRZY4'),
                ("P30350%20(2013-10-11%2013_49_14).dna.txt", 'JpBmJCXMcxOWtW7CK8Pema5c_OA'),
                ("ApE_example.gb", 'xn0tBiibm4SByUC29Byc7QiKOeQ'),
                ("VectorNTI_example.gb", 'mqcIyInPI3ZbwBsUdEUN8zQ7zaU'), ]

    for file_, seg in files:        
        #print('("'+file_+'", ', end="")        
        with open("tests/broken_genbank_files/"+file_, "r") as f:
            infile = f.read()
        assert pydna.Dseqrecord( pydna.gbtext_clean(infile).gbtext ).seguid() == seg
        #print(repr(pydna.Dseqrecord( pydna.gbtext_clean(infile).gbtext ).seguid())+"),")
        #print(gbclean_text(infile).jseq)

if __name__ == '__main__':
    nose.runmodule(argv=[sys.argv[0], '--nocapture'])
    

    
