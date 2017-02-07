#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pytest
from pydna.readers      import read
from pydna.genbankfixer import gbtext_clean

def test_pydna_gbtext_clean():

    files = [    ("sequence.gb", "j2yAlBCZ-txSTCkakAmykAielRI"),
                 ("NCBI_example.gb", "j2yAlBCZ-txSTCkakAmykAielRI"),
                 ("YEplac181.txt", "lbnQtxi5LyDONRswRdG88-l8NF0"),
                 ("pGADT7-Rec.gb", "rhnYE78wGKdqAZWiyVJzQ7HXXys"),
                 ("P30350%20(2013-10-11%2013_49_14).dna.txt", "_aEPoGLctHcOZdQdZIh-KyBt5WY"),
                 ("ApE_example.gb", "c47i2ifiNZVvvnLQbX5anTVVoPE"),
                 ("VectorNTI_example.gb", "bDPbx5P4yigGWh1zK7FiG_SF8qQ"),
                 ("hej.txt", "lbnQtxi5LyDONRswRdG88-l8NF0"), ]

    for file_, seg in files:        
        #print('("'+file_+'", ', end="")        
        with open("broken_genbank_files/"+file_, "r") as f:
            infile = f.read()
        #print(file_, pydna.read( gbtext_clean(infile).gbtext ).seguid(), seg)
        assert read( gbtext_clean(infile).gbtext ).seguid() == seg
        #calcseg = pydna.read( pydna.gbtext_clean(infile).gbtext ).seguid()
        #print('"'+calcseg+'"),')

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
