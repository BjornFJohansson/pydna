#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test pGUP1
'''
import pytest

from pydna.readers  import read
from pydna.assembly import Assembly
from pydna.utils    import eq
from pydna.amplify  import pcr

def test_assemble_pGUP1():
    ''' test pGUP1'''
   
    GUP1rec1sens = read("GUP1rec1sens.txt")
    GUP1rec2AS = read("GUP1rec2AS.txt")
    GUP1_locus = read("GUP1_locus.gb")
    pGREG505 = read("pGREG505.gb")

    insert = pcr(GUP1rec1sens, GUP1rec2AS, GUP1_locus)

    from Bio.Restriction import SalI

    lin_vect, his3 = pGREG505.cut(SalI)

    a = Assembly([insert, lin_vect], limit=28)
    
    pGUP1 = a.circular_products[0]
   
    pGUP1 = pGUP1.synced(pGREG505.seq[:50])    
    
    pGUP1_correct = read("pGUP1_correct.gb")        
    
    assert len(pGUP1_correct) == 9981
    assert len(pGUP1) == 9981
    assert eq(pGUP1, pGUP1_correct)     
    assert pGUP1_correct.seguid() == "42wIByERn2kSe_Exn405RYwhffU"       
    assert pGUP1.seguid() == "42wIByERn2kSe_Exn405RYwhffU"    

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
