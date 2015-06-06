#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test pGUP1
'''

import unittest

class test_empty(unittest.TestCase):

    def test_empty(self):
        ''' test pGUP1'''

        import os

        cwd = os.getcwd()

        os.chdir("../docs/cookbook/")
        
        import pydna

        GUP1rec1sens = pydna.read("GUP1rec1sens.txt")
        GUP1rec2AS = pydna.read("GUP1rec2AS.txt")
        GUP1_locus = pydna.read("GUP1_locus.gb")
        pGREG505 = pydna.read("pGREG505.gb")
        
        os.chdir(cwd)

        insert = pydna.pcr(GUP1rec1sens, GUP1rec2AS, GUP1_locus)

        from Bio.Restriction import SalI

        lin_vect, his3 = pGREG505.cut(SalI)

        a = pydna.Assembly([insert, lin_vect], limit=28)
        
        pGUP1 = a.circular_products[0]
       
        pGUP1 = pGUP1.synced(pGREG505.seq[:50])    
        
        pGUP1_correct = pydna.read("pGUP1_correct.gb")        
        
        self.assertEqual(len(pGUP1_correct), 9981)
        self.assertEqual(len(pGUP1), 9981)
        self.assertTrue( pydna.eq(pGUP1, pGUP1_correct) )        
        self.assertEqual(pGUP1_correct.seguid(), "42wIByERn2kSe_Exn405RYwhffU")        
        self.assertEqual(pGUP1.seguid(), "42wIByERn2kSe_Exn405RYwhffU")       

if __name__ == '__main__':
    unittest.main()









