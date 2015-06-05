#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test YEp24PGK_XK
'''

import unittest

class test_empty(unittest.TestCase):

    def test_empty(self):
        ''' test YEp24PGK_XK'''

        import os
        import pydna
        
        cwd = os.getcwd()
        YEp24PGK_XK_correct = pydna.read("YEp24PGK_XK_correct.gb")
        os.chdir("../docs/cookbook/")        
        p1 =   pydna.read("primer1.txt", ds = False)
        p3 =   pydna.read("primer3.txt", ds = False)
        XKS1 = pydna.read("XKS1_orf.txt")
        YEp24PGK = pydna.read("YEp24PGK.txt")
        
        os.chdir(cwd)

        PCR_prod = pydna.pcr(p1, p3, XKS1)

        from Bio.Restriction import BamHI
        
        stuffer1, insert, stuffer2 = PCR_prod.cut(BamHI)

        from Bio.Restriction import BglII

        YEp24PGK_BglII = YEp24PGK.cut(BglII).pop()

        YEp24PGK_XK = YEp24PGK_BglII + insert

        YEp24PGK_XK=YEp24PGK_XK.looped()

        YEp24PGK_XK = YEp24PGK_XK.synced("gaattctgaaccagtcctaaaacgagtaaataggaccggcaattc") #YEp24PGK)
        
        self.assertTrue( pydna.eq(YEp24PGK_XK, YEp24PGK_XK_correct))
        self.assertEqual( YEp24PGK_XK_correct.seguid() ,"HRVpCEKWcFsKhw_W-25ednUfldI" )
        self.assertEqual( YEp24PGK_XK.seguid() ,"HRVpCEKWcFsKhw_W-25ednUfldI" )

if __name__ == '__main__':
    unittest.main()









