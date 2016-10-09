#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test YEp24PGK_XK
'''

import nose, sys
import pydna

def test_empty():
    ''' test YEp24PGK_XK'''

    YEp24PGK_XK_correct = pydna.read("tests/YEp24PGK_XK_correct.gb")
       
    p1 =   pydna.read("docs/cookbook/primer1.txt", ds = False)
    p3 =   pydna.read("docs/cookbook/primer3.txt", ds = False)
    XKS1 = pydna.read("docs/cookbook/XKS1_orf.txt")
    YEp24PGK = pydna.read("docs/cookbook/YEp24PGK.txt")

    PCR_prod = pydna.pcr(p1, p3, XKS1)

    from Bio.Restriction import BamHI
    
    stuffer1, insert, stuffer2 = PCR_prod.cut(BamHI)

    from Bio.Restriction import BglII

    YEp24PGK_BglII = YEp24PGK.cut(BglII).pop()

    YEp24PGK_XK = YEp24PGK_BglII + insert

    YEp24PGK_XK=YEp24PGK_XK.looped()

    YEp24PGK_XK = YEp24PGK_XK.synced("gaattctgaaccagtcctaaaacgagtaaataggaccggcaattc") #YEp24PGK)
    
    assert pydna.eq(YEp24PGK_XK, YEp24PGK_XK_correct)
    assert YEp24PGK_XK_correct.seguid() == "HRVpCEKWcFsKhw_W-25ednUfldI"
    assert YEp24PGK_XK.seguid() == "HRVpCEKWcFsKhw_W-25ednUfldI"

if __name__ == '__main__':
    nose.runmodule(argv=[sys.argv[0], '--nocapture'])










