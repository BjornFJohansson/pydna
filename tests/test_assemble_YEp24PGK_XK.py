#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test YEp24PGK_XK
'''
import pytest
import sys

from pydna.readers  import read
from pydna.utils    import eq
from pydna.amplify  import pcr

def test_assemble_YEp24PGK_XK():
    ''' test YEp24PGK_XK'''       
    p1 =   read("primer1.txt", ds = False)
    p3 =   read("primer3.txt", ds = False)
    XKS1 = read("XKS1_orf.txt")
    YEp24PGK = read("YEp24PGK.txt")

    PCR_prod = pcr(p1, p3, XKS1)

    from Bio.Restriction import BamHI
    
    stuffer1, insert, stuffer2 = PCR_prod.cut(BamHI)

    from Bio.Restriction import BglII

    YEp24PGK_BglII = YEp24PGK.cut(BglII).pop()

    YEp24PGK_XK = YEp24PGK_BglII + insert

    YEp24PGK_XK=YEp24PGK_XK.looped()

    YEp24PGK_XK = YEp24PGK_XK.synced("gaattctgaaccagtcctaaaacgagtaaataggaccggcaattc") #YEp24PGK)

    assert YEp24PGK_XK.seguid()  == "HRVpCEKWcFsKhw_W-25ednUfldI"
    assert YEp24PGK_XK.cseguid() == "t9fs_9UvEuD-Ankyy8XEr1hD5DQ"

    YEp24PGK_XK_correct = read("YEp24PGK_XK_manually_assembled.txt")
    assert YEp24PGK_XK_correct.cseguid() == "t9fs_9UvEuD-Ankyy8XEr1hD5DQ"
    assert eq(YEp24PGK_XK, YEp24PGK_XK_correct)

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])










