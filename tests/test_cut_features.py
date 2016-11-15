#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test cut
'''
import pytest
import sys

from pydna import pcr, cloning_primers, read, Dseqrecord

from Bio.Restriction import EcoRI

def test_cut_feat():
    puc19 = read('PUC19_MarkBudde.gb')
    pf, pr = cloning_primers(puc19)
    pcrProd = pcr(pf, pr, puc19)
    assert len(pcrProd.features) == 23
    #print len(pcrProd.cut(EcoRI)[1].features)
    assert len(pcrProd.cut(EcoRI)[1].features) == 17

    def amplicon_to_dseqrecord(a):
        d = Dseqrecord(a.seq)
        d.features = a.features
        return d

    pcrProdDseqrecord = amplicon_to_dseqrecord(pcrProd)
    assert len(pcrProdDseqrecord.cut(EcoRI)[1].features) == 17

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
