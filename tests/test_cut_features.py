#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test cut
'''
import unittest

from pydna import pcr, cloning_primers, read, Dseqrecord

from Bio.Restriction import EcoRI

class test_cut_features(unittest.TestCase):

    def test_cut_feat(self):
        puc19 = read('PUC19_MarkBudde.gb')
        pf, pr = cloning_primers(puc19)
        pcrProd = pcr(pf, pr, puc19)
        self.assertEqual(23, len(pcrProd.features))
        #print len(pcrProd.cut(EcoRI)[1].features)
        self.assertEqual(17, len(pcrProd.cut(EcoRI)[1].features))

        def amplicon_to_dseqrecord(a):
            d = Dseqrecord(a.seq)
            d.features = a.features
            return d

        pcrProdDseqrecord = amplicon_to_dseqrecord(pcrProd)
        self.assertEqual(17, len(pcrProdDseqrecord.cut(EcoRI)[1].features))

if __name__ == '__main__':
    unittest.main()
