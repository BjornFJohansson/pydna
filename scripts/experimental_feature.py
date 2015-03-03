#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''Experimental'''

from Bio.SeqFeature import SeqFeature

class Feature(SeqFeature):

    def __init__(self, seq="", *args,  **kwargs):
        SeqFeature.__init__(self, location = None, type = '', location_operator = '',
                 strand = None, id = "<unknown id>",
                 qualifiers = None, sub_features = None,
                 ref = None, ref_db = None, *args, **kwargs)
        self.seq = seq

if __name__=="__main__":
    a = Feature(seq = "x")





