#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2013 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

from .dseqrecord            import Dseqrecord

class GenbankRecord(Dseqrecord):

    def __init__(self, record, *args, item = "", start=None, stop=None, strand=1,**kwargs):
        super().__init__(record, *args, **kwargs)
        self.item = item
        self.start = start
        self.stop = stop
        self.strand = strand

    def __repr__(self):
        '''returns a short string representation of the object'''
        return "Genbank({})({}{})".format(self.id, {True:"-", False:"o"}[self.linear],len(self))
        
    def _repr_pretty_(self, p, cycle):
        '''returns a short string representation of the object'''
        p.text("Genbank({})({}{})".format(self.id, {True:"-", False:"o"}[self.linear],len(self)))
            
    def _repr_html_(self):
        if not self.item:
            return self.__repr__()
        linktext = self.item
        if self.start != None and self.stop != None:
            linktext += " {}-{}".format(self.start, self.stop)
        return "<a href='https://www.ncbi.nlm.nih.gov/nuccore/{item}?from={start}&to={stop}&strand={strand}' target='_blank'>{linktext}</a>".format(item=self.item, start=self.start or "", stop=self.stop or "", strand=self.strand, linktext=linktext)
        
