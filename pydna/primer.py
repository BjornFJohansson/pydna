#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2013, 2014 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''
This module contain functions for primer design.

'''
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .tm import tmbresluc

class Primer(SeqRecord):
    '''This class holds information about a primer and its position on a template '''
    def __init__(self, record, *args,
                 position=None, footprint=None, tail=None, concentration = 1000.0, **kwargs):
        if hasattr(record, "features"):            
            for key, value in list(record.__dict__.items()):
                setattr(self, key, value )
        elif hasattr(record, "alphabet"):
            super().__init__(record, *args, **kwargs)            
        else:        
            super().__init__(Seq(record, IUPACAmbiguousDNA), *args, **kwargs)
        self.position  = position
        self.footprint = footprint
        self.tail      = tail
        self.concentration = concentration
    def __repr__(self):        
        return "{id} {len}-mer:5'{seq}-3'".format(id=self.id,len=len(self),seq=self.seq)
    def __radd__(self, other):
        new = super().__radd__(other)
        return Primer(new.seq)
    def tm(self, saltc=50.0, formula=tmbresluc):
        return formula( str(self.seq).upper(), primerc=self.concentration, saltc=saltc )
