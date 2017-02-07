#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2013 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

from Bio.SeqRecord import SeqRecord as _SeqRecord
import datetime as _datetime
import pickle as _pickle
import shelve as _shelve

import copy     as _copy
import datetime as _datetime
import operator as _operator
import os       as _os
import re       as _re
import colorsys as _colorsys

from warnings import warn as _warn

import logging as _logging
_module_logger = _logging.getLogger("pydna."+__name__)

from prettytable import PrettyTable as _PrettyTable

from Bio.Seq                import Seq as _Seq
from Bio.Seq                import translate as _translate

from Bio.SeqRecord          import SeqRecord as _SeqRecord
from Bio.SeqFeature         import SeqFeature as _SeqFeature
from Bio.SeqFeature         import FeatureLocation as _FeatureLocation
from Bio.SeqFeature         import CompoundLocation as _CompoundLocation
from Bio.SeqUtils           import GC as _GC
from Bio.Data.CodonTable    import TranslationError as _TranslationError

from pydna._sequencetrace         import SequenceTraceFactory as _SequenceTraceFactory
from pydna.common_sub_strings import common_sub_strings as _common_sub_strings
from pydna.utils  import seguid  as _seg
from pydna.utils  import cseguid as _cseg
from pydna._pretty import pretty_str as _pretty_str
from pydna.dseq import Dseq as _Dseq
from pydna.utils import rc as _rc


class SeqRecord(_SeqRecord):

    def __init__(self,*args, **kwargs):
        super().__init__(*args, **kwargs)

        if len(self.name)>16:
            short_name = self.name[:16]
            _warn("name property {} truncated to 16 chars {}".format(self.name, short_name))
            self.name = short_name

        if self.name == "<unknown name>":
            self.name = "name?"

        if self.id == "<unknown id>":
            self.id = "id?"

        if self.description =="<unknown description>":
            self.description = "description?"

        if not 'date' in self.annotations:
            self.annotations.update({"date": _datetime.date.today().strftime("%d-%b-%Y").upper()})

        self.map_target = None

    @property
    def locus(self):
        ''' alias for name property '''
        return self.name

    @locus.setter
    def locus(self, value):
        ''' alias for name property '''
        if len(value)>16:
            raise Exception()
        self.name = value
        return

    @property
    def accession(self):
        ''' alias for id property '''
        return self.id

    @accession.setter
    def accession(self, value):
        ''' alias for id property '''
        self.id = value
        return

    @property
    def definition(self):
        ''' alias for description property '''
        return self.description

    @definition.setter
    def definition(self, value):
        ''' alias for id property '''
        self.description = value
        return

    def isorf(self, table=1):
        '''Detects if sequence is an open reading frame (orf) in the 5'-3' direction.
        Translation tables are numbers according to the NCBI numbering [#]_.

        Parameters
        ----------
        table  : int
            Sets the translation table, default is 1 (standard code)

        Returns
        -------
        bool
            True if sequence is an orf, False otherwise.


        Examples
        --------

        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("atgtaa")
        >>> a.isorf()
        True
        >>> b=Dseqrecord("atgaaa")
        >>> b.isorf()
        False
        >>> c=Dseqrecord("atttaa")
        >>> c.isorf()
        False

        References
        ----------

        .. [#] http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c

        '''

        try:
            self.seq.translate(table=table, cds=True)
        except _TranslationError:
            return False
        else:
            return True

    def spread_ape_colors(self):
        ''' This method assigns random colors compatible with the ApE editor
        to features.
        '''

        def get_N_HexCol(N):
            HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in range(N)]
            hex_out = []
            for rgb in HSV_tuples:
                rgb = [int(x*255) for x in _colorsys.hsv_to_rgb(*rgb)]
                hex_out.append("".join([chr(x).encode('hex') for x in rgb]))
            return hex_out

        for i, color in enumerate(get_N_HexCol(len(self.features))):
            self.features[i].qualifiers['ApEinfo_fwdcolor'] = "#"+color
            self.features[i].qualifiers['ApEinfo_revcolor'] = "#"+color

            
    def list_features(self):
        '''Prints an ASCII table with all features.

        Examples
        --------

        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("atgtaa")
        >>> a.add_feature(2,4)
        >>> print(a.list_features())
        +----------+-----------+-------+-----+--------+--------------+------+------+
        | Feature# | Direction | Start | End | Length | id           | type | orf? |
        +----------+-----------+-------+-----+--------+--------------+------+------+
        | 0        |    None   |   2   |  4  |      2 | <unknown id> | misc |  no  |
        +----------+-----------+-------+-----+--------+--------------+------+------+
        >>>
        '''

        x = _PrettyTable(["Feature#", "Direction", "Start", "End", "Length", "id", "type", "orf?"])
        x.align["Feature#"] = "l" # Left align
        x.align["Length"] = "r"
        x.align["id"] = "l"
        x.align["type"] = "l"
        x.padding_width = 1 # One space between column edges and contents
        for i, sf in enumerate(self.features):
            x.add_row([ i,
                        {1:"-->", -1:"<--", 0:"---", None:"None"}[sf.strand],
                        sf.location.start,
                        sf.location.end,
                        len(sf), sf.id, sf.type,
                        {True:"yes",False:"no"}[self.extract_feature(i).isorf() or self.extract_feature(i).rc().isorf()]])
        return _pretty_str(x)

    def gc(self):
        '''Returns GC content '''
        return round(_GC(str(self.seq)), 1)

            
        
        
        
if __name__=="__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
