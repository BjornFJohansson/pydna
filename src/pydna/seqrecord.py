#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""This module provide a subclass of the Biopython SeqRecord class.
It has a number of extra methods and uses
the :class:`pydna._pretty_str.pretty_str` class instread of str for a
nicer output in the IPython shell."""


from pydna.seqfeature import SeqFeature as _SeqFeature
from pydna._pretty import pretty_str as _pretty_str
from pydna.utils import seguid as _seg
from pydna.common_sub_strings import common_sub_strings as _common_sub_strings

# from Bio.Alphabet import generic_dna as _generic_dna
from Bio.Data.CodonTable import TranslationError as _TranslationError
from Bio.SeqUtils import GC as _GC
from Bio.SeqRecord import SeqRecord as _SeqRecord
from Bio.SeqFeature import FeatureLocation as _FeatureLocation
from Bio.Seq import Seq as _Seq
from prettytable import PrettyTable as _PrettyTable
import os as _os
import re as _re

from pydna import _PydnaWarning
from warnings import warn as _warn

import logging as _logging

_module_logger = _logging.getLogger("pydna." + __name__)


class SeqRecord(_SeqRecord):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.annotations.update({"molecule_type": "DNA"})
        if len(self.name) > 16:
            short_name = self.name[:16]
            _warn(
                "name property {} truncated to 16 chars {}".format(
                    self.name, short_name
                ),
                _PydnaWarning,
                stacklevel=2,
            )
            self.name = short_name

        if self.name == "<unknown name>":
            self.name = "name"

        if self.id == "<unknown id>":
            self.id = "id"

        if self.description == "<unknown description>":
            self.description = "description"

        self.map_target = None

        if not hasattr(self.seq, "transcribe"):
            self.seq = _Seq(self.seq)

        self.seq._data = "".join(self.seq._data.split())  # remove whitespaces
        # self.seq.alphabet = _generic_dna
        self.id = _pretty_str(self.id)
        self.name = _pretty_str(self.name)
        self.description = _pretty_str(self.description)
        self.annotations = {
            _pretty_str(k): _pretty_str(v) for k, v in self.annotations.items()
        }

    @property
    def locus(self):
        """ alias for name property """
        return self.name

    @locus.setter
    def locus(self, value):
        """ alias for name property """
        if len(value) > 16:
            shortvalue = value[:16]
            _warn(
                "locus property {} truncated to 16 chars {}".format(value, shortvalue),
                _PydnaWarning,
                stacklevel=2,
            )
            value = shortvalue
        self.name = value
        return

    @property
    def accession(self):
        """ alias for id property """
        return self.id

    @accession.setter
    def accession(self, value):
        """ alias for id property """
        self.id = value
        return

    @property
    def definition(self):
        """ alias for description property """
        return self.description

    @definition.setter
    def definition(self, value):
        """ alias for id property """
        self.description = value
        return

    def reverse_complement(self, *args, **kwargs):
        answer = super().reverse_complement(*args, **kwargs)
        answer.__class__ = type(self)
        # https://stackoverflow.com/questions/15404256/changing-the-class-of-a-python-object-casting
        # answer = type(self)(super().reverse_complement(*args,**kwargs).seq, *args,**kwargs)
        return answer

    def isorf(self, table=1):
        """Detects if sequence is an open reading frame (orf) in the 5'-3' direction.
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

        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("atgtaa")
        >>> a.isorf()
        True
        >>> b=SeqRecord("atgaaa")
        >>> b.isorf()
        False
        >>> c=SeqRecord("atttaa")
        >>> c.isorf()
        False

        References
        ----------

        .. [#] http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c

        """

        try:
            self.seq.translate(table=table, cds=True)
        except _TranslationError:
            return False
        else:
            return True

    def add_colors_to_features_for_ape(self):
        """This method assigns colors to features compatible with the
        `ApE editor <http://jorgensen.biology.utah.edu/wayned/ape/>`_"""

        cols = (
            "#66ffa3",
            "#84ff66",
            "#e0ff66",
            "#ffc166",
            "#ff6666",
            "#ff99d6",
            "#ea99ff",
            "#ad99ff",
            "#99c1ff",
            "#99ffff",
            "#99ffc1",
            "#adff99",
            "#eaff99",
            "#ffd699",
            "#ff9999",
            "#ffccea",
            "#f4ccff",
            "#d6ccff",
            "#cce0ff",
            "#ccffff",
            "#ccffe0",
            "#d6ffcc",
            "#f4ffcc",
            "#ffeacc",
            "#ffcccc",
            "#ff66c1",
            "#e066ff",
            "#8466ff",
            "#66a3ff",
            "#66ffff",
        )

        for i, f in enumerate(self.features):
            f.qualifiers["ApEinfo_fwdcolor"] = [cols[i % len(cols)]]
            f.qualifiers["ApEinfo_revcolor"] = [cols[::-1][i % len(cols)]]

    def add_feature(
        self, x=None, y=None, seq=None, type="misc", strand=1, *args, **kwargs
    ):

        #         location=None,
        #         type='',
        #         location_operator='',
        #         strand=None,
        #         id="<unknown id>",
        #         qualifiers=None,
        #         sub_features=None,
        #         ref=None,
        #         ref_db=None
        """Adds a feature of type misc to the feature list of the sequence.

        Parameters
        ----------
        x  : int
            Indicates start of the feature
        y  : int
            Indicates end of the feature

        Examples
        --------

        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("atgtaa")
        >>> a.features
        []
        >>> a.add_feature(2,4)
        >>> a.features
        [SeqFeature(FeatureLocation(ExactPosition(2), ExactPosition(4), strand=1), type='misc')]
        """
        qualifiers = {}
        qualifiers.update(kwargs)

        if seq:
            if hasattr(seq, "seq"):
                seq = seq.seq
                if hasattr(seq, "watson"):
                    seq = str(seq.watson).lower()
                else:
                    seq = str(seq).lower()
            else:
                seq = str(seq).lower()
            x = self.seq.lower().find(seq)
            if x == -1:
                raise TypeError("Could not find {}".format(seq))
            y = x + len(seq)
        else:
            x = x or 0
            y = y or len(self)

        if "label" not in qualifiers:
            qualifiers["label"] = ["ft{}".format(y - x)]

            if self[x:y].isorf() or self[x:y].reverse_complement().isorf():
                qualifiers["label"] = ["orf{}".format(y - x)]

        sf = _SeqFeature(
            _FeatureLocation(x, y, strand=strand), type=type, qualifiers=qualifiers
        )

        self.features.append(sf)

        """

        In [11]: a.seq.translate()
        Out[11]: Seq('K', ExtendedIUPACProtein())

        In [12]:
        """

    def list_features(self):
        """Prints an ASCII table with all features.

        Examples
        --------
        >>> from Bio.Seq import Seq
        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord(Seq("atgtaa"))
        >>> a.add_feature(2,4)
        >>> print(a.list_features())
        +-----+---------------+-----+-----+-----+-----+------+------+
        | Ft# | Label or Note | Dir | Sta | End | Len | type | orf? |
        +-----+---------------+-----+-----+-----+-----+------+------+
        |   0 | L:ft2         | --> | 2   | 4   |   2 | misc |  no  |
        +-----+---------------+-----+-----+-----+-----+------+------+"""

        x = _PrettyTable(
            ["Ft#", "Label or Note", "Dir", "Sta", "End", "Len", "type", "orf?"]
        )
        x.align["Ft#"] = "r"  # Left align
        x.align["Label or Note"] = "l"  # Left align
        x.align["Len"] = "r"
        x.align["Sta"] = "l"
        x.align["End"] = "l"
        x.align["type"] = "l"
        x.padding_width = 1  # One space between column edges and contents
        for i, sf in enumerate(self.features):
            try:
                lbl = sf.qualifiers["label"]
            except KeyError:
                try:
                    lbl = sf.qualifiers["note"]
                except KeyError:
                    lbl = "nd"
                else:
                    lbl = "N:{}".format(" ".join(lbl).strip())
            else:
                lbl = "L:{}".format(" ".join(lbl).strip())
            x.add_row(
                [
                    i,
                    lbl[:16],
                    {1: "-->", -1: "<--", 0: "---", None: "---"}[sf.strand],
                    sf.location.start,
                    sf.location.end,
                    len(sf),
                    sf.type,
                    {True: "yes", False: "no"}[
                        self.extract_feature(i).isorf()
                        or self.extract_feature(i).reverse_complement().isorf()
                    ],
                ]
            )
        return _pretty_str(x)

    def extract_feature(self, n):
        """Extracts a feature and returns a new SeqRecord object.

        Parameters
        ----------
        n  : int
            Indicates the feature to extract

        Examples
        --------

        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("atgtaa")
        >>> a.add_feature(2,4)
        >>> b=a.extract_feature(0)
        >>> b
        SeqRecord(seq=Seq('gt'), id='ft2', name='ft2', description='description', dbxrefs=[])
        """
        return self.features[n].extract(self)

    def sorted_features(self):
        """Returns a list of the features sorted by start position.

        Examples
        --------

        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("atgtaa")
        >>> a.add_feature(3,4)
        >>> a.add_feature(2,4)
        >>> print(a.features)
        [SeqFeature(FeatureLocation(ExactPosition(3), ExactPosition(4), strand=1), type='misc'), SeqFeature(FeatureLocation(ExactPosition(2), ExactPosition(4), strand=1), type='misc')]
        >>> b
        print(a.sorted_features())
        [SeqFeature(FeatureLocation(ExactPosition(2), ExactPosition(4), strand=1), type='misc'), SeqFeature(FeatureLocation(ExactPosition(3), ExactPosition(4), strand=1), type='misc')]
        """
        return sorted(self.features, key=lambda x: x.location.start)

    def stamp(self):
        """Adds a SEGUID or cSEGUID checksum to the description property.
        This will show in the genbank format.

        For linear sequences:

        ``SEGUID_<seguid>``

        For circular sequences:

        ``cSEGUID_<seguid>``



        Examples
        --------

        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("aaa")
        >>> a.stamp()
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE...'
        >>> a.description
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE...'


        """

        try:
            blunt = self.seq.isblunt()
        except AttributeError:
            blunt = True

        try:
            linear = self.seq.linear
        except AttributeError:
            linear = True

        if (not blunt) and linear:
            return _pretty_str("Sequence is not blunt nor circular,"
                               " so it can not be stamped.")

        algorithm = {True: "SEGUID", False: "cSEGUID"}[linear]
        chksum = getattr(self, algorithm.lower())()
        pattern = r"((?P<algorithm>c?SEGUID)(?:_|\s){1,5}(?P<sha1>\S{27}))"
        oldstamp = _re.search(pattern, self.description)

        if oldstamp:
            old_stamp, old_algorithm, old_chksum = oldstamp.groups()
            newstamp = _pretty_str("{}_{}".format(algorithm, chksum))
            if chksum == old_chksum and algorithm == old_algorithm:
                return newstamp
            else:
                raise ValueError("Stamp is wrong.\n"
                                 f"Old: {old_stamp}\n"
                                 f"New: {newstamp}")
        else:
            newstamp = "{}_{}".format(algorithm, chksum)
            if not self.description or self.description == "description":
                self.description = newstamp
            else:
                self.description += " " + newstamp
        return _pretty_str("{}_{}".format(algorithm, chksum))

    def seguid(self):
        """Returns the url safe SEGUID [#]_ for the sequence.
        This checksum is the same as seguid but with base64.urlsafe
        encoding [#]_ instead of the normal base 64. This means that
        the characters + and / are replaced with - and _ so that
        the checksum can be a part of and URL or a filename.

        Examples
        --------
        >>> from pydna.seqrecord import SeqRecord
        >>> a=SeqRecord("aaaaaaa")
        >>> a.seguid() # original seguid is +bKGnebMkia5kNg/gF7IORXMnIU
        '-bKGnebMkia5kNg_gF7IORXMnIU'

        References
        ----------

        .. [#] http://wiki.christophchamp.com/index.php/SEGUID"""
        return _seg(str(self.seq))


    def lcs(self, other, *args, limit=25, **kwargs):
        """Returns the longest common substring between the sequence
        and another sequence (other). The other sequence can be a string,
        Seq, SeqRecord, Dseq or DseqRecord.

        The method returns a SeqFeature with type "read" as this method
        is mostly used to map sequence reads to the sequence. This can be
        changed by passing a type as keyword with some other string value.

        Examples
        --------
        >>> from pydna.seqrecord import SeqRecord
        >>> a = SeqRecord("GGATCC")
        >>> a.lcs("GGATCC", limit=6)
        SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(6), strand=1), type='read')
        >>> a.lcs("GATC", limit=4)
        SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(6), strand=1), type='read')
        >>> a = SeqRecord("CCCCC")
        >>> a.lcs("GGATCC", limit=6)
        SeqFeature(None)

        """

        # longest_common_substring
        # https://biopython.org/wiki/ABI_traces
        if hasattr(other, "seq"):
            r = other.seq
            if hasattr(r, "watson"):
                r = str(r.watson).lower()
            else:
                r = str(r).lower()
        else:
            r = str(other.lower())

        olaps = _common_sub_strings(str(self.seq).lower(), r, limit, **kwargs)

        try:
            start_in_self, start_in_other, length = olaps.pop(0)
        except IndexError:
            result = _SeqFeature()
        else:
            label = "sequence" if not hasattr(other, "name") else other.name
            result = _SeqFeature(_FeatureLocation(start_in_self,
                                                  start_in_self+length),
                                  type = kwargs.get("type") or "read",
                                  strand=1,
                                  qualifiers={
                                  "label": [label],
                                  "ApEinfo_fwdcolor": ["#DAFFCF"],
                                  "ApEinfo_revcolor": ["#DFFDFF"],
                                  },)
        return result


    def gc(self):
        """Returns GC content"""
        return round(_GC(str(self.seq)), 1)

    def __lt__(self, other):
        try:
            return str(self.seq) < str(other.seq)
        except AttributeError:
            # I don't know how to compare to other
            return NotImplemented

    def __gt__(self, other):
        try:
            return str(self.seq) > str(other.seq)
        except AttributeError:
            # I don't know how to compare to other
            return NotImplemented

    def __eq__(self, other):
        try:
            if self.seq == other.seq and str(self.__dict__) == str(other.__dict__):
                return True
        except AttributeError:
            pass
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        """__hash__ must be based on __eq__"""
        return hash((str(self.seq).lower(), str(tuple(sorted(self.__dict__.items())))))

    def __str__(self):
        return _pretty_str(super().__str__())

    def __repr__(self):
        return _pretty_str(super().__repr__())

    def __format__(self, format):
        return _pretty_str(super().__format__(format))

    def __add__(self, other):
        answer = super().__add__(other)
        if answer.name == "<unknown name>":
            answer.name = "name"
        return answer

    def __getitem__(self, index):
        from pydna.utils import (
            identifier_from_string as _identifier_from_string,
        )  # TODO: clean this up

        answer = super().__getitem__(index)
        if len(answer) < 2:
            return answer
        identifier = "part_{id}".format(id=self.id)
        if answer.features:
            sf = max(answer.features, key=len)  # default
            if "label" in sf.qualifiers:
                identifier = " ".join(sf.qualifiers["label"])
            elif "note" in sf.qualifiers:
                identifier = " ".join(sf.qualifiers["note"])
        answer.id = _identifier_from_string(identifier)[:16]
        answer.name = _identifier_from_string("part_{name}".format(name=self.name))[:16]
        return answer


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
