#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2013 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''Provides Dseqrecord, for handling double stranded
DNA sequences. Dseq and Dseqrecord are subclasses of Biopythons
Seq and SeqRecord classes, respectively. These classes support the
notion of circular and linear DNA.

'''
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
from pydna.utils  import seguid   as _seg
from pydna.utils  import cseguid  as _cseg
from pydna.utils  import rc       as _rc
from pydna.utils  import memorize as _memorize

from pydna._pretty import pretty_str as _pretty_str
from pydna.dseq import Dseq as _Dseq


try:
    from IPython.display import display as _display
    from IPython.display import display_html as _display_html
    from IPython.display import HTML    as _HTML
except ImportError:
    def _display(item): return item
    def _HTML(item): return item

class Dseqrecord(_SeqRecord):
    '''Dseqrecord is a double stranded version of the Biopython SeqRecord [#]_ class.
    The Dseqrecord object holds a Dseq object describing the sequence.
    Additionally, Dseqrecord hold meta information about the sequence in the
    from of a list of SeqFeatures, in the same way as the SeqRecord does.
    The Dseqrecord can be initialized with a string, Seq, Dseq, SeqRecord
    or another Dseqrecord. The sequence information will be stored in a
    Dseq object in all cases. Dseqrecord objects can be read or parsed
    from sequences in Fasta, Embl or Genbank format.

    There is a short representation associated with the Dseqrecord.
    ``Dseqrecord(-3)`` represents a linear sequence of length 2
    while ``Dseqrecord(o7)``
    represents a circular sequence of length 7.

    Dseqrecord and Dseq share the same concept of length

    ::

        <-- length -->
        GATCCTTT
             AAAGCCTAG




    Parameters
    ----------
    record  : string, Seq, SeqRecord, Dseq or other Dseqrecord object
        This data will be used to form the seq property

    circular : bool, optional
        True or False reflecting the shape of the DNA molecule

    linear : bool, optional
        True or False reflecting the shape of the DNA molecule


    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> a=Dseqrecord("aaa")
    >>> a
    Dseqrecord(-3)
    >>> a.seq
    Dseq(-3)
    aaa
    ttt
    >>> from Bio.Seq import Seq
    >>> b=Dseqrecord(Seq("aaa"))
    >>> b
    Dseqrecord(-3)
    >>> b.seq
    Dseq(-3)
    aaa
    ttt
    >>> from Bio.SeqRecord import SeqRecord
    >>> c=Dseqrecord(SeqRecord(Seq("aaa")))
    >>> c
    Dseqrecord(-3)
    >>> c.seq
    Dseq(-3)
    aaa
    ttt
    >>> a.seq.alphabet
    IUPACAmbiguousDNA()
    >>> b.seq.alphabet
    IUPACAmbiguousDNA()
    >>> c.seq.alphabet
    IUPACAmbiguousDNA()
    >>>

    References
    ----------

    .. [#] http://biopython.org/wiki/SeqRecord

    '''

    def __init__(self, record,
                       *args,
                       circular               = None,
                       linear                 = None,
                       n                      = 5E-14, # mol ( = 0.05 pmol)
                       **kwargs):
        self.n = n
        if circular == None and linear in (True, False,):
            circular = not linear
        elif linear == None and circular in (True, False,):
            linear   = not circular

        # record is Dseq object ?
        if hasattr(record, "watson"): 
            if record.circular and linear:
                record = record.tolinear()
            if record.linear and circular:
                record = record.looped()
            super().__init__(record, *args, **kwargs)
        # record is a Bio.Seq object ?
        elif hasattr(record, "alphabet"):          
            super().__init__(_Dseq(str(record),
                                  str(record.reverse_complement()),
                             *args,
                             ovhg=0 ,
                             linear=linear,
                             circular=circular),                                          
                             **kwargs)
        # record is aBio.SeqRecord or Dseqrecord object ?           
        elif hasattr(record, "features"):            
            for key, value in list(record.__dict__.items()):
                setattr(self, key, value )
            record.letter_annotations = {}
            # record.seq is a Dseq object ?
            if hasattr(record.seq, "watson"):          
                new_seq = _copy.copy(record.seq)
                if new_seq.circular and linear:
                    new_seq = new_seq.tolinear()
                if new_seq.linear and circular:
                    new_seq = new_seq.looped()
                self.seq=new_seq
            # record is Bio.SeqRecord object ?
            else:                                   
                self.seq=_Dseq(str(self.seq),
                              _rc(str(self.seq)),
                              ovhg=0 ,
                              linear=linear,
                              circular=circular)
        # assume that record is a string
        else:                                   
            super().__init__(_Dseq(record,
                                  _rc(record),
                                  ovhg=0 ,
                                  linear=linear,
                                  circular=circular),
                             *args,
                             **kwargs)

        if len(self.name)>16:
            short_name = self.name[:16]            
            _module_logger.warning("name property %s truncated to 16 chars: %s", self.name, short_name)
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
        #self.key = str( self.__hash__() )
    
    @property
    def linear(self):
        '''The linear property'''
        return self.seq.linear

    @property
    def circular(self):
        '''The circular property'''
        return self.seq.circular

    @property
    def locus(self):
        ''' alias for name property '''
        return self.name

    @locus.setter
    def locus(self, value):
        ''' alias for name property, max 16 letters'''
        self.name = value[:16]
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

    def seguid(self):
        '''Returns the url safe SEGUID [#]_ for the sequence.
           This checksum is the same as seguid but with base64.urlsafe
           encoding [#]_ instead of the normal base 64. This means that
           the characters + and / are replaced with - and _ so that
           the checksum can be a pert of and URL or a filename.

           Examples
           --------
           >>> from pydna.dseqrecord import Dseqrecord
           >>> a=Dseqrecord("aaaaaaa")
           >>> a.seguid() # original seguid is +bKGnebMkia5kNg/gF7IORXMnIU
           '-bKGnebMkia5kNg_gF7IORXMnIU'

           References
           ----------

       .. [#] http://wiki.christophchamp.com/index.php/SEGUID

       '''
        return _seg(self.seq)

    def m(self):
        return self.seq.mw() * self.n # Da(g/mol) * mol = g

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

    def add_feature(self, x=None, y=None, seq=None, label=None, type="misc", **kwargs):
        '''Adds a feature of type misc to the feature list of the sequence.

        Parameters
        ----------
        x  : int
            Indicates start of the feature
        y  : int
            Indicates end of the feature

        Examples
        --------

        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("atgtaa")
        >>> a.features
        []
        >>> a.add_feature(2,4)
        >>> a.features
        [SeqFeature(FeatureLocation(ExactPosition(2), ExactPosition(4)), type='misc')]
        '''
        qualifiers = {"label": label}
        if seq:
            seq = Dseqrecord(seq)
            x = self.seq.lower().find(seq.seq.lower())
            if x==-1:
                return
            y = x + len(seq)
        self.features.append( _SeqFeature(_FeatureLocation(x, y), type=type, qualifiers = qualifiers, **kwargs))

        '''

        In [11]: a.seq.translate()
        Out[11]: Seq('K', ExtendedIUPACProtein())

        In [12]:
        '''


    def extract_feature(self, n):
        '''Extracts a feature and creates a new Dseqrecord object.

        Parameters
        ----------
        n  : int
            Indicates the feature to extract

        Examples
        --------

        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("atgtaa")
        >>> a.add_feature(2,4)
        >>> b=a.extract_feature(0)
        >>> b
        Dseqrecord(-2)
        >>> b.seq
        Dseq(-2)
        gt
        ca

        '''
        return self.features[n].extract(self)

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

    def olaps(self, other, *args, **kwargs):
        other = Dseqrecord(other)
        olaps = _common_sub_strings(str(self.seq).lower(), str(other.seq).lower(), **kwargs)
        return [ self[olap[0]:olap[0]+olap[2]] for olap in olaps ]

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

    def cseguid(self):
        '''Returns the url safe cSEGUID for the sequence.

        Only defined for circular sequences.

        The cSEGUID checksum uniquely identifies a circular
        sequence regardless of where the origin is set.
        The two Dseqrecord objects below are circular
        permutations.

        Examples
        --------

        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("agtatcgtacatg", circular=True)
        >>> a.cseguid() # cseguid is CTJbs6Fat8kLQxHj+/SC0kGEiYs
        'CTJbs6Fat8kLQxHj-_SC0kGEiYs'

        >>> a=Dseqrecord("gagtatcgtacat", circular=True)
        >>> a.cseguid()
        'CTJbs6Fat8kLQxHj-_SC0kGEiYs'

       '''
        if self.linear:
            raise Exception("cseguid is only defined for circular sequences.")
        return _cseg(self.seq)

    def lseguid(self):
        '''Returns the url safe lSEGUID for the sequence.

        Only defined for linear double stranded sequences.

        The lSEGUID checksum uniquely identifies a linear
        sequence independent of the direction.
        The two Dseqrecord objects below are each others
        reverse complements, so they do in fact refer to
        the same molecule.

        Examples
        --------

        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("agtatcgtacatg")
        >>> a.lseguid()
        'DPshMN4KeAjMovEjGEV4Kzj18lU'

        >>> b=Dseqrecord("catgtacgatact")
        >>> a.lseguid()
        'DPshMN4KeAjMovEjGEV4Kzj18lU'

       '''
        if self.circular:
            raise Exception("lseguid is only defined for linear sequences.")
        return self.seq.seguid()

    def stamp(self):
        '''Adds a SEGUID or cSEGUID checksum and a datestring to the description property.
        This will show in the genbank format. 
        
        For linear sequences:

        ``SEGUID_<seguid>_<datestring>``

        For circular sequences:

        ``cSEGUID_<seguid>_<datestring>``

        If there is already a stamp,  
        
        The stamp can be verified with :func:`verify_stamp`

        Examples
        --------

        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("aaa")
        >>> a.stamp()
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE...'
        >>> a.description
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE...'
        >>> a.verify_stamp()
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE'

        See also
        --------
        pydna.dseqrecord.Dseqrecord.verify_stamp
        '''
     
        alg = {True:"SEGUID", False:"cSEGUID"}[self.linear]        
        chksum = getattr(self, alg.lower())()        
        pattern = "(SEGUID|cSEGUID)_\s*(\S{27})_(\S{26})"
        oldstamp = _re.search(pattern, self.description)
        
        if oldstamp:
            old_alg, old_chksum, old_datestring = oldstamp.groups()
            if alg==old_alg and chksum==old_chksum:
                return _pretty_str("{}_{}".format(alg,chksum))
            else:
                raise Exception("Stamp incorrect.")
        else:
            datestring = _datetime.datetime.utcnow().isoformat("T")
            newstamp = "{}_{}_{}".format(alg, chksum, datestring)
            if not self.description or self.description=="description?":
                self.description = newstamp
            else:
                self.description += " "+newstamp
        return _pretty_str("{}_{}".format(alg, chksum))

    def verify_stamp(self):
        '''Verifies the SEGUID stamp in the description property is
       valid. True if stamp match the sequid calculated from the sequence.
       Exception raised if no stamp can be found.

        >>> from pydna.dseqrecord import Dseqrecord
        >>> from pydna.readers import read
        >>> b=read(">a\\naaa")
        >>> b.annotations['date'] = '02-FEB-2013'
        >>> b.seguid()
        'YG7G6b2Kj_KtFOX63j8mRHHoIlE'
        >>> print(b.format("gb"))
        LOCUS       a                          3 bp    DNA     linear   UNK 02-FEB-2013
        DEFINITION  a
        ACCESSION   a
        VERSION     a
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 aaa
        //
        >>> b.stamp()
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE'
        >>> b
        Dseqrecord(-3)
        >>> print(b.format("gb"))
        LOCUS       a                          3 bp    DNA     linear   UNK 02-FEB-2013
        DEFINITION  a SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE_...
        ACCESSION   a
        VERSION     a
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 aaa
        //
        >>> b.verify_stamp()
        'SEGUID_YG7G6b2Kj_KtFOX63j8mRHHoIlE'
        >>>

       See also
       --------
       pydna.dseqrecord.Dseqrecord.stamp

       '''
        name, alg = {True:("SEGUID", _seg), False:("cSEGUID", _cseg)}[self.linear]
        pattern = "{name}_{chksum}".format(name=name, chksum=alg(self.seq))

        if not pattern in self.description:
            raise Exception("No stamp present in the description property.")
        else:
            return _pretty_str(pattern)

    def looped(self):
        '''
        Returns a circular version of the Dseqrecord object. The
        underlying Dseq object has to have compatible ends.


        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("aaa")
        >>> a
        Dseqrecord(-3)
        >>> b=a.looped()
        >>> b
        Dseqrecord(o3)
        >>>

        See also
        --------
        pydna.dseq.Dseq.looped
        '''
        new = _copy.copy(self)
        for key, value in list(self.__dict__.items()):
            setattr(new, key, value )
        new._seq = self.seq.looped()
        for fn, fo in zip(new.features, self.features):
            fn.qualifiers = fo.qualifiers

        return new

    def tolinear(self):
        '''
        Returns a linear, blunt copy of a circular Dseqrecord object. The
        underlying Dseq object has to be circular.

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("aaa", circular = True)
        >>> a
        Dseqrecord(o3)
        >>> b=a.tolinear()
        >>> b
        Dseqrecord(-3)
        >>>

        '''

        new = _copy.copy(self)
        for key, value in list(self.__dict__.items()):
            setattr(new, key, value )
        new._seq = self.seq.tolinear()
        for fn, fo in zip(new.features, self.features):
            fn.qualifiers = fo.qualifiers

        return new

    def format(self, f="gb"):
        '''Returns the sequence as a string using a format supported by Biopython
        SeqIO [#]_. Default is "gb" which is short for Genbank.

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> x=Dseqrecord("aaa")
        >>> x.annotations['date'] = '02-FEB-2013'
        >>> x
        Dseqrecord(-3)
        >>> print(x.format("gb"))
        LOCUS       name?                      3 bp    DNA     linear   UNK 02-FEB-2013
        DEFINITION  description?
        ACCESSION   id?
        VERSION     id?
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 aaa
        //


        References
        ----------

        .. [#] http://biopython.org/wiki/SeqIO


        '''

        s = _SeqRecord.format(self, f).strip()

        if f in ("genbank","gb"):
            if self.circular:
                return _pretty_str(s[:55]+"circular"+s[63:])
            else:
                return _pretty_str(s[:55]+"linear"+s[61:])
        else:
            return _pretty_str(s)

    def write(self, filename=None, f="gb"):
        '''Writes the Dseqrecord to a file using the format f, which must
        be a format supported by Biopython SeqIO for writing [#]_. Default
        is "gb" which is short for Genbank. Note that Biopython SeqIO reads
        more formats than it writes.

        Filename is the path to the file where the sequece is to be
        written. The filename is optional, if it is not given, the
        description property (string) is used together with the format.

        If obj is the Dseqrecord object, the default file name will be:

        ``<obj.locus>.<f>``

        Where <f> is "gb" by default. If the filename already exists and
        AND the sequence it contains is different, a new file name will be
        used so that the old file is not lost:

        ``<obj.locus>_NEW.<f>``

        References
        ----------
        .. [#] http://biopython.org/wiki/SeqIO

        '''
        msg=""
        if not filename:
            filename = "{name}.{type}".format(name=self.locus, type=f)
            # generate a name if no name was given
        if str(filename)==filename:                 # is filename a string???
            name, ext = _os.path.splitext(filename)
            msg = "<font face=monospace><a href='{filename}' target='_blank'>{filename}</a></font><br>".format(filename=filename)
            if not _os.path.isfile(filename):
                with open(filename, "w") as fp: fp.write(self.format(f))
            else:
                from pydna.readers import read
                old_file = read(filename)
                if self.seq != old_file.seq:
                    # If new sequence is different, the old file is renamed with "OLD" suffix:
                    old_filename = "{}_OLD{}".format(name, ext)
                    _os.rename(filename, old_filename)
                    
                    msg = ("<font color='DarkOrange ' face=monospace>"
                           "Sequence changed.<br>"
                           "</font>"
                           "<font color='red' face=monospace>"
                           "new: <a href='{filename}' target='_blank'>{filename}</a> &nbsp&nbsp&nbsp size: {nlen}bp topology: {ntop} SEGUID: {ns}<br>"
                           "</font>"
                           "<font color='green' face=monospace>"
                           "old: <a href='{oldfname}' target='_blank'>{oldfname}</a> size: {olen}bp topology: {otop} SEGUID: {os}<br>"
                           "</font>").format(filename=filename, 
                                             oldfname=old_filename, 
                                             nlen=len(self),
                                             olen=len(old_file),
                                             ns=self.seguid(),
                                             os=old_file.seguid(),
                                             ntop={True:"-", False:"o"}[self.linear],
                                             otop={True:"-", False:"o"}[old_file.linear])

                    with open(filename, "w") as fp: fp.write(self.format(f))
        else:
            raise Exception("filename has to be a string, got", type(filename))
        return _display_html(_HTML(msg))
    
    def find(self, other):
        # TODO allow strings, seqs, seqrecords or Dseqrecords 
        # TODO check for linearity of other, raise exception if not
        # TODO add tests and docstring for this method   
        o = str(other.seq).upper()
        
        if self.linear:
            s = str(self.seq).upper()
        else:
            s = str(self.seq).upper()+str(self.seq).upper()[:len(other)-1] #allow wrapping around origin
        return s.find(o)


    def __str__(self):
        return ( "Dseqrecord\n"
                 "circular: {}\n"
                 "size: {}\n").format(self.circular, len(self))+_SeqRecord.__str__(self)
    def __contains__(self, other):
        if other.lower() in str(self.seq).lower():
            return True
        else:
            s = self.seq.watson.replace(" ","")
            ln  = len(s)
            spc = 3-ln%3 if ln%3 else 0
            s = "n" * spc + s + "nnn"
            for frame in range(3):
                if other.lower() in _translate(s[frame:frame+spc+ln]).lower():
                    return True
        return False

    def find_aa(self, other):
        return self.find_aminoacids(other)

    def find_aminoacids(self, other):
        '''
        >>> from pydna.dseqrecord import Dseqrecord
        >>> s=Dseqrecord("atgtacgatcgtatgctggttatattttag")
        >>> s.seq.translate()
        Seq('MYDRMLVIF*', HasStopCodon(ExtendedIUPACProtein(), '*'))
        >>> "RML" in s
        True
        >>> "MMM" in s
        False
        >>> s.seq.rc().translate()
        Seq('LKYNQHTIVH', ExtendedIUPACProtein())
        >>> "QHT" in s.rc()
        True
        >>> "QHT" in s
        False
        >>> slc = s.find_aa("RML")
        >>> slc
        slice(9, 18, None)
        >>> s[slc]
        Dseqrecord(-9)
        >>> code = s[slc].seq
        >>> code
        Dseq(-9)
        cgtatgctg
        gcatacgac
        >>> code.translate()
        Seq('RML', ExtendedIUPACProtein())
        '''
        other = str(other).lower()
        assert self.seq.watson == "".join(self.seq.watson.split())
        s = self.seq.watson
        ln  = len(s)
        spc = 3-ln%3 if ln%3 else 0
        s = s + "n"*spc + "nnn"
        start = None
        for frame in range(3):
            try:
                start = _translate(s[frame:frame+ln+spc]).lower().index(other)
                break
            except ValueError:
                pass
        oh = self.seq.ovhg if self.seq.ovhg>0 else 0
        if start == None:
            return None
        else:
            return slice(frame+start*3+oh, frame+(start+len(other))*3+oh)

    def map_trace_files(self, pth):
        import glob
        traces = []
        stf = _SequenceTraceFactory()
        for name in glob.glob(pth):
            traces.append( stf.loadTraceFile( name ))
        if not traces:
            raise Exception
        if hasattr( self.map_target, "step" ):
            area = self.map_target
        elif hasattr( self.map_target, "extract" ):
            area = slice(self.map_target.location.start, self.map_target.location.end)
        else:
            area = None

        if area:
            self.matching_reads     = []
            self.not_matching_reads = []
            target    = str(self[area].seq).lower()
            target_rc = str(self[area].seq.rc()).lower()
            for trace in traces:
                if target in str(trace.basecalls).lower() or target_rc in str(trace.basecalls).lower():
                   self.matching_reads.append(trace)
                else:
                    self.not_matching_reads.append(trace)
            reads = self.matching_reads
        else:
            self.matching_reads     = None
            self.not_matching_reads = None
            reads = traces

        for read_ in reads:
            matches = _common_sub_strings(str(self.seq).lower(), read_.basecalls.lower(), 25)
            if len(matches)>1:
                newmatches = [matches[0],]
                for i, x in enumerate(matches[1:]):
                    g,f,h = matches[i]
                    if g+h < x[0] and f+h < x[1]:
                        newmatches.append(x)
            elif len(matches)==1:
                newmatches = matches
            else:
                continue

            if len(newmatches)>1:
                ms = []
                for m in newmatches:
                    ms.append(_FeatureLocation(m[0], m[0]+m[2]))
                loc = _CompoundLocation(ms)
            else:
                a,b,c = newmatches[0]
                loc = _FeatureLocation(a,a+c)

            self.features.append( _SeqFeature(loc,
                                             qualifiers = {"label": read_.getFileName()},
                                             type="trace") )

        return [x.getFileName() for x in reads]

    def __repr__(self):
        return "Dseqrecord({}{})".format({True:"-", False:"o"}[self.linear],len(self))

    def _repr_pretty_(self, p, cycle):
            p.text("Dseqrecord({}{})".format({True:"-", False:"o"}[self.linear],len(self)))

    def __add__(self, other):
        if hasattr(other, "seq") and hasattr(other.seq, "watson"):
            offset = other.seq.ovhg
            other.features = [f._shift(offset) for f in other.features]
            #answer = self.__class__(_SeqRecord.__add__(self, other))
            answer = Dseqrecord(_SeqRecord.__add__(self, other))
            answer.n = min(self.n, other.n)
        else:
            #answer = self.__class__(_SeqRecord.__add__(self, Dseqrecord(other)))
            answer = Dseqrecord(_SeqRecord.__add__(self, Dseqrecord(other)))
            answer.n = self.n
        return answer

    def __mul__(self, number):
        if not isinstance(number, int):
            raise TypeError("TypeError: can't multiply Dseqrecord by non-int of type {}".format(type(number)))
        if self.circular:
            raise TypeError("TypeError: can't multiply circular Dseqrecord")
        if number>0:
            new = _copy.copy(self)
            for i in range(1, number):
                new += self
            return new
        else:
            return self.__class__("")

    def __getitem__(self, sl):
        answer = Dseqrecord(_copy.copy(self))        
        answer.seq = answer.seq.__getitem__(sl)
        answer.seq.alphabet = self.seq.alphabet

        sl_start = sl.start if sl.start is not None else 0
        sl_stop = sl.stop if sl.stop is not None else len(answer.seq)
        
        if self.linear or sl_start<sl_stop:
            answer.features = _SeqRecord.__getitem__(self, sl).features
        else:
            try:
                answer.features = self.shifted(sl.stop).features
            except Exception:
                answer.features = self.features
            answer.features = [f for f in answer.features if f.location.parts == sorted(f.location.parts)]
        return answer

    def __eq__( self, other ):
        try:
            if self.seq == other.seq and str(self.__dict__) == str(other.__dict__):
                return True
        except AttributeError:
            pass
        return False

    def __ne__( self, other ):
        return not self.__eq__(other)

    def __hash__(self):
        """__hash__ must be based on __eq__"""
        return hash( (str(self.seq).lower(), str(tuple(sorted(self.__dict__.items())))))

    def linearize(self, *enzymes):
        '''This method is similar to :func:`cut` but throws an exception if there
        is not excactly on cut i.e. none or more than one digestion products.

        '''

        if self.seq._linear:
            raise Exception("Can only linearize circular molecules!")
        fragments = self.cut(*enzymes)
        if len(fragments)>1:
            raise Exception("More than one fragment is formed!")
        if not fragments:
            raise Exception("The enzyme(s) do not cut!")
        return fragments.pop()

    def cut(self, *enzymes):
        '''Digest the Dseqrecord object with one or more restriction enzymes.
        returns a list of linear Dseqrecords. If there are no cuts, an empty
        list is returned.

        See also :func:`Dseq.cut`

        Parameters
        ----------

        enzymes : enzyme object or iterable of such objects
            A Bio.Restriction.XXX restriction object or iterable of such.

        Returns
        -------
        Dseqrecord_frags : list
            list of Dseqrecord objects formed by the digestion

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("ggatcc")
        >>> from Bio.Restriction import BamHI
        >>> a.cut(BamHI)
        [Dseqrecord(-5), Dseqrecord(-5)]
        >>> frag1, frag2 = a.cut(BamHI)
        >>> frag1.seq
        Dseq(-5)
        g
        cctag
        >>> frag2.seq
        Dseq(-5)
        gatcc
            g


        '''

        output, stack = [], []
        stack.extend(reversed(enzymes))
        while stack:
            top = stack.pop()
            if hasattr(top, "__iter__"):
                stack.extend(reversed(top))
            else:
                output.append(top)
        enzymes = output
        if not hasattr(enzymes, '__iter__'):
            enzymes = (enzymes,)

        frags = self.seq.cut(enzymes)

        if not frags:
            return []

        if self.linear:
            last_pos=0
            #template = self.__class__(self, linear=True)
            #template = _copy.copy(self)
            template = self
        else:
            last_pos = [p.pop(0)-max(0,e.ovhg)-1 for (p,e) in
                         sorted([(enzyme.search(_Seq(self.seq.dsdata),
                                                linear = self.linear)[::-1],
                                   enzyme) for enzyme in enzymes]) if p]
            if not last_pos:
                return [self]
            if 0 in last_pos:
                last_pos=0
            else:
                last_pos = last_pos.pop()
            template = self._multiply_circular(3)

        Dseqrecord_frags = []
        start = last_pos

        for f in frags:

            end = start + len(str(f))
            Dseqrecord_frag = Dseqrecord(f, linear=True, n=self.n)

            Dseqrecord_frag.features = template[start:end].features
            Dseqrecord_frag.annotations         = _copy.copy(self[start:end].annotations)
            Dseqrecord_frag.name                = _copy.copy(self.name)
            Dseqrecord_frag.dbxrefs             = _copy.copy(self[start:end].dbxrefs)
            Dseqrecord_frag.id                  = _copy.copy(self.id)
            Dseqrecord_frag.letter_annotations  = _copy.copy(self[start:end].letter_annotations)

            Dseqrecord_frag.description = self.description+"_"+"_".join(str(e) for e in enzymes)

            Dseqrecord_frags.append(Dseqrecord_frag)
            start = end
            start-= len(f.three_prime_end()[1])

        return Dseqrecord_frags

    def number_of_cuts(self, *enzymes):
        output = []
        stack = []
        stack.extend(reversed(enzymes))
        while stack:
            top = stack.pop()
            if hasattr(top, "__iter__"):
                stack.extend(reversed(top))
            else:
                output.append(top)
        enzymes = output
        if not hasattr(enzymes, '__iter__'):
            enzymes = (enzymes,)
        return sum([len(enzyme.search(self.seq)) for enzyme in enzymes])

    def reverse_complement(self):
        '''Returns the reverse complement.

        Examples
        --------
        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("ggaatt")
        >>> a
        Dseqrecord(-6)
        >>> a.seq
        Dseq(-6)
        ggaatt
        ccttaa
        >>> a.reverse_complement().seq
        Dseq(-6)
        aattcc
        ttaagg
        >>>

        See also
        --------
        pydna.dseq.Dseq.reverse_complement

        '''

        return self.rc()

    def rc(self):
        '''alias of the reverse_complement method'''
        answer = Dseqrecord(super(Dseqrecord, self).reverse_complement())
        assert answer.circular == self.circular
        answer.name       = "{}_rc".format(self.name[:13])
        answer.description= self.description+"_rc"
        answer.id         = self.id+"_rc"
        return answer
        #return Dseqrecord(self.seq.rc())


    def _multiply_circular(self, multiplier):
        '''returns a linearised version of a circular sequence multiplied by
       multiplier '''

        if self.linear:
            raise TypeError("sequence has to be circular!")
        if not isinstance(multiplier, int):
            raise TypeError("TypeError: can't multiply Dseq by non-int of type {}".format(type(multiplier)))
        if multiplier<=0:
            return self.__class__("")

        new_features = []

        for feature in self.features:
            new_feature = _copy.deepcopy(feature)
            if len(new_feature.location.parts)>1:    # CompoundFeature
                j=0
                while (j+1)<=len(new_feature.location.parts):
                    if new_feature.location.parts[j].end == len(self) and new_feature.location.parts[j+1].start==0:
                        new_feature.location.parts[j] = _FeatureLocation(new_feature.location.parts[j].start,
                                                                        new_feature.location.parts[j].end+len(new_feature.location.parts[j+1]))
                        del new_feature.location.parts[j+1]
                    j+=1
                slask = [new_feature.location.parts.pop(0)]
                for fl in new_feature.location.parts:
                    if fl.start < slask[-1].start:
                        slask.append(fl+len(self))
                    else:
                        slask.append(fl)
                if len(slask)>1:
                    new_feature.location.parts=slask
                else:
                    new_feature.location=slask[0]
            new_features.append(new_feature)

        sequence = self.tolinear()
        sequence.features = new_features
        sequence = sequence * multiplier

        sequence.features = [f for f in sequence.features if f.location.end <= len(sequence)]
        sequence.features.sort(key = _operator.attrgetter('location.start'))
        return sequence

    def shifted(self, shift):
        '''Returns a circular Dseqrecord with a new origin <shift>.
         This only works on circular Dseqrecords. If we consider the following
         circular sequence:


         | ``GAAAT   <-- watson strand``
         | ``CTTTA   <-- crick strand``

         The T and the G on the watson strand are linked together as well
         as the A and the C of the of the crick strand.

         if ``shift`` is 1, this indicates a new origin at position 1:

         |    new origin at the | symbol:
         |
         | ``G|AAAT``
         | ``C|TTTA``

         new sequence:

         | ``AAATG``
         | ``TTTAC``

         Shift is always positive and 0<shift<length, so in the example
         below, permissible values of shift are 1,2 and 3

         Examples
         --------

         >>> from pydna.dseqrecord import Dseqrecord
         >>> a=Dseqrecord("aaat",circular=True)
         >>> a
         Dseqrecord(o4)
         >>> a.seq
         Dseq(o4)
         aaat
         ttta
         >>> b=a.shifted(1)
         >>> b
         Dseqrecord(o4)
         >>> b.seq
         Dseq(o4)
         aata
         ttat

        '''
        if self.linear:
            raise Exception("Sequence is linear.\n"
                             "The origin can only be\n"
                             "shifted for a circular sequence!\n")

        length=len(self)

        if not 0<shift<length:
            raise Exception("shift is {}, has to be 0<shift<{}".format(shift, length))

        new = self._multiply_circular(3)[shift:]

        features_to_fold = [f for f in new.features if f.location.start<length<f.location.end]

        folded_features = []

        for feature in features_to_fold:

            if len(feature.location.parts)>1: # CompoundFeature
                nps=[]
                for part in feature.location.parts:

                    if part.start<part.end<=length:
                        nps.append(part)

                    elif part.start<length<part.end:
                        nps.append(_FeatureLocation(part.start,length))
                        nps.append(_FeatureLocation(0, part.end-length))

                    elif length<=part.start<part.end:
                        nps.append(_FeatureLocation(part.start-length, part.end-length))

                folded_features.append(_SeqFeature(_CompoundLocation(nps),
                                                  qualifiers = feature.qualifiers,
                                                  type=feature.type))

            else:
                folded_features.append(_SeqFeature(_CompoundLocation([_FeatureLocation(feature.location.start, length),
                                                                    _FeatureLocation(0, feature.location.end-length)]),
                                                  qualifiers = feature.qualifiers,
                                                  type=feature.type))

        new = new[:length].looped()
        new.features.extend(folded_features)
        new.features.sort(key = _operator.attrgetter('location.start'))
        new.description = self.description #!
        return new

    @_memorize("Dseqrecord_synced")
    def synced(self, ref, limit = 25):
        '''This function returns a new circular sequence (Dseqrecord object), which has been rotated
        in such a way that there is maximum overlap between the sequence and
        ref, which may be a string, Biopython Seq, SeqRecord object or
        another Dseqrecord object.

        The reason for using this could be to rotate a recombinant plasmid so
        that it starts at the same position after cloning. See the example below:


        Examples
        --------

        >>> from pydna.dseqrecord import Dseqrecord
        >>> a=Dseqrecord("gaat",circular=True)
        >>> a.seq
        Dseq(o4)
        gaat
        ctta
        >>> d = a[2:] + a[:2]
        >>> d.seq
        Dseq(-4)
        atga
        tact
        >>> insert=Dseqrecord("CCC")
        >>> recombinant = (d+insert).looped()
        >>> recombinant.seq
        Dseq(o7)
        atgaCCC
        tactGGG
        >>> recombinant.synced(a).seq
        Dseq(o7)
        gaCCCat
        ctGGGta

        '''

        if self.linear:
            raise Exception("Only circular DNA can be synced!")
        try:
            rs = ref.seguid()
        except AttributeError:
            rs = _seg(ref)

        newseq = _copy.copy(self)

        s    = str(self.seq.watson).lower()
        s_rc = str(self.seq.crick).lower()

        if hasattr(ref, "seq"):
            r=ref.seq
            if hasattr(ref, "watson"):
                r = str(r.watson).lower()
            else:
                r = str(r).lower()
        else:
            r = str(ref.lower())

        lim = min(limit, limit*(len(s)//limit)+1)

        c = _common_sub_strings(s+s,       r, limit = lim)
        d = _common_sub_strings(s_rc+s_rc, r, limit = lim)

        c =  [(x[0],x[2]) for x in c if x[1]==0]
        d =  [(x[0],x[2]) for x in d if x[1]==0]

        if not c and not d:
            raise Exception("There is no overlap between sequences!")

        if c:
            start, length = c.pop(0)
        else:
            start, length = 0,0

        if d:
            start_rc, length_rc = d.pop(0)
        else:
            start_rc, length_rc = 0,0

        if length_rc>length:
            start = start_rc
            newseq = newseq.rc()

        if start == 0:
            result = newseq
        else:
            result = newseq.shifted(start)
        _module_logger.info("synced")
        return result


    
if __name__=="__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"]=cache
