#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""Provides two functions, parse and parse_primers"""

import os as _os
import re as _re
import io as _io
import textwrap as _textwrap

# import glob      as _glob

from Bio import SeqIO as _SeqIO

# from Bio.Alphabet import generic_dna as _generic_dna
from pydna.genbankfile import GenbankFile as _GenbankFile
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.primer import Primer as _Primer


def parse(data, ds=True):
    """This function returns *all* DNA sequences found in data. If no
    sequences are found, an empty list is returned. This is a greedy
    function, use carefully.

    Parameters
    ----------
    data : string or iterable
        The data parameter is a string containing:

        1. an absolute path to a local file.
           The file will be read in text
           mode and parsed for EMBL, FASTA
           and Genbank sequences. Can be
           a string or a Path object.

        2. a string containing one or more
           sequences in EMBL, GENBANK,
           or FASTA format. Mixed formats
           are allowed.

        3. data can be a list or other iterable where the elements are 1 or 2

    ds : bool
        If True double stranded :class:`Dseqrecord` objects are returned.
        If False single stranded :class:`Bio.SeqRecord` [#]_ objects are returned.

    Returns
    -------
    list
        contains Dseqrecord or SeqRecord objects

    References
    ----------

    .. [#] http://biopython.org/wiki/SeqRecord

    See Also
    --------
    read

    """

    def embl_gb_fasta(raw, ds, path=None):

        pattern = r"(?:>.+\n^(?:^[^>]+?)(?=\n\n|>|LOCUS|ID))|(?:(?:LOCUS|ID)(?:(?:.|\n)+?)^//)"

        result_list = []

        rawseqs = _re.findall(
            pattern, _textwrap.dedent(raw + "\n\n"), flags=_re.MULTILINE
        )

        for rawseq in rawseqs:
            handle = _io.StringIO(rawseq)
            circular = False
            try:
                parsed = _SeqIO.read(handle, "embl")
            except ValueError:
                handle.seek(0)
                try:
                    parsed = _SeqIO.read(handle, "genbank")
                    if "circular" in str(parsed.annotations.get("topology")).lower():
                        circular = True
                except ValueError:
                    handle.seek(0)
                    try:
                        parsed = _SeqIO.read(handle, "fasta")
                    except ValueError:
                        parsed = ""
            handle.close()
            if (
                "circular" in rawseq.splitlines()[0].lower().split()
            ):  # hack to pick up topology from malformed files
                circular = True
            if parsed:
                from copy import deepcopy as _deepcopy  # TODO: clean up !
                from pydna.seqfeature import SeqFeature as _SeqFeature

                nfs = [_SeqFeature() for f in parsed.features]
                for f, nf in zip(parsed.features, nfs):
                    nf.__dict__ = _deepcopy(f.__dict__)
                parsed.features = nfs
                if ds and path:
                    result_list.append(
                        _GenbankFile.from_SeqRecord(
                            parsed, linear=not circular, circular=circular, path=path
                        )
                    )
                elif ds:
                    result_list.append(
                        _Dseqrecord.from_SeqRecord(
                            parsed, linear=not circular, circular=circular
                        )
                    )
                else:
                    result_list.append(parsed)

        return result_list

    # a string is an iterable datatype but on Python2.x
    # it doesn't have an __iter__ method.
    if not hasattr(data, "__iter__") or isinstance(data, (str, bytes)):
        data = (data,)

    sequences = []

    for item in data:
        try:
            # item is a path to a utf-8 encoded text file?
            with open(item, "r", encoding="utf-8") as f:
                raw = f.read()
        except IOError:
            # item was not a path, add sequences parsed from item
            raw = item
            path = None
        else:
            # item was a readable text file, seqences are parsed from the file
            path = item
        finally:
            sequences.extend(embl_gb_fasta(raw, ds, path))
    return sequences


def parse_primers(data):
    """ """
    return [_Primer(x) for x in parse(data, ds=False)]


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
