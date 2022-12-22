#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Miscellaneous functions."""

from Bio.Data.IUPACData import (ambiguous_dna_complement as
                                _ambiguous_dna_complement)
from Bio.Seq import _maketrans
from pydna._pretty import pretty_str as _pretty_str
from Bio.SeqUtils.CheckSum import seguid as _seguid
import shelve as _shelve
import os as _os
import re as _re
import logging as _logging
import base64 as _base64
import pickle as _pickle
import hashlib as _hashlib
import keyword as _keyword
import collections as _collections
import itertools as _itertools

import sys as _sys
import re
import textwrap
import itertools
import math
import random
import subprocess as _subprocess

from pydna.codon import weights as _weights
from pydna.codon import rare_codons as _rare_codons
from pydna.codon import n_end as _n_end
from Bio.SeqUtils import seq3 as _seq3
from pydna._pretty import PrettyTable as _PrettyTable

_module_logger = _logging.getLogger("pydna." + __name__)
_ambiguous_dna_complement.update({"U": "A"})
_complement_table = _maketrans(_ambiguous_dna_complement)


def smallest_rotation(s):
    """Smallest rotation of a string.

    Algorithm described in Pierre Duval, Jean. 1983. Factorizing Words
    over an Ordered Alphabet. Journal of Algorithms & Computational Technology
    4 (4) (December 1): 363–381. and Algorithms on strings and sequences based
    on Lyndon words, David Eppstein 2011.
    https://gist.github.com/dvberkel/1950267

    Examples
    --------
    >>> from pydna.utils import smallest_rotation
    >>> smallest_rotation("taaa")
    'aaat'

    """
    prev, rep = None, 0
    ds = 2 * s
    lens = len(s)
    lends = len(ds)
    old = 0
    k = 0
    w = ""
    while k < lends:
        i, j = k, k + 1
        while j < lends and ds[i] <= ds[j]:
            i = (ds[i] == ds[j]) and i + 1 or k
            j += 1
        while k < i + 1:
            k += j - i
            prev = w
            w = ds[old:k]
            old = k
            if w == prev:
                rep += 1
            else:
                prev, rep = w, 1
            if len(w) * rep == lens:
                return w * rep


def cai(seq: str,
        organism: str = "sce"):
    """docstring."""
    from CAI import CAI as _CAI
    return round(_CAI(seq.upper(), weights=_weights[organism]), 3)


def rarecodons(seq: str,
               organism="sce"):
    """docstring."""
    rare = _rare_codons[organism]
    s = seq.upper()
    slices = []
    for i in range(0, len(seq)//3):
        x, y = i*3, i*3+3
        trip = s[x:y]
        if trip in rare:
            slices.append(slice(x, y, 1))
    return slices


def express(seq: str, organism="sce"):
    """docstring.


    **NOT IMPLEMENTED YET**
    """
    x = _PrettyTable(["cds", "len", "cai", "gc", "sta", "stp",
                      "n-end"]+_rare_codons[organism]+["rare"])
    val = []

    val.append(f"{self._data.upper().decode('ASCII')[:3]}..."
               f"{self._data.upper().decode('ASCII')[-3:]}")
    val.append(len(self)/3)
    val.append(cai(organism))
    val.append(gc())
    val.append(startcodon())
    val.append(stopcodon())
    val.append(_n_end[organism].get(_seq3(self[3:6].translate())))
    s = self._data.upper().decode("ASCII")
    trps = [s[i*3:i*3+3] for i in range(0, len(s)//3)]
    tot = 0
    for cdn in _rare_codons[organism]:
        cnt = trps.count(cdn)
        tot += cnt
        val.append(cnt)
    val.append(round(tot/len(trps), 3))
    x.add_row(val)
    return x


def open_folder(pth):
    """docstring."""
    if _sys.platform == "win32":
        _subprocess.run(["start", pth], shell=True)
    elif _sys.platform == "darwin":
        _subprocess.run(["open", pth])
    else:
        try:
            _subprocess.run(["xdg-open", pth])
        except OSError:
            return "no cache to open."


def rc(sequence: str):
    """Reverse complement.

    accepts mixed DNA/RNA
    """
    return sequence.translate(_complement_table)[::-1]


def complement(sequence: str):
    """Complement.

    accepts mixed DNA/RNA
    """
    return sequence.translate(_complement_table)


def memorize(filename):
    """Cache functions and classes.

    see pydna.download and pydna.Assembly for use
    """

    def decorator(f):
        def wrappee(*args, **kwargs):
            _module_logger.info("#### memorizer ####")
            _module_logger.info("cache filename                   = %s", filename)
            _module_logger.info(
                "os.environ['pydna_cached_funcs'] = %s",
                _os.environ["pydna_cached_funcs"],
            )
            if filename not in _os.environ["pydna_cached_funcs"]:
                _module_logger.info(
                    "cache filename not among cached functions, made it new!"
                )
                return f(*args, **kwargs)
            key = _base64.urlsafe_b64encode(
                _hashlib.sha1(_pickle.dumps((args, kwargs))).digest()
            ).decode("ascii")
            _module_logger.info("key = %s", key)
            cache = _shelve.open(
                _os.path.join(
                    _os.environ["pydna_data_dir"], identifier_from_string(filename)
                ),
                writeback=False,
            )
            try:
                result = cache[key]
            except KeyError:
                _module_logger.info(
                    "no result for key %s in shelve %s",
                    key,
                    identifier_from_string(filename),
                )
                result = f(*args, **kwargs)
                _module_logger.info("made it new!")
                cache[key] = result
                _module_logger.info("saved result under key %s", key)
            else:
                _module_logger.info("found %s in cache", key)
            cache.close()
            return result

        return wrappee

    return decorator


def identifier_from_string(s: str) -> str:
    """Return a valid python identifier.

    based on the argument s or an empty string
    """
    s = s.strip()
    s = _re.sub(r"\s+", r"_", s)
    s.replace("-", "_")
    s = _re.sub("[^0-9a-zA-Z_]", "", s)
    if s and not s[0].isidentifier() or _keyword.iskeyword(s):
        s = "_{s}".format(s=s)
    assert s == "" or s.isidentifier()
    return s


def seguid(seq: str) -> _pretty_str:
    """SEGUID checksum for a string representing a biological sequence.

    This is the SEGUID checksum with the standard Base64
    encoding that can contain '+' and '/' as originally
    defined by Babnigg and Giometti:

    Babnigg, Giometti 2006. “A Database of Unique
    Protein Sequence Identifiers for Proteome Studies.” *Proteomics* 6
    (16) (August): 4514–4522.

    Examples
    --------
    >>> from pydna.utils import seguid
    >>> seguid("aaa")
    'YG7G6b2Kj/KtFOX63j8mRHHoIlE'
    """
    return _pretty_str(_seguid(seq.upper()))

def useguid(seq: str) -> _pretty_str:
    """Returns the url safe SEGUID checksum for the sequence.
    This is the SEGUID checksum with the '+' and '/' characters of standard
    Base64 encoding are respectively replaced by '-' and '_'.

    Examples
    --------
    >>> from pydna.utils import useguid
    >>> useguid("aaa")
    'YG7G6b2Kj_KtFOX63j8mRHHoIlE'
    """
    return seguid(seq).replace("+", "-").replace("/", "_")


def lseguid_blunt(seq: str) -> _pretty_str:
    """lSEGUID checksum.

    for a string representing a blunt double stranded DNA molecule.

    Examples
    --------
    >>> from pydna.utils import lseguid_blunt
    >>> lseguid_blunt("ttt")
    'YG7G6b2Kj_KtFOX63j8mRHHoIlE'
    >>> lseguid_blunt("aaa")
    'YG7G6b2Kj_KtFOX63j8mRHHoIlE'
    """
    return useguid(min(seq.upper(), str(rc(seq)).upper()))


def lseguid_sticky(watson: str, crick: str, overhang: int) -> _pretty_str:
    """Linear SEGUID (lSEGUID) checksum.

    Calculates the lSEGUID checksum for a double stranded DNA sequence
    described by two strings (watson and crick) representing the two
    complementary DNA strands and an integer describing the stagger
    between the two strands in the 5' end.

    The overhang is defined as the amount of 3' overhang in the 5'
    side of the molecule. A molecule with 5' overhang has a negative
    value.

        dsDNA    ovhg

          nnn...    2
        nnnnn...

          nnnn...    1
        nnnnn...

        nnnnn...    0
        nnnnn...

        nnnnn...   -1
          nnnn...

        nnnnn...   -2
          nnn...


    """
    watson = watson.upper()
    crick = crick.upper()
    lw = len(watson)
    lc = len(crick)
    if overhang == 0 and lw == lc:
        return lseguid_blunt(watson)
    else:
        w, c, o = min(((watson, crick, overhang),
                       (crick, watson, lw - lc + overhang)))

    return useguid(f"{o*chr(32)}{w}\n{-o*chr(32)}{c[::-1]}")


def cseguid(seq: str) -> _pretty_str:
    """Url safe cSEGUID for a string representing a circular double stranded
    DNA molecule.

    The cSEGUID is the SEGUID checksum calculated for the lexicographically
    minimal string rotation of a DNA sequence. Only defined for circular
    sequences.

    Examples
    --------
    >>> from pydna.utils import cseguid
    >>> cseguid("attt")
    'oopV-6158nHJqedi8lsshIfcqYA'
    >>> cseguid("ttta")
    'oopV-6158nHJqedi8lsshIfcqYA'
    """
    return useguid(min(smallest_rotation(seq.upper()),
                       smallest_rotation(str(rc(seq)).upper())))


def flatten(*args):
    """Flattens an iterable of iterables.

    Down to str, bytes, bytearray or any of the pydna or Biopython seq objects
    """
    output = []
    args = list(args)
    while args:
        top = args.pop()
        if (
                isinstance(top, _collections.abc.Iterable)
                and not isinstance(top, (str, bytes, bytearray))
                and not hasattr(top, "reverse_complement")
                ):
            args.extend(top)
        else:
            output.append(top)
    return output[::-1]


def seq31(seq):
    """Turn a three letter code protein sequence into one with one letter code.

    The single input argument 'seq' should be a protein sequence using single
    letter codes, as a python string.

    This function returns the amino acid sequence as a string using the one
    letter amino acid codes. Output follows the IUPAC standard (including
    ambiguous characters B for "Asx", J for "Xle" and X for "Xaa", and also U
    for "Sel" and O for "Pyl") plus "Ter" for a terminator given as an
    asterisk.

    Any unknown
    character (including possible gap characters), is changed into 'Xaa'.

    Examples
    --------
    >>> from Bio.SeqUtils import seq3
    >>> seq3("MAIVMGRWKGAR*")
    'MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer'
    >>> from pydna.utils import seq31
    >>> seq31('MetAlaIleValMetGlyArgTrpLysGlyAlaArgTer')
    'M  A  I  V  M  G  R  W  K  G  A  R  *'
    """
    threecode = {
        "Ala": "A",
        "Asx": "B",
        "Cys": "C",
        "Asp": "D",
        "Glu": "E",
        "Phe": "F",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Lys": "K",
        "Leu": "L",
        "Met": "M",
        "Asn": "N",
        "Pro": "P",
        "Gln": "Q",
        "Arg": "R",
        "Ser": "S",
        "Thr": "T",
        "Val": "V",
        "Trp": "W",
        "Tyr": "Y",
        "Glx": "Z",
        "Xaa": "X",
        "Ter": "*",
        "Sel": "U",
        "Pyl": "O",
        "Xle": "J",
    }

    nr_of_codons = int(len(seq) / 3)
    sequence = [seq[i * 3 : i * 3 + 3].title() for i in range(nr_of_codons)]
    padding = " " * 2
    return padding.join([threecode.get(aa, "X") for aa in sequence])


def parse_text_table(rawtable, tabs=4):

    table = textwrap.dedent(rawtable.expandtabs(tabs)).strip()
    max_row_length = max([len(row.strip()) for row in table.splitlines()])
    rows = [row.ljust(max_row_length) for row in table.splitlines()]
    table = "\n".join(rows)
    empty_column_regex = r"(?:\n?\s{%s}\n)+" % len(rows)
    transposed_table = "\n".join(["".join(c) for c in zip(*rows)])
    cols = re.split(empty_column_regex, transposed_table)
    list_of_lists_cr = []

    for col in cols:
        columnlist = [
            "".join(c).strip() for c in zip(*[a for a in col.splitlines() if a])
        ]
        maxlen = max([len(c) for c in columnlist])
        columnlist = [c.ljust(maxlen) for c in columnlist]
        list_of_lists_cr.append(columnlist)

    list_of_lists_rc = [list(i) for i in zip(*list_of_lists_cr)]

    formatted = _pretty_str("\n".join(" ".join(cell) for cell in
                            list_of_lists_rc))

    columnsplit = _pretty_str("\n|||\n".join(
        [
            "\n".join(
                [
                    x.strip()
                    for x in ["".join(c) for c in zip(*col.strip("\n").splitlines())]
                ]
            )
            for col in cols
        ]
    ))

    rowsplit = "\n---\n".join(["\n".join(a).strip() for a in zip(*list_of_lists_cr)])
    rowsplit = _pretty_str("\n".join(row.strip() for row in rowsplit.splitlines()))

    return (formatted,
            columnsplit,
            rowsplit,
            list_of_lists_rc,
            list_of_lists_cr,)


def join_list_to_table(rawlist):

    if "|||\n" in rawlist:
        raw_columns = rawlist.split("|||\n")
        cols = [col.splitlines() for col in raw_columns]
    elif "---\n" in rawlist:
        rawrows = rawlist.split("---\n")
        rows = [row.splitlines() for row in rawrows]
        cols = list(itertools.zip_longest(*rows, fillvalue=""))
    else:
        return None

    number_of_rows = max([len(col) for col in cols])
    formatted_cols = []

    for col in cols:
        # print col
        rows = [row.strip() or "\"" for row in col]
        width = max([len(row) for row in rows])

        rows = [row.ljust(width) for row in rows]
        rows += [" " * width] * (number_of_rows - len(rows))
        formatted_cols.append(rows)

    rows = list(zip(*formatted_cols))

    combinedlist = []

    for row in rows:
        combinedlist.append(" ".join(row))

    new_text = "\n".join(combinedlist)

    return _pretty_str(new_text)


def expandtolist(content):
    """docstring."""
    resultlist = []
    for line in re.finditer(r"(?P<item>[^\(\)]*?)(?P<brack>\[.*?\])", content):
        text2rep = line.group("item")
        bracket = line.group("brack")
        padding = max(
            [len(str(x).strip()) for x in re.split(r"\.\.|,",
                                                   bracket.strip("[ ]"))]
        )
        inbracket = [item.strip("[ ]") for item in bracket.split(",")]
        expanded = []

        for item in inbracket:
            if re.match(r"(\d+\.\.\d+)|([a-z]+\.\."
                        r"[a-z]+)|([A-Z]+\.\.[A-Z]+)", item):
                low, high = item.split(
                    "..",
                )
                if low.isdigit() and high.isdigit():
                    r = [
                        "{:{}d}".format(x, padding)
                        for x in range(int(low), 1 + int(high))
                    ]
                if (low.islower() and high.islower()) or (
                    low.isupper() and high.isupper()
                ):
                    r = [chr(a) for a in range(ord(low), 1 + ord(high))]
                expanded.extend(r)
            else:
                expanded.append(item.strip())

        resultlist.append([text2rep + x for x in expanded])

    ml = max([len(x) for x in resultlist])

    norm = []
    for r in resultlist:
        mp = int(math.ceil(float(ml) / float(len(r))))
        norm.append(list(itertools.chain.from_iterable(list(zip(*(r,) * mp)))))

    rt = ""
    for a in range(ml):
        rt += "".join([b[a] for b in norm]) + "\n"
    return _pretty_str(rt)


def randomRNA(length, maxlength=None):
    """docstring."""
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength)))
    return "".join([random.choice("GAUC") for x in range(length)])


def randomDNA(length, maxlength=None):
    """docstring."""
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength)))
    return "".join([random.choice("GATC") for x in range(length)])


def randomORF(length, maxlength=None):
    """docstring."""
    length -= 2
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength - 2)))

    cdns = (
        "TTT",
        "TTC",
        "TTA",
        "TTG",
        "TCT",
        "TCC",
        "TCA",
        "TCG",
        "TAT",
        "TAC",
        "TGT",
        "TGC",
        "TGG",
        "CTT",
        "CTC",
        "CTA",
        "CTG",
        "CCT",
        "CCC",
        "CCA",
        "CCG",
        "CAT",
        "CAC",
        "CAA",
        "CAG",
        "CGT",
        "CGC",
        "CGA",
        "CGG",
        "ATT",
        "ATC",
        "ATA",
        "ATG",
        "ACT",
        "ACC",
        "ACA",
        "ACG",
        "AAT",
        "AAC",
        "AAA",
        "AAG",
        "AGT",
        "AGC",
        "AGA",
        "AGG",
        "GTT",
        "GTC",
        "GTA",
        "GTG",
        "GCT",
        "GCC",
        "GCA",
        "GCG",
        "GAT",
        "GAC",
        "GAA",
        "GAG",
        "GGT",
        "GGC",
        "GGA",
        "GGG",
    )

    starts = ("ATG",)
    stops = ("TAA", "TAG", "TGA")

    return (
        random.choice(starts)
        + "".join([random.choice(cdns) for x in range(length)])
        + random.choice(stops)
    )


def randomprot(length, maxlength=None):
    """docstring."""
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength)))
    return "".join([random.choice("ACDEFGHIKLMNPQRSTVWY")
                    for x in range(length)])


def eq(*args, **kwargs):
    """Compare two or more DNA sequences for equality.

    Compares two or more DNA sequences for equality i.e. if they
    represent the same double stranded DNA molecule.

    Parameters
    ----------
    args : iterable
        iterable containing sequences
        args can be strings, Biopython Seq or SeqRecord, Dseqrecord
        or dsDNA objects.
    circular : bool, optional
        Consider all molecules circular or linear
    linear : bool, optional
        Consider all molecules circular or linear

    Returns
    -------
    eq : bool
        Returns True or False

    Notes
    -----
    Compares two or more DNA sequences for equality i.e. if they
    represent the same DNA molecule.

    Two linear sequences are considiered equal if either:

    1. They have the same sequence (case insensitive)
    2. One sequence is the reverse complement of the other

    Two circular sequences are considered equal if they are circular
    permutations meaning that they have the same length and:

    1. One sequence can be found in the concatenation of the other sequence with itself.
    2. The reverse complement of one sequence can be found in the concatenation of the other sequence with itself.

    The topology for the comparison can be set using one of the keywords
    linear or circular to True or False.

    If circular or linear is not set, it will be deduced from the topology of
    each sequence for sequences that have a linear or circular attribute
    (like Dseq and Dseqrecord).

    Examples
    --------
    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.utils import eq
    >>> eq("aaa","AAA")
    True
    >>> eq("aaa","AAA","TTT")
    True
    >>> eq("aaa","AAA","TTT","tTt")
    True
    >>> eq("aaa","AAA","TTT","tTt", linear=True)
    True
    >>> eq("Taaa","aTaa", linear = True)
    False
    >>> eq("Taaa","aTaa", circular = True)
    True
    >>> a=Dseqrecord("Taaa")
    >>> b=Dseqrecord("aTaa")
    >>> eq(a,b)
    False
    >>> eq(a,b,circular=True)
    True
    >>> a=a.looped()
    >>> b=b.looped()
    >>> eq(a,b)
    True
    >>> eq(a,b,circular=False)
    False
    >>> eq(a,b,linear=True)
    False
    >>> eq(a,b,linear=False)
    True
    >>> eq("ggatcc","GGATCC")
    True
    >>> eq("ggatcca","GGATCCa")
    True
    >>> eq("ggatcca","tGGATCC")
    True
    """
    args = flatten(args)  # flatten

    topology = None

    if "linear" in kwargs:
        if kwargs["linear"] is True:
            topology = "linear"
        if kwargs["linear"] is False:
            topology = "circular"
    elif "circular" in kwargs:
        if kwargs["circular"] is True:
            topology = "circular"
        if kwargs["circular"] is False:
            topology = "linear"
    else:
        topology = set([arg.circular if
                        hasattr(arg, "circular") else None for arg in args])

        if len(topology) != 1:
            raise ValueError("sequences have different topologies")
        topology = topology.pop()
        if topology in (False, None):
            topology = "linear"
        elif topology is True:
            topology = "circular"

    args = [arg.seq if hasattr(arg, "seq") else arg for arg in args]
    args_string_list = [
        arg.watson.lower() if hasattr(arg, "watson") else str(arg).lower()
        for arg in args
    ]

    length = set((len(s) for s in args_string_list))

    if len(length) != 1:
        return False
    same = True

    if topology == "circular":
        # force circular comparison of all given sequences
        for s1, s2 in _itertools.combinations(args_string_list, 2):
            if not (s1 in s2 + s2 or rc(s1) in s2 + s2):
                same = False
    elif topology == "linear":
        # force linear comparison of all given sequences
        for s1, s2 in _itertools.combinations(args_string_list, 2):
            if not (s1 == s2 or s1 == rc(s2)):
                same = False
    return same

if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
