#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""This module provides miscellaneous functions."""

from Bio.Data.IUPACData import ambiguous_dna_complement as _ambiguous_dna_complement
from Bio.Seq import _maketrans
from pydna._pretty import pretty_str as _pretty_str
from Bio.SeqUtils.CheckSum import seguid as _base64_seguid
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

import re
import textwrap
import itertools
import math
import random

_module_logger = _logging.getLogger("pydna." + __name__)


# from Bio.Seq             import reverse_complement as _reverse_complement

# from Bio.Seq             import reverse_complement as _rc


_ambiguous_dna_complement.update({"U": "A"})
_complement_table = _maketrans(_ambiguous_dna_complement)


# Zhang, S. P., Zubay, G., & Goldman, E. (1991). Low-usage codons in Escherichia
# coli, yeast, fruit fly and primates. Gene, 105(1), 61–72.
# https://www.embl.de/pepcore/pepcore_services/cloning/choice_expression_systems/codons8

rare_codons = {
    "E. coli": ["AGG", "AGA", "ATA", "CTA", "CGA", "CGG", "CCC", "TCG"],
    "S. cerevisiae": ["AGG", "CGA", "CGG", "CGC", "CCG", "CTC", "GCG", "ACG"],
    "D. melanogaster": ["AGA", "ATA", "CGA", "CGG", "TTA", "GGG", "AGT", "TGT"],
    "Primates": ["CGA", "CGG", "TCG", "CGC", "CCG", "GCG", "ACG", "CGT"],
}


def rc(sequence: str):
    """returns the reverse complement of sequence (string)
    accepts mixed DNA/RNA
    """
    return sequence.translate(_complement_table)[::-1]


def complement(sequence: str):
    """returns the complement of sequence (string)
    accepts mixed DNA/RNA
    """
    return sequence.translate(_complement_table)


def memorize(filename):
    """Decorator for caching fucntions and classes, see pydna.download and pydna.Assembly for use"""

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


def eq(*args, **kwargs):
    """Compares two or more DNA sequences for equality i.e. they
    represent the same DNA molecule. Comparisons are case insensitive.

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

    * They have the same sequence (case insensitive)
    * One sequence is the reverse complement of the other (case insensitive)

    Two circular sequences are considered equal if they are circular permutations:

    1. They have the same lengt, AND
    2. One sequence or can be found in the concatenation of the other sequence with it    , OR
    3. The reverse complement can be found in the concatenation of the other sequence with itself.

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
        if kwargs["linear"] == True:
            topology = "linear"
        if kwargs["linear"] == False:
            topology = "circular"
    elif "circular" in kwargs:
        if kwargs["circular"] == True:
            topology = "circular"
        if kwargs["circular"] == False:
            topology = "linear"
    else:
        # topology keyword not set, look for topology associated to each sequence
        # otherwise raise exception
        topology = set(
            [arg.circular if hasattr(arg, "circular") else None for arg in args]
        )

        if len(topology) != 1:
            raise ValueError("sequences have different topologies")
        topology = topology.pop()
        if topology in (False, None):
            topology = "linear"
        elif topology == True:
            topology = "circular"

    # args_string_list    = [str(arg.seq).lower() if hasattr(arg,"seq") else str(arg).lower() for arg in args]

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


def SmallestRotation(s):
    prev, rep = None, 0
    ds = 2 * s
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
            if len(w) * rep == len(s):
                return w * rep


# try:
#    import pyximport
# except ImportError:
#    pass
# else:
#    pyximport.install()
#    from pydna._smallest import SmallestRotation


def identifier_from_string(s: str) -> str:
    """This function returns a string that is a valid python identifier
    based on the argument s or an empty string"""
    s = s.strip()
    s = _re.sub(r"\s+", r"_", s)
    s.replace("-", "_")
    s = _re.sub("[^0-9a-zA-Z_]", "", s)
    if s and not s[0].isidentifier() or _keyword.iskeyword(s):
        s = "_{s}".format(s=s)
    assert s == "" or s.isidentifier()
    return s


def seguid(seq: str) -> _pretty_str:
    """Returns the url safe SEGUID checksum for the sequence.
    This is the SEGUID checksum with the '+' and '/' characters of standard
    Base64 encoding are respectively replaced by '-' and '_'.

    Examples
    --------
    >>> from pydna.utils import seguid
    >>> seguid("a")
    'bc1M4j2I4u6VaLpUbAB8Y9kTHBs'
    """
    return (
        _pretty_str(_base64_seguid(seq.upper())
        .replace("+", "-")
        .replace("/", "_"))
    )


def lseguid(seq: str) -> _pretty_str:
    """Returns the url safe lSEGUID checksum for the sequence (seq).
    This is the SEGUID checksum with the '+' and '/' characters of standard
    Base64 encoding are respectively replaced by '-' and '_'.

    Examples
    --------
    >>> from pydna.utils import lseguid
    >>> lseguid("a")
    'bc1M4j2I4u6VaLpUbAB8Y9kTHBs'
    >>> lseguid("t")
    'bc1M4j2I4u6VaLpUbAB8Y9kTHBs'
    """
    return (
        seguid(min(seq.upper(), str(rc(seq)).upper()))
        .replace("+", "-")
        .replace("/", "_")
    )


def cseguid(seq: str) -> _pretty_str:
    """Returns the url safe cSEGUID for the sequence.
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
    return seguid(
        min(SmallestRotation(seq.upper()),
            SmallestRotation(str(rc(seq)).upper()))
    )


def flatten(*args):  # flatten
    """Flattens an iterable of iterables down to str, bytes, bytearray or any of the pydna or Biopython seq objects"""
    output = []
    args = list(args)
    while args:
        top = args.pop()
        if (
            isinstance(top, _collections.abc.Iterable)
            and not isinstance(top, (str, bytes, bytearray))
            and not hasattr(top, "features")
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
    for "Sel" and O for "Pyl") plus "Ter" for a terminator given as an asterisk.

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


# def guess_alphabet(sequence: str):
#     """
#     This function guesses the alphabet of a string representing a
#     biological sequence.

#     """

#     import string

#     from Bio.Alphabet import SingleLetterAlphabet
#     from Bio.Alphabet import NucleotideAlphabet
#     from Bio.Alphabet import ProteinAlphabet
#     from Bio.Alphabet.IUPAC import extended_protein
#     from Bio.Alphabet.IUPAC import protein
#     from Bio.Alphabet.IUPAC import ambiguous_dna
#     from Bio.Alphabet.IUPAC import unambiguous_dna
#     from Bio.Alphabet.IUPAC import extended_dna
#     from Bio.Alphabet.IUPAC import ambiguous_rna
#     from Bio.Alphabet.IUPAC import unambiguous_rna

#     if len(sequence) < 1:
#         return SingleLetterAlphabet()

#     for c in sequence:
#         if c not in string.printable:
#             return SingleLetterAlphabet()

#     xp = set(extended_protein.letters)
#     pr = set(protein.letters)

#     ad = set(ambiguous_dna.letters)
#     ud = set(unambiguous_dna.letters)
#     ed = set(extended_dna.letters)

#     ar = set(ambiguous_rna.letters)
#     ur = set(unambiguous_rna.letters)

#     all = xp | pr | ad | ud | ed | ar | ur

#     sequence_chars = set(sequence.upper())

#     if sequence_chars - all - set(string.punctuation + string.whitespace):
#         return SingleLetterAlphabet()

#     nucleic_count = 0

#     for letter in "GATCUNgatcun":
#         nucleic_count += sequence.count(letter)

#     if float(nucleic_count) / float(len(sequence)) >= 0.9:  # DNA or RNA
#         if "T" in sequence_chars and "U" in sequence_chars:
#             alphabet = NucleotideAlphabet()
#         elif not sequence_chars - ud:
#             alphabet = unambiguous_dna
#         elif not sequence_chars - ad:
#             alphabet = ambiguous_dna
#         elif not sequence_chars - ed:
#             alphabet = extended_dna
#         elif not sequence_chars - ur:
#             alphabet = unambiguous_rna
#         elif not sequence_chars - ar:
#             alphabet = ambiguous_rna
#         else:
#             alphabet = NucleotideAlphabet()
#     else:
#         threecode = [
#             "ALA",
#             "ASX",
#             "CYS",
#             "ASP",
#             "GLU",
#             "PHE",
#             "GLY",
#             "HIS",
#             "ILE",
#             "LYS",
#             "LEU",
#             "MET",
#             "ASN",
#             "PRO",
#             "GLN",
#             "ARG",
#             "SER",
#             "THR",
#             "VAL",
#             "TRP",
#             "TYR",
#             "GLX",
#             "XAA",
#             "TER",
#             "SEL",
#             "PYL",
#             "XLE",
#         ]
#         tc = set(threecode)
#         three_letter_alphabet = set(
#             [sequence[i: i + 3] for i in range(0, len(sequence), 3)]
#         )
#         if not three_letter_alphabet - tc:
#             alphabet = "three letter code"
#         elif sequence_chars - pr:
#             alphabet = protein
#         elif sequence_chars - xp:
#             alphabet = extended_protein
#         else:
#             alphabet = ProteinAlphabet()
#     return alphabet


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

    formatted = "\n".join(" ".join(cell) for cell in list_of_lists_rc)

    columnsplit = "\n|||\n".join(
        [
            "\n".join([x.strip() for x in ["".join(c) for c in zip(*col.strip("\n").splitlines())]])
            for col in cols
        ]
    )

    rowsplit = "\n---\n".join(["\n".join(a).strip() for a in zip(*list_of_lists_cr)])
    rowsplit ="\n".join(row.strip() for row in rowsplit.splitlines())

    return (_pretty_str(item) for item in ( formatted,
                                            columnsplit,
                                            rowsplit,
                                            list_of_lists_rc,
                                            list_of_lists_cr))


def join_list_to_table(rawlist):

    if "|||\n" in rawlist:
        raw_columns = rawlist.split("|||\n")
        cols = [col.splitlines() for col in raw_columns]
        if "" in [item for sublist in cols for item in sublist]:
            cols = [col.split("\n\n") for col in raw_columns]
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
        rows = [row.strip() for row in col]
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

    resultlist = []
    for line in re.finditer("(?P<item>[^\(\)]*?)(?P<brack>\[.*?\])", content):
        text2rep = line.group("item")
        bracket = line.group("brack")
        padding = max(
            [len(str(x).strip()) for x in re.split("\.\.|,", bracket.strip("[ ]"))]
        )
        inbracket = [item.strip("[ ]") for item in bracket.split(",")]
        expanded = []

        for item in inbracket:
            if re.match("(\d+\.\.\d+)|([a-z]+\.\.[a-z]+)|([A-Z]+\.\.[A-Z]+)", item):
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
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength)))
    return "".join([random.choice("GAUC") for x in range(length)])


def randomDNA(length, maxlength=None):
    """ string! """
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength)))
    return "".join([random.choice("GATC") for x in range(length)])


def randomORF(length, maxlength=None):
    length-=2
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength-2)))

    cdns = (
        "TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC",
        "TGT", "TGC", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA",
        "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT",
        "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA",
        "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT",
        "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA",
        "GGG")

    starts = ("ATG",)
    stops = ("TAA", "TAG", "TGA")

    return (
        random.choice(starts)
        + "".join([random.choice(cdns) for x in range(length)])
        + random.choice(stops)
    )


def randomprot(length, maxlength=None):
    if maxlength and maxlength > length:
        length = int(round(random.triangular(length, maxlength)))
    return "".join([random.choice("ACDEFGHIKLMNPQRSTVWY") for x in range(length)])


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
