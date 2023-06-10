#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
# doctest: +NORMALIZE_WHITESPACE
# doctest: +SKIP
"""This module provide the :class:`Anneal` class and the :func:`pcr` function
for PCR simulation. The pcr function is simpler to use, but expects only one
PCR product. The Anneal class should be used if more flexibility is required.

Primers with 5' tails as well as inverse PCR on circular templates are handled
correctly."""

from pydna._pretty import pretty_str as _pretty_str
from pydna.utils import flatten as _flatten

# from pydna.utils import memorize as _memorize
from pydna.utils import rc as _rc
from pydna.amplicon import Amplicon as _Amplicon
from pydna.primer import Primer as _Primer
from pydna.seqrecord import SeqRecord as _SeqRecord
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from Bio.SeqFeature import SeqFeature as _SeqFeature
from Bio.SeqFeature import SimpleLocation as _SimpleLocation
from Bio.SeqFeature import CompoundLocation as _CompoundLocation
from pydna.seq import Seq as _Seq
import itertools as _itertools
import re as _re
import copy as _copy
import operator as _operator
import os as _os
import logging as _logging

_module_logger = _logging.getLogger("pydna." + __name__)

_table = {  # IUPAC Ambiguity Codes for Nucleotide Degeneracy
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "(A|G)",
    "Y": "(C|T)",
    "S": "(G|C)",
    "W": "(A|T)",
    "K": "(G|T)",
    "M": "(A|C)",
    "B": "(C|G|T)",
    "D": "(A|G|T)",
    "H": "(A|C|T)",
    "V": "(A|C|G)",
    "N": "(A|G|C|T)",
}


def _annealing_positions(primer, template, limit=15):
    """Finds the annealing position(s) for a primer on a template where the
    primer anneals perfectly with at least limit nucleotides in the 3' part.
    The primer is the lower strand in the figure below.

    start is a position (integer)

    footprint and tail are strings.

    ::

        <- - - - - - - - - - template - - - - - - - - - - - - - >

        <------- start (int) ------>
     5'-...gctactacacacgtactgactgcctccaagatagagtcagtaaccacactcgat...3'
           ||||||||||||||||||||||||||||||||||||||||||||||||
                                  3'-gttctatctcagtcattggtgtATAGTG-5'

                                     <-footprint length -->

    Parameters
    ----------
    primer : string
        The primer sequence 5'-3'

    template : string
        The template sequence 5'-3'

    limit : int = 15, optional
        footprint needs to be at least of length limit.

    Returns
    -------
    describe : list of tuples (int, int)
        [ (start1, footprint1), (start2, footprint2) ,..., ]
    """

    # return empty list if primer too short
    if len(primer) < limit:
        return []

    prc = _rc(primer)

    # head is minimum part of primer that must anneal
    head = prc[:limit].upper()

    # Make regex pattern that reflects extended IUPAC DNA code
    head = "".join(_table[key] for key in head)

    positions = [m.start() for m in _re.finditer(f"(?={head})", template, _re.I)]

    if positions:
        tail = prc[limit:].lower()
        length = len(tail)
        results = []
        for match_start in positions:
            tm = template[match_start + limit : match_start + limit + length].lower()
            footprint = len(list(_itertools.takewhile(lambda x: x[0] == x[1], zip(tail, tm))))
            results.append((match_start, footprint + limit))
        return results
    return []


# class _Memoize(type):
#     @_memorize("pydna.amplify.Anneal")
#     def __call__(cls, *args, **kwargs):
#         return super().__call__(*args, **kwargs)


class Anneal(object):  # ), metaclass=_Memoize):
    """The Anneal class has the following important attributes:

    Attributes
    ----------
    forward_primers : list
        Description of `forward_primers`.
    reverse_primers : list
        Description of `reverse_primers`.
    template : Dseqrecord
        A copy of the template argument. Primers annealing sites has been
        added as features that can be visualized in a seqence editor such as
        ApE.
    limit : int, optional
        The limit of PCR primer annealing, default is 13 bp."""

    def __init__(self, primers, template, limit=13, **kwargs):
        r"""The Anneal class has to be initiated with at least an iterable of
        primers and a template.



        Parameters
        ----------
        primers : iterable of :class:`Primer` or Biopython SeqRecord like
                  objects Primer sequences 5'-3'.

        template : Dseqrecord
            The template sequence 5'-3'.

        limit : int, optional
            limit length of the annealing part of the primers.

        Attributes
        ----------
        products: list
            A list of Amplicon objects, one for each primer pair that may
            form a PCR product.


        Examples
        --------
        >>> from pydna.readers import read
        >>> from pydna.amplify import Anneal
        >>> from pydna.dseqrecord import Dseqrecord as Ds
        >>> t = Ds("tacactcaccgtctatcattatcta" +
        ...        "gatc"*240 +
        ...        "ctatcgactgtatcatctgatagcac")
        >>> from Bio.SeqRecord import SeqRecord
        >>> p1 = read(">p1\ntacactcaccgtctatcattatc", ds = False)
        >>> p2 = read(">p2\ngtgctatcagatgatacagtcg", ds = False)
        >>> ann = Anneal((p1, p2), t)
        >>> print(ann.report())
        Template name 1011 bp linear limit=13:
        p1 anneals forward (--->) at 23
        p2 anneals reverse (<---) at 989
        >>> ann.products
        [Amplicon(1011)]
        >>> amplicon_list = ann.products
        >>> amplicon = amplicon_list.pop()
        >>> amplicon
        Amplicon(1011)
        >>> print(amplicon.figure())
        5tacactcaccgtctatcattatc...cgactgtatcatctgatagcac3
                                   ||||||||||||||||||||||
                                  3gctgacatagtagactatcgtg5
        5tacactcaccgtctatcattatc3
         |||||||||||||||||||||||
        3atgtgagtggcagatagtaatag...gctgacatagtagactatcgtg5
        >>> print(amplicon)
        Dseqrecord
        circular: False
        size: 1011
        ID: 1011bp_8SyDnG-azV61tx-z8PalCWZoVDo
        Name: 1011bp_PCR_prod
        Description: pcr_product_p1_p2
        Number of features: 2
        /molecule_type=DNA
        Dseq(-1011)
        taca..gcac
        atgt..cgtg
        >>> print(amplicon.program())
        |95°C|95°C               |    |tmf:59.5
        |____|_____          72°C|72°C|tmr:59.7
        |3min|30s  \ 58.6°C _____|____|45s/kb
        |    |      \______/ 0:45|5min|GC 49%
        |    |       30s         |    |1011bp
        >>>

        """
        self.primers = primers
        self.template = _copy.deepcopy(template)

        self.limit = limit
        self.kwargs = kwargs

        self._products = None

        self.forward_primers = []
        self.reverse_primers = []

        twl = len(self.template.seq.watson)
        tcl = len(self.template.seq.crick)

        if self.template.circular:
            tw = self.template.seq.watson + self.template.seq.watson
            tc = self.template.seq.crick + self.template.seq.crick
        else:
            tw = self.template.seq.watson
            tc = self.template.seq.crick

        for p in self.primers:
            self.forward_primers.extend(
                (
                    _Primer(
                        p,
                        #          template = self.template,
                        position=tcl - pos - min(self.template.seq.ovhg, 0),
                        footprint=fp,
                    )
                    for pos, fp in _annealing_positions(str(p.seq), tc, self.limit)
                    if pos < tcl
                )
            )
            self.reverse_primers.extend(
                (
                    _Primer(
                        p,
                        #          template = self.template,
                        position=pos + max(0, self.template.seq.ovhg),
                        footprint=fp,
                    )
                    for pos, fp in _annealing_positions(str(p.seq), tw, self.limit)
                    if pos < twl
                )
            )

        self.forward_primers.sort(key=_operator.attrgetter("position"))
        self.reverse_primers.sort(key=_operator.attrgetter("position"), reverse=True)

        for fp in self.forward_primers:
            if fp.position - fp._fp >= 0:
                start = fp.position - fp._fp
                end = fp.position
                self.template.features.append(
                    _SeqFeature(
                        _SimpleLocation(start, end, strand=1),
                        type="primer_bind",
                        qualifiers={
                            "label": [fp.name],
                            "PCR_conditions": [f"primer sequence:{fp.seq}"],
                            "ApEinfo_fwdcolor": ["#baffa3"],
                            "ApEinfo_revcolor": ["#ffbaba"],
                        },
                    )
                )
            else:
                start = len(self.template) - fp._fp + fp.position
                end = start + fp._fp - len(self.template)
                sf = _SeqFeature(
                    _CompoundLocation(
                        [
                            _SimpleLocation(start, len(self.template)),
                            _SimpleLocation(0, end),
                        ]
                    ),
                    type="primer_bind",
                    qualifiers={
                        "label": [fp.name],
                        "PCR_conditions": [f"primer sequence:{fp.seq}"],
                        "ApEinfo_fwdcolor": ["#baffa3"],
                        "ApEinfo_revcolor": ["#ffbaba"],
                    },
                )
                self.template.features.append(sf)

        for rp in self.reverse_primers:
            if rp.position + rp._fp <= len(self.template):
                start = rp.position
                end = rp.position + rp._fp
                self.template.features.append(
                    _SeqFeature(
                        _SimpleLocation(start, end, strand=-1),
                        type="primer_bind",
                        qualifiers={
                            "label": [rp.name],
                            "PCR_conditions": [f"primer sequence:{rp.seq}"],
                            "ApEinfo_fwdcolor": ["#baffa3"],
                            "ApEinfo_revcolor": ["#ffbaba"],
                        },
                    )
                )
            else:
                start = rp.position
                end = rp.position + rp._fp - len(self.template)
                self.template.features.append(
                    _SeqFeature(
                        _CompoundLocation(
                            [
                                _SimpleLocation(0, end, strand=-1),
                                _SimpleLocation(start, len(self.template), strand=-1),
                            ],
                        ),
                        type="primer_bind",
                        qualifiers={"label": [rp.name]},
                    )
                )

    @property
    def products(self):
        if self._products:
            return self._products

        self._products = []

        for fp in self.forward_primers:
            for rp in self.reverse_primers:
                if self.template.circular:
                    tmpl = self.template.shifted(fp.position - fp._fp)
                    tmpl = tmpl[:] * 2
                    for f in tmpl.features:
                        for x, y in zip(f.location.parts, f.location.parts[1:]):
                            if x.end == y.start + len(self.template):
                                f.location = _SimpleLocation(
                                    x.start,
                                    y.end + len(self.template),
                                    strand=f.location.strand,
                                )
                    if fp.position > rp.position:
                        tmpl = tmpl[: len(self.template) - fp.position + rp.position + rp._fp + fp._fp]
                    else:
                        tmpl = tmpl[: rp.position + rp._fp - (fp.position - fp._fp)]
                elif fp.position <= rp.position:
                    tmpl = self.template[fp.position - fp._fp : rp.position + rp._fp]
                else:
                    continue
                prd = _Dseqrecord(fp.tail) + tmpl + _Dseqrecord(rp.tail).reverse_complement()

                full_tmpl_features = [f for f in tmpl.features if f.location.start == 0 and f.location.end == len(tmpl)]

                new_identifier = ""
                if full_tmpl_features:
                    ft = full_tmpl_features[0]
                    if "label" in ft.qualifiers:
                        new_identifier = " ".join(ft.qualifiers["label"])
                    elif "note" in ft.qualifiers:
                        new_identifier = " ".join(ft.qualifiers["note"])

                from pydna.utils import (
                    identifier_from_string as _identifier_from_string,
                )  # TODO:  clean this up

                prd.name = (
                    _identifier_from_string(new_identifier)[:16]
                    or self.kwargs.get("name")
                    or "{}bp_PCR_prod".format(len(prd))[:16]
                )
                prd.id = (
                    _identifier_from_string(new_identifier)[:16]
                    or self.kwargs.get("id")
                    or "{}bp_{}".format(str(len(prd))[:14], prd.useguid())
                )
                prd.description = self.kwargs.get("description") or "pcr_product_{}_{}".format(
                    fp.description, rp.description
                )

                amplicon = _Amplicon(
                    prd,
                    template=self.template,
                    forward_primer=fp,
                    reverse_primer=rp,
                    **self.kwargs,
                )

                # amplicon.forward_primer.amplicon = amplicon
                # amplicon.reverse_primer.amplicon = amplicon

                self._products.append(amplicon)

        return self._products

    def __repr__(self):
        """returns a short string representation"""
        return "Reaction(products = {})".format(len(self.forward_primers * len(self.reverse_primers)))

    def __str__(self):
        """returns a short report describing if or where primer
        anneal on the template."""

        mystring = "Template {name} {size} bp {top} limit={limit}:\n".format(
            name=self.template.name,
            size=len(self.template),
            top={True: "circular", False: "linear"}[self.template.circular],
            limit=self.limit,
        )
        if self.forward_primers:
            for p in self.forward_primers:
                mystring += "{name} anneals forward (--->) at {pos}\n".format(name=p.name, pos=p.position)
        else:
            mystring += "No forward primers anneal...\n"
        # mystring +="\n"
        if self.reverse_primers:
            for p in self.reverse_primers:
                mystring += "{name} anneals reverse (<---) at {pos}\n".format(name=p.name, pos=p.position)
        else:
            mystring += "No reverse primers anneal...\n"
        return _pretty_str(mystring.strip())

    report = __str__


def pcr(*args, **kwargs):
    """pcr is a convenience function for the Anneal class to simplify its
    usage, especially from the command line. If more than one or no PCR
    product is formed, a ValueError is raised.

    args is any iterable of Dseqrecords or an iterable of iterables of
    Dseqrecords. args will be greedily flattened.

    Parameters
    ----------

    args : iterable containing sequence objects
        Several arguments are also accepted.

    limit : int = 13, optional
        limit length of the annealing part of the primers.

    Notes
    -----

    sequences in args could be of type:

    * string
    * Seq
    * SeqRecord (or subclass)
    * Dseqrecord (or sublcass)

    The last sequence will be assumed to be the template while
    all preceeding sequences will be assumed to be primers.

    This is a powerful function, use with care!

    Returns
    -------

    product : Amplicon
        An :class:`pydna.amplicon.Amplicon` object representing the PCR
        product. The direction of the PCR product will be the same as for
        the template sequence.

    Examples
    --------

    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.readers import read
    >>> from pydna.amplify import pcr
    >>> from pydna.primer import Primer
    >>> template = Dseqrecord("tacactcaccgtctatcattatctac\
tatcgactgtatcatctgatagcac")
    >>> from Bio.SeqRecord import SeqRecord
    >>> p1 = Primer("tacactcaccgtctatcattatc")
    >>> p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()
    >>> pcr(p1, p2, template)
    Amplicon(51)
    >>> pcr([p1, p2], template)
    Amplicon(51)
    >>> pcr((p1,p2,), template)
    Amplicon(51)
    >>>

    """

    output = _flatten(args)  # flatten
    new = []
    for s in output:
        if hasattr(s, "watson"):
            s = _SeqRecord(_Seq(s.watson))
        elif hasattr(s, "transcribe"):
            s = _SeqRecord(s)
        elif isinstance(s, str):
            s = _SeqRecord(_Seq(s))
        elif hasattr(s, "features"):
            pass
        else:
            raise TypeError(
                "arguments need to be a string, Bio.Seq, SeqRecord" ", Primer, Dseqrecord or Amplicon object"
            )
        new.append(s)

    # A single Amplicon object
    if len(new) == 1 and hasattr(new[0], "forward_primer"):
        new = [new[0].forward_primer, new[0].reverse_primer, new[0].template]

    if not hasattr(new[-1].seq, "watson"):
        new[-1] = _Dseqrecord(s)

    anneal_primers = Anneal(new[:-1], new[-1], **kwargs)

    if len(anneal_primers.products) == 1:
        return anneal_primers.products[0]
    elif len(anneal_primers.products) == 0:
        raise ValueError(f"No PCR product! {anneal_primers.report()}")
    raise ValueError("PCR not specific! {format(anneal_primers.report()}")


if __name__ == "__main__":
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
