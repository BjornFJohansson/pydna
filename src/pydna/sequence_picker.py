#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

from pydna.dseqrecord import Dseqrecord
import os as _os
import logging as _logging
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

_module_logger = _logging.getLogger("pydna." + __name__)


email = _os.getenv("pydna_email")
tool = "pydna"


def genbank_accession(s: str):
    """docstring."""
    s = Dseqrecord(s)

    NCBIWWW.email = email
    NCBIWWW.tool = tool

    result_handle = NCBIWWW.qblast(
        "blastn",
        "nt",
        str(s.seq),
        hitlist_size=1,
        alignments=1,
        descriptions=1,
        expect=1e-8,
        megablast=True,
        service="megablast",
        ungapped_alignment=True,
    )

    blast_records = NCBIXML.read(result_handle)
    best_alignment, *rest = blast_records.alignments
    best_hsp, *rest = best_alignment.hsps
    dbs = best_hsp.sbjct
    start, stop = sorted((best_hsp.sbjct_start, best_hsp.sbjct_end))
    result = Dseqrecord(
        dbs,
        circular=False,
        id=s.name,
        name=s.name,
        description=(f"{best_alignment.accession} " f"REGION: {start}..{stop}"),
    )
    return result


if __name__ == "__main__":
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
    pass
