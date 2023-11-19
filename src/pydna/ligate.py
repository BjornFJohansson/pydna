#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.


from operator import add
from functools import reduce
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
import networkx as _nx
from itertools import permutations
import logging as _logging

_module_logger = _logging.getLogger("pydna." + __name__)


def ligate(fragments: list):
    """docstring."""
    G = _nx.DiGraph()
    G.add_nodes_from(["begin", "end"])

    for node in fragments:
        G.add_edge("begin", node)
        G.add_edge(node, "end")

    for seq1, seq2 in permutations(fragments, 2):
        try:
            seq1 + seq2
        except TypeError as err:
            if str(err) != "sticky ends not compatible!":
                raise
        else:
            if seq1.seq.three_prime_end() != (
                "blunt",
                "",
            ) and seq2.seq.five_prime_end() != ("blunt", ""):
                G.add_edge(seq1, seq2)
                G.remove_edge("begin", seq2)
                G.remove_edge(seq1, "end")

    cpaths = sorted(_nx.simple_cycles(G), key=len, reverse=True)

    csequences = [reduce(add, x).looped() for x in cpaths]

    lpaths = sorted(_nx.all_simple_paths(G, "begin", "end"), key=len, reverse=True)

    lsequences = [reduce(add, lp[1:-1]) for lp in lpaths]

    return csequences, lsequences


if __name__ == "__main__":
    a = Dseqrecord(
        Dseq.from_representation(
            """
                                            GATCaaa
                                                tttTTCC"""
        )
    )
    b = Dseqrecord(
        Dseq.from_representation(
            """
                                            AAGGanna
                                                ttntAGGA"""
        )
    )
    c = Dseqrecord(
        Dseq.from_representation(
            """
                                            TCCTcnnnn
                                                gnnnnCTAG"""
        )
    )

    d = Dseqrecord(
        Dseq.from_representation(
            """
                                               Tcnnnn
                                                gnnnnC"""
        )
    )

    e = Dseqrecord(
        Dseq.from_representation(
            """
                                                Gcnnnn
                                                 gnnnn"""
        )
    )

    fragments = [a, b, c, d, e]

    csequences, lsequences = ligate(fragments)

    for s in csequences + lsequences:
        print(repr(s.seq))
        print()
