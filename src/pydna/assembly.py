#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""Assembly of sequences by homologous recombination.

Should also be useful for related techniques such as Gibson assembly and fusion
PCR. Given a list of sequences (Dseqrecords), all sequences are analyzed for
shared homology longer than the set limit.

A graph is constructed where each overlapping region form a node and
sequences separating the overlapping regions form edges.


::


                -- A --
    catgatctacgtatcgtgt     -- B --
                atcgtgtactgtcatattc
                            catattcaaagttct



    --x--> A --y--> B --z-->   (Graph)

    Nodes:

    A : atcgtgt
    B : catattc

    Edges:

    x : catgatctacgt
    y : actgt
    z : aaagttct


The NetworkX package is used to trace linear and circular paths through the
graph.
"""
from Bio.SeqFeature import ExactPosition as _ExactPosition
from Bio.SeqFeature import FeatureLocation as _FeatureLocation
from Bio.SeqFeature import CompoundLocation as _CompoundLocation
from pydna.utils import rc as _rc
from pydna.utils import memorize as _memorize
from pydna._pretty import pretty_str as _pretty_str
from pydna.contig import Contig as _Contig
from pydna.common_sub_strings import common_sub_strings
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
import networkx as _nx
from copy import deepcopy as _deepcopy
import itertools as _itertools
import logging as _logging

_module_logger = _logging.getLogger("pydna." + __name__)


# TODO use quicker inits for contig
# TODO remove maxnodes for init


class _Memoize(type):
    @_memorize("pydna.assembly.Assembly")
    def __call__(cls, *args, **kwargs):
        return super().__call__(*args, **kwargs)


class Assembly(object, metaclass=_Memoize):
    """Assembly of a list of linear DNA fragments into linear or circular
    constructs. The Assembly is meant to replace the Assembly method as it
    is easier to use. Accepts a list of Dseqrecords (source fragments) to
    initiate an Assembly object. Several methods are available for analysis
    of overlapping sequences, graph construction and assembly.

    Parameters
    ----------

    fragments : list
        a list of Dseqrecord objects.
    limit : int, optional
        The shortest shared homology to be considered
    algorithm : function, optional
        The algorithm used to determine the shared sequences.
    max_nodes : int
        The maximum number of nodes in the graph. This can be tweaked to
        manage sequences with a high number of shared sub sequences.



    Examples
    --------

    >>> from pydna.assembly import Assembly
    >>> from pydna.dseqrecord import Dseqrecord
    >>> a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta")
    >>> b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc")
    >>> c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg")
    >>> x = Assembly((a,b,c), limit=14)
    >>> x
    Assembly
    fragments....: 33bp 34bp 35bp
    limit(bp)....: 14
    G.nodes......: 6
    algorithm....: common_sub_strings
    >>> x.assemble_circular()
    [Contig(o59), Contig(o59)]
    >>> x.assemble_circular()[0].seq.watson
    'acgatgctatactgCCCCCtgtgctgtgctctaTTTTTtattctggctgtatcGGGGGt'


    """

    def __init__(self, frags=None, limit=25, algorithm=common_sub_strings):

        # Fragments is a string subclass with some extra properties
        # The order of the fragments has significance
        fragments = []
        for f in frags:
            fragments.append(
                {
                    "upper": str(f.seq).upper(),
                    "mixed": str(f.seq),
                    "name": f.name,
                    "features": f.features,
                    "nodes": [],
                }
            )

        # rcfragments is a dict with fragments as keys and the reverse
        # complement as value
        rcfragments = dict(
            (
                f["mixed"],
                {
                    "upper": str(frc.seq).upper(),
                    "mixed": str(frc.seq),
                    "name": frc.name,
                    "features": frc.features,
                    "nodes": [],
                },
            )
            for f, frc in zip(fragments, (f.rc() for f in frags))
        )
        # The nodemap dict holds nodes and their reverse complements
        nodemap = {
            "begin": "end",
            "end": "begin",
            "begin_rc": "end_rc",
            "end_rc": "begin_rc",
        }

        # all combinations of fragments are compared.
        # see https://docs.python.org/3.10/library/itertools.html
        # itertools.combinations('ABCD', 2)-->  AB AC AD BC BD CD
        for first, secnd in _itertools.combinations(fragments, 2):

            if first["upper"] == secnd["upper"]:
                continue

            firrc = rcfragments[first["mixed"]]
            secrc = rcfragments[secnd["mixed"]]

            # matches is a list of tuples of three integers describing
            # overlapping sequences:
            # (start position in first, start position in secnd, length)
            # This comparison is done using uppercase strings, see _
            # Fragment class
            matches = algorithm(first["upper"], secnd["upper"], limit)

            for start_in_first, start_in_secnd, length in matches:
                # node is a string and represent the shared sequence in upper
                # case.
                node = first["upper"][start_in_first : start_in_first + length]

                first["nodes"].append((start_in_first, length, node))
                secnd["nodes"].append((start_in_secnd, length, node))

                # The same node exists between the reverse complements of
                # first and secnd
                # The new positions are calculated from the length of the
                # fragment and
                # the overlapping sequence
                start_in_firrc = len(first["upper"]) - start_in_first - length
                start_in_secrc = len(secnd["upper"]) - start_in_secnd - length
                # noderc is the reverse complement of node
                noderc = firrc["upper"][start_in_firrc : start_in_firrc + length]
                firrc["nodes"].append((start_in_firrc, length, noderc))
                secrc["nodes"].append((start_in_secrc, length, noderc))
                nodemap[node] = noderc

            # first is also compared to the rc of secnd
            matches = algorithm(first["upper"], secrc["upper"], limit)

            for start_in_first, start_in_secrc, length in matches:
                node = first["upper"][start_in_first : start_in_first + length]
                first["nodes"].append((start_in_first, length, node))
                secrc["nodes"].append((start_in_secrc, length, node))

                start_in_firrc, start_in_secnd = (
                    len(first["upper"]) - start_in_first - length,
                    len(secnd["upper"]) - start_in_secrc - length,
                )
                noderc = firrc["upper"][start_in_firrc : start_in_firrc + length]
                firrc["nodes"].append((start_in_firrc, length, noderc))
                secnd["nodes"].append((start_in_secnd, length, noderc))
                nodemap[node] = noderc

        # A directed graph class that can store multiedges.
        # Multiedges are multiple edges between two nodes. Each edge can hold
        # optional data or attributes.
        # https://networkx.github.io/documentation/stable/reference/classes/
        # multidigraph.html

        order = 0
        G = _nx.MultiDiGraph()
        # loop through all fragments their and reverse complements

        for f in fragments:
            f["nodes"] = sorted(set(f["nodes"]))

        for f in rcfragments.values():
            f["nodes"] = sorted(set(f["nodes"]))

        for f in _itertools.chain(fragments, rcfragments.values()):

            # nodes are sorted in place in the order of their position
            # duplicates are removed (same position and sequence)
            # along the fragment since nodes are a tuple (position(int),
            # sequence(str))

            before = G.order()
            G.add_nodes_from(
                (node, {"order": order + od, "length": length})
                for od, (start, length, node) in enumerate(
                    n for n in f["nodes"] if n[2] not in G
                )
            )
            order += G.order() - before

            for (start1, length1, node1), (
                start2,
                length2,
                node2,
            ) in _itertools.combinations(f["nodes"], 2):

                feats = [
                    ft
                    for ft in f["features"]
                    if start1 <= ft.location.start
                    and start2 + G.nodes[node2]["length"] >= ft.location.end
                ]

                for feat in feats:
                    feat.location += -start1

                G.add_edge(
                    node1,
                    node2,  # nodes (strings)
                    piece=slice(start1, start2),  # slice
                    features=feats,  # features
                    seq=f["mixed"],  # mixed case string
                    name=f["name"],
                )  # string

        self.G = _nx.create_empty_copy(G)
        self.G.add_edges_from(
            sorted(
                G.edges(data=True), key=lambda t: len(t[2].get("seq", 1)), reverse=True
            )
        )
        self.nodemap = {**nodemap, **{nodemap[i]: i for i in nodemap}}
        self.limit = limit
        self.fragments = fragments
        self.rcfragments = rcfragments
        self.algorithm = algorithm

    def assemble_linear(self, start=None, end=None, max_nodes=None):

        G = _nx.MultiDiGraph(self.G)

        G.add_nodes_from(["begin", "begin_rc", "end", "end_rc"], length=0)

        # add edges from "begin" to nodes in the first
        # sequence in self.fragments
        firstfragment = self.fragments[0]
        for start, length, node in firstfragment["nodes"][::-1]:
            G.add_edge(
                "begin",
                node,
                piece=slice(0, start),
                features=[
                    f
                    for f in firstfragment["features"]
                    if start + length >= f.location.end
                ],
                seq=firstfragment["mixed"],
                name=firstfragment["name"],
            )

        # add edges from "begin_rc" to nodes in the reverse complement of the
        # first sequence
        firstfragmentrc = self.rcfragments[firstfragment["mixed"]]
        for start, length, node in firstfragmentrc["nodes"][::-1]:
            G.add_edge(
                "begin_rc",
                node,
                piece=slice(0, start),
                features=[
                    f
                    for f in firstfragmentrc["features"]
                    if start + length >= f.location.end
                ],
                seq=firstfragmentrc["mixed"],
                name=firstfragmentrc["name"],
            )

        # add edges from nodes in last sequence to "end"
        lastfragment = self.fragments[-1]
        for start, length, node in lastfragment["nodes"]:
            G.add_edge(
                node,
                "end",
                piece=slice(start, len(lastfragment["mixed"])),
                features=[
                    f for f in lastfragment["features"] if start <= f.location.end
                ],
                seq=lastfragment["mixed"],
                name=lastfragment["name"],
            )

        # add edges from nodes in last reverse complement sequence to "end_rc"
        lastfragmentrc = self.rcfragments[lastfragment["mixed"]]
        for start, length, node in lastfragmentrc["nodes"]:
            G.add_edge(
                node,
                "end_rc",
                piece=slice(start, len(lastfragmentrc["mixed"])),
                features=[
                    f for f in lastfragmentrc["features"] if start <= f.location.end
                ],
                seq=lastfragmentrc["mixed"],
                name=lastfragmentrc["name"],
            )

        max_nodes = max_nodes or len(self.fragments)

        linearpaths = list(
            _itertools.chain(
                _nx.all_simple_paths(_nx.DiGraph(G), "begin", "end", cutoff=max_nodes),
                _nx.all_simple_paths(
                    _nx.DiGraph(G), "begin", "end_rc", cutoff=max_nodes
                ),
                _nx.all_simple_paths(
                    _nx.DiGraph(G), "begin_rc", "end", cutoff=max_nodes
                ),
                _nx.all_simple_paths(
                    _nx.DiGraph(G), "begin_rc", "end_rc", cutoff=max_nodes
                ),
            )
        )

        lps = {}

        for lp in linearpaths:
            edgelol = []

            for u, v in zip(lp, lp[1:]):
                e = []
                for d in G[u][v].values():
                    e.append((u, v, d))
                edgelol.append(e)

            for edges in _itertools.product(*edgelol):
                # TODO explain
                if [
                    True
                    for ((u, v, e), (x, y, z)) in zip(edges, edges[1:])
                    if ((e["seq"], e["piece"].stop) == (z["seq"], z["piece"].start))
                ]:
                    continue
                ct = "".join(e["seq"][e["piece"]] for u, v, e in edges)
                key = ct.upper()

                if key in lps:
                    continue    # TODO: is this test needed?
                sg = _nx.DiGraph()
                sg.add_edges_from(edges)
                sg.add_nodes_from((n, d) for n, d in G.nodes(data=True) if n in lp)

                edgefeatures = []
                offset = 0
                for u, v, e in edges:
                    feats = _deepcopy(e["features"])
                    for f in feats:
                        f.location += offset - e["piece"].start
                    edgefeatures.extend(feats)
                    offset += e["piece"].stop - e["piece"].start

                lps[key] = ct, edgefeatures, sg, {n: self.nodemap[n] for n in lp}

        return sorted(
            (
                _Contig.from_string(
                    lp[0],
                    features=lp[1],
                    graph=lp[2],
                    nodemap=lp[3],
                    linear=True,
                    circular=False,
                )
                for lp in lps.values()
            ),
            key=len,
            reverse=True,
        )

    def assemble_circular(self):
        cps = {}  # circular assembly
        cpsrc = {}
        cpaths = sorted(_nx.simple_cycles(self.G), key=len)
        cpaths_sorted = []
        for cpath in cpaths:
            order, node = min((self.G.nodes[node]["order"], node) for node in cpath)
            i = cpath.index(node)
            cpaths_sorted.append((order, cpath[i:] + cpath[:i]))
        cpaths_sorted.sort()

        for (
            _,
            cp,
        ) in (
            cpaths_sorted
        ):  # cpaths is a list of nodes representing a circular assembly
            edgelol = []  # edgelol is a list of lists of all edges along cp
            cp += cp[0:1]
            for u, v in zip(cp, cp[1:]):
                e = []
                for d in self.G[u][v].values():
                    e.append((u, v, d))
                edgelol.append(e)

            for edges in _itertools.product(*edgelol):
                if [
                    True
                    for ((u, v, e), (x, y, z)) in zip(edges, edges[1:])
                    if ((e["seq"], e["piece"].stop) == (z["seq"], z["piece"].start))
                ]:
                    continue
                ct = "".join(e["seq"][e["piece"]] for u, v, e in edges)
                key = ct.upper()

                if key in cps or key in cpsrc:
                    continue  # TODO: is test in cpsrc needed?
                sg = _nx.DiGraph()
                sg.add_edges_from(edges)
                sg.add_nodes_from((n, d) for n, d in self.G.nodes(data=True) if n in cp)

                edgefeatures = []
                offset = 0

                for u, v, e in edges:
                    feats = _deepcopy(e["features"])
                    for feat in feats:
                        feat.location += offset
                    edgefeatures.extend(feats)
                    offset += e["piece"].stop - e["piece"].start
                    for f in edgefeatures:
                        if f.location.start > len(ct) and f.location.end > len(ct):
                            f.location += -len(ct)
                        elif f.location.end > len(ct):
                            f.location = _CompoundLocation(
                                (
                                    _FeatureLocation(
                                        f.location.start, _ExactPosition(len(ct))
                                    ),
                                    _FeatureLocation(
                                        _ExactPosition(0), f.location.end - len(ct)
                                    ),
                                )
                            )

                cps[key] = cpsrc[_rc(key)] = (
                    ct,
                    edgefeatures,
                    sg,
                    {n: self.nodemap[n] for n in cp[:-1]},
                    cp,
                )

        return sorted(
            (
                _Contig.from_string(
                    cp[0],
                    features=cp[1],
                    graph=cp[2],
                    nodemap=cp[3],
                    linear=False,
                    circular=True,
                )
                for cp in cps.values()
            ),
            key=len,
            reverse=True,
        )

    def __repr__(self):
        # https://pyformat.info
        return _pretty_str(
            "Assembly\n"
            "fragments..: {sequences}\n"
            "limit(bp)..: {limit}\n"
            "G.nodes....: {nodes}\n"
            "algorithm..: {al}".format(
                sequences=" ".join(
                    "{}bp".format(len(x["mixed"])) for x in self.fragments
                ),
                limit=self.limit,
                nodes=self.G.order(),
                al=self.algorithm.__name__,
            )
        )


example_fragments = (
    _Dseqrecord("AacgatCAtgctcc", name="a"),
    _Dseqrecord("TtgctccTAAattctgc", name="b"),
    _Dseqrecord("CattctgcGAGGacgatG", name="c"),
)


linear_results = (
    _Dseqrecord("AacgatCAtgctccTAAattctgcGAGGacgatG", name="abc"),
    _Dseqrecord("ggagcaTGatcgtCCTCgcagaatG", name="ac_rc"),
    _Dseqrecord("AacgatG", name="ac"),
)


circular_results = (
    _Dseqrecord("acgatCAtgctccTAAattctgcGAGG", name="abc", circular=True),
    _Dseqrecord("ggagcaTGatcgtCCTCgcagaatTTA", name="abc_rc", circular=True),
)


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
