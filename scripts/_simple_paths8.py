#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#    Copyright (C) 2012 by
#    Sergio Nery Simoes <sergionery@gmail.com>
#    All rights reserved.
#    BSD license (see second license in LICENSE.txt).

import networkx as _nx

__author__ = """\n""".join(
    [
        "Sérgio Nery Simões <sergionery@gmail.com>",
        "Aric Hagberg <aric.hagberg@gmail.com>",
    ]
)
__all__ = ["all_circular_paths_edges", "all_circular_paths_edges"]


def _all_simple_paths_graph(G, source, target, cutoff=None):
    if source not in G:
        raise _nx.NetworkXError("source node %s not in graph" % source)
    if target not in G:
        raise _nx.NetworkXError("target node %s not in graph" % source)
    if cutoff is None:
        cutoff = len(G) - 1

    if cutoff < 1:
        return
    visited = [source]
    stack = [iter(G[source])]
    while stack:
        children = stack[-1]
        child = next(children, None)
        if child is None:
            stack.pop()
            visited.pop()
        elif len(visited) < cutoff:
            if child == target:
                yield visited + [target]
            elif child not in visited:
                visited.append(child)
                stack.append(iter(G[child]))
        else:  # len(visited) == cutoff:
            if child == target or target in children:
                yield visited + [target]
            stack.pop()
            visited.pop()


def all_simple_paths_edges(G, source, target, cutoff=None, data=False):
    if data is True:

        def edge_data(u, v, n, E, i):
            return (u, v, dict([E[n][i[n]]]))

    else:

        def edge_data(u, v, n, E, i):
            return (u, v)

    for path in _all_simple_paths_graph(G, source, target, cutoff=cutoff):
        edges = list(zip(path[:-1], path[1:]))
        E = []  # list: items of each edge
        N = []  # list: number of items of each edge
        for u, v in edges:
            edge_items = list(G[u][v].items())
            E += [edge_items]
            N += [len(edge_items)]
        i_list = [0 for n in N]
        idx = [i for i in reversed(list(range(len(i_list))))]
        while True:
            path_edges = []
            for n, (u, v) in enumerate(edges):
                path_edges += [edge_data(u, v, n, E, i_list)]
            yield path_edges
            for i in idx:
                i_list[i] = (i_list[i] + 1) % N[i]
                if i_list[i] != 0:
                    break
            if i == 0 and i_list[0] == 0:
                break


def all_circular_paths_edges(G):
    for path in sorted(_nx.simple_cycles(G), key=len, reverse=True):
        edges = list(zip(path, path[1:] + [path[0]]))
        N = []
        for u, v in edges:
            n = len(G[u][v])
            N += [n]
        i_list = [0 for n in N]
        idx = [i for i in reversed(list(range(len(i_list))))]
        while True:
            path_edges = []
            for i, (u, v) in enumerate(edges):
                path_edges += [(u, v, G[u][v][i_list[i]])]
            yield path_edges
            for i in idx:
                i_list[i] = (i_list[i] + 1) % N[i]
                if i_list[i] != 0:
                    break
            if i == 0 and i_list[0] == 0:
                break


if __name__ == "__main__":
    import os as _os

    cache = _os.getenv("pydna_cache", "nocache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"] = cache

    G = _nx.MultiDiGraph()
    G.add_edge("a", "b", weight=0.1)
    G.add_edge("a", "b", weight=1.0)
    G.add_edge("b", "c", weight=0.2)
    G.add_edge("b", "c", weight=2.0)
    G.add_edge("c", "d", weight=0.3)
    G.add_edge("c", "d", weight=3.0)
    G.add_edge("a", "c", weight=10, key="77")

    source = "a"
    target = "d"

    # MULTIDIGRAPH
    print("MULTIDIGRAPH (data=False)")
    for path in all_simple_paths_edges(G, source, target):
        print(path)

    print()
    print("MULTIDIGRAPH (data=True)")
    for path in all_simple_paths_edges(G, source, target, data=True):
        # print path
        total_weight = sum([(list(I.values())[0]["weight"]) for u, v, I in path])
        print(total_weight, "\t", path)

    print()
    # DIGRAPH
    H = _nx.DiGraph(G)
    print("DIGRAPH (data=False)")
    for path in all_simple_paths_edges(H, source, target):
        print(path)

    print()
    print("DIGRAPH (data=True)")
    for path in all_simple_paths_edges(H, source, target, data=True):
        total_weight = sum([w["weight"] for u, v, w in path])
        print(total_weight, "\t", path)
