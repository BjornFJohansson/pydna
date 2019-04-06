#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''This module provides functions for assembly of sequences by homologous 
recombination and other related techniques. Given a list of sequences 
(Dseqrecords), all sequences are analyzed for shared homology longer than the 
set limit.

A graph is constructed where each overlapping region form a node and
sequences separating the overlapping regions form edges.

::
           
                 -- A --
     catgatctacgtatcgtgt     -- B --
                 atcgtgtactgtcatattc
                             catattcaaagttct
                 
     Graph:
                 --x--> A --y--> B --z-->
                 
                 Nodes:
                     
                 A : atcgtgt
                 B : catattc
                 
                 Edges:

                 x : catgatctacgt
                 y : actgt
                 z : aaagttct

The NetworkX package is used to trace linear and circular paths through the
graph.

'''
import sys
    
if sys.version_info < (3, 6):
    from collections import OrderedDict as _od
else:
    _od = dict

import logging as _logging
_module_logger = _logging.getLogger("pydna."+__name__)

import itertools as _itertools
from   copy     import deepcopy as _deepcopy
import networkx as _nx

from pydna.dseqrecord import Dseqrecord as _Dseqrecord

from pydna.common_sub_strings import common_sub_strings
from pydna.common_sub_strings import terminal_overlap
from pydna.contig  import Contig as _Contig
from pydna._pretty import pretty_str as _pretty_str
from pydna.utils   import memorize   as _memorize
from pydna.utils   import lseguid as _lseguid
from pydna.utils   import cseguid as _cseguid

from Bio.SeqFeature import CompoundLocation as _CompoundLocation
from Bio.SeqFeature import FeatureLocation  as _FeatureLocation
from Bio.SeqFeature import ExactPosition    as _ExactPosition
from collections import UserString


class _Fragment(UserString):
    '''This class holds information about a DNA fragment in an assembly.
    This class is instantiated by the :class:`Assembly` class and is not
    meant to be instantiated directly.
    '''
    def __init__(self, record, *args, nodes=None, **kwargs):
        self.data   = record.seq.todata.upper()
        self.original_case = record.seq.todata
        self.name   = record.name
        self.record = record
        self.nodes  = nodes or []
        #(upper, originalcase, name, record, nodes)

    def __add__(self, other):
        return self.data + str(other)
    
    def __radd__(self, other):
        return str(other) + self.data
    
    def __getitem__(self, index): 
        return self.data[index]
   
# class _Fragment(str):

#     def __new__(cls, record, *args, **kwargs):
#         return super(_Fragment, cls).__new__(cls, record.seq.todata.upper())

#     def __init__(self, record, nodes=None):
#         self.original_case = record.seq.todata
#         self.name   = record.name
#         self.record = record
#         self.nodes  = nodes or []

class _Memoize(type):
    @_memorize("pydna.assembly.Assembly")
    def __call__(cls, *args, **kwargs):
        return super().__call__(*args, **kwargs)


class Assembly(object, metaclass = _Memoize):
    '''Assembly of a list of linear DNA fragments into linear or circular constructs.
    The Assembly is meant to replace the Assembly method as it is easier to use.
    Accepts a list of Dseqrecords (source fragments) to initiate an Assembly object.
    Several methods are available for analysis of overlapping sequences, graph construction
    and assembly.

    Parameters
    ----------

    fragments : list
        a list of Dseqrecord objects.
    limit : int, optional
        The shortest shared homology to be considered
    algorithm : function, optional
        The algorithm used to determine the shared sequences.
    max_nodes : int
          The maximum number of nodes in the graph. This can be tweaked to manage 
        sequences with a high number of shared sub sequences.
    
    

    Examples
    --------

    >>> from pydna.assembly import Assembly
    >>> from pydna.dseqrecord import Dseqrecord
    >>> a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta")
    >>> b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc")
    >>> c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg")
    >>> x = Assembly((a,b,c), limit=14)
    >>> x
    Assembly (max_nodes=3)
    fragments....: 33bp 34bp 35bp
    limit(bp)....: 14
    G.nodes......: 6
    algorithm....: common_sub_strings
    >>> x.assemble_circular()
    [Contig(o59)]
    >>> x.assemble_circular()[0].seq.watson
    'acgatgctatactgCCCCCtgtgctgtgctctaTTTTTtattctggctgtatcGGGGGt'

    '''
    
    
    def __init__(self, frags=[],  limit = 25, algorithm=common_sub_strings, max_nodes=None, **attr):
        
        fragments   = [_Fragment(f) for f in frags]
        rcfragments = _od( (f, _Fragment(f.record.rc())) for f in fragments)
        
        
        nodemap = {"begin":"end","end":"begin"}
        G = _nx.MultiDiGraph()

        for first, secnd in _itertools.combinations(fragments, 2):

                if first == secnd:
                    continue
            
                firrc = rcfragments[first]
                secrc = rcfragments[secnd]

                matches = algorithm( first,
                                     secnd,
                                     limit)

                for start_in_first, start_in_secnd, length in matches:
                    node    = first[start_in_first:start_in_first+length]
                    G.add_node(node, length = length)
                    first.nodes.append((start_in_first, node))
                    secnd.nodes.append((start_in_secnd, node))

                    start_in_firrc, start_in_secrc = len(first) - start_in_first - length, len(secnd) - start_in_secnd - length
                    noderc  = firrc[start_in_firrc:start_in_firrc+length]
                    G.add_node(noderc, length = length)
                    firrc.nodes.append((start_in_firrc, noderc) )
                    secrc.nodes.append((start_in_secrc, noderc) )
                    nodemap[node]=noderc
                matches = algorithm( first,
                                     secrc,
                                     limit)
                
                for start_in_first, start_in_secrc, length in matches:
                    node    = first[start_in_first:start_in_first+length]
                    G.add_node(node, length = length)
                    first.nodes.append( (start_in_first, node))
                    secrc.nodes.append( (start_in_secrc, node))

                    start_in_firrc, start_in_secnd = len(first) - start_in_first - length, len(secnd) - start_in_secrc - length
                    noderc  = firrc[start_in_firrc:start_in_firrc+length]
                    G.add_node(noderc, length = length)
                    firrc.nodes.append( (start_in_firrc, noderc) )
                    secnd.nodes.append( (start_in_secnd, noderc) )
                    nodemap[node]=noderc


        for f in _itertools.chain(fragments, rcfragments.values()):            
            f.nodes.sort()

            for (s1, n1),(s2, n2) in _itertools.combinations(f.nodes,2):
                if n1==n2:
                    continue
                if not G.has_edge(n1, n2) or (f[s1:s2].lower() not in (e["fragment_original_case"].lower() for e in G[n1][n2].values())):

                    feats = [f for f in f.record.features if s1<=int(f.location.start) and s2+G.node[n2]["length"]>=int(f.location.end)]
                    
                    for feat in feats:
                        feat.location+=(-s1)
                    
                    G.add_edge(n1, n2,
                               fragment_original_case=f.original_case[s1:s2],
                               feats=feats, 
                               length=s2-s1, 
                               start=s1, 
                               end=s2, 
                               seq=f.original_case,
                               name=f.name)
        self.G=G
        self.nodemap=nodemap
        self.limit = limit
        self.fragments = fragments
        self.rcfragments = rcfragments
        self.algorithm = algorithm
        self.max_nodes = max_nodes or len(self.fragments)


    def assemble_linear(self, start=None,end=None, max_nodes=0):
        
        self.G.add_node("begin",    length=0)
        self.G.add_node("begin_rc", length=0)
        self.G.add_node("end",      length=0)
        self.G.add_node("end_rc",   length=0)
        
        frst = self.fragments[0].original_case
        frst_name = self.fragments[0].name
        last = self.fragments[-1].original_case
        last_name = self.fragments[-1].name
        frstrc = str(list(self.rcfragments.values())[0].original_case)
        frstrc_name = list(self.rcfragments.values())[0].name
        lastrc = str(list(self.rcfragments.values())[-1].original_case)
        lastrc_name = list(self.rcfragments.values())[-1].name
        # add edges from "begin" to nodes in the first sequence
        for start, node in self.fragments[0].nodes:
            feats = [f for f in self.fragments[0].record.features if start+self.G.node[node]["length"]>=int(f.location.end)]
            self.G.add_edge("begin", node, 
                            fragment_original_case=frst[0:start], 
                            feats  = feats, 
                            length = start, 
                            start  = 0, 
                            end    = start, 
                            seq    = frst,
                            name   = frst_name)

        # add edges from "begin_rc" to nodes in the first reverse complement sequence
        for start, node_id in list(self.rcfragments.values())[0].nodes:
            feats = [f for f in list(self.rcfragments.values())[0].record.features if start+self.G.node[node_id]["length"]>=int(f.location.end)]
            self.G.add_edge("begin_rc", node_id, 
                            fragment_original_case = frstrc[0:start], 
                            feats = feats, length =start, 
                            start = 0, 
                            end   = start, 
                            seq   = frstrc,
                            name  = frstrc_name)

        # add edges from nodes in last sequence to "end"
        for start_in_last, last_node_id in self.fragments[-1].nodes:
            feats = [f for f in self.fragments[-1].record.features if start_in_last<=int(f.location.start)]
            self.G.add_edge(last_node_id, "end", 
                            fragment_original_case = last[start_in_last:len(self.fragments[-1])],
                            feats  = feats, 
                            length = start, 
                            start  = start_in_last, 
                            end    = len(self.fragments[-1]), 
                            seq    = last,
                            name   = last_name)

        # add edges from nodes in last reverse complement sequence to "end_rc"
        for start_in_last, last_node_id in list(self.rcfragments.values())[-1].nodes:
            feats = [f for f in list(self.rcfragments.values())[-1].record.features if start_in_last<=int(f.location.start)]
            self.G.add_edge(last_node_id, "end_rc", 
                            fragment_original_case = lastrc[start_in_last:len(list(self.rcfragments.values())[-1])],
                            feats  = feats, 
                            length = start, 
                            start  = start_in_last, 
                            end    = len(self.fragments[-1]), 
                            seq    = lastrc,
                            name   = lastrc_name)
        
        lps={}        #linear assembly
        lps=_od(lps)  # fix
        
        lpths = _itertools.chain( _nx.all_simple_paths(_nx.DiGraph(self.G),"begin",    "end",      cutoff=self.max_nodes),
                                  _nx.all_simple_paths(_nx.DiGraph(self.G),"begin",    "end_rc",   cutoff=self.max_nodes),
                                  _nx.all_simple_paths(_nx.DiGraph(self.G),"begin_rc", "end",      cutoff=self.max_nodes),
                                  _nx.all_simple_paths(_nx.DiGraph(self.G),"begin_rc", "end_rc",   cutoff=self.max_nodes), 
                                  )

        for lpth in lpths:
            e1=[]
            for u,v in zip(lpth,lpth[1:]):
                e2=[]
                for d in self.G[u][v].values():
                    e2.append((u,v,d))
                e1.append(e2)

            for edges in _itertools.product(*e1):
                sg=_nx.DiGraph(self.G.subgraph(lpth).copy())
                sg.add_edges_from(edges)
                ct = "".join(e[2]["fragment_original_case"] for e in edges)
                if ct in lps: 
                    continue
                edgefeatures=[]
                offset=0
                for e in edges:
                    feats = _deepcopy(e[2]["feats"])
                    for f in feats:
                        f.location+=offset
                    edgefeatures.extend(feats)
                    offset+=e[2]["length"]
                lps[ct]= _Contig( ct, features=edgefeatures, graph=sg, path=lpth, nodemap=self.nodemap)

        return sorted(lps.values(), key=len, reverse=True)
    

    def assemble_circular(self):
        cps = {} # circular assembly
        cps = _od(cps)
        nodes  = list(_itertools.chain.from_iterable([f.nodes for f in self.fragments]))
        nodes  = list(_od.fromkeys([n[1] for n in nodes]))
        cpaths = [list(x) for x in _nx.simple_cycles(self.G)]
        first_cpaths  = []
        second_cpaths = []
        
        for cpath in cpaths:
            for j, node in enumerate(nodes):
                cp=[]
                try:
                    i = cpath.index(node)
                except ValueError:
                    pass
                else:
                    cp = cpath[i:]+cpath[:i]
                    first_cpaths.append((j, cp))
                if cp:
                    break
            if not cp:
                second_cpaths.append(cpath)

        first_cpaths.sort()
        first_cpaths = [cp for (j,cp) in first_cpaths]

        cpaths = first_cpaths + second_cpaths
        cpaths.sort(key=len)
        

        for cp in cpaths:
            e1=[]
            cp+=cp[0:1]
            for u,v in zip(cp, cp[1:]):
                e2=[]
                for d in self.G[u][v].values():
                    e2.append((u,v,d))
                e1.append(e2)

            for edges in _itertools.product(*e1):
                ct = "".join(e[2]["fragment_original_case"] for e in edges)
                if ct in [str(s.seq) for s in cps.values()]:
                    continue
                cseguid = _cseguid(ct)                   
                if cseguid in cps:
                    continue           
                sg=_nx.DiGraph(self.G.subgraph(cp).copy())
                sg.add_edges_from(edges)
                edgefeatures=[]
                offset=0
                for e in edges:
                    feats = _deepcopy(e[2]["feats"])
                    for feat in feats:
                        feat.location+=offset                    
                    edgefeatures.extend(feats)
                    offset+=e[2]["length"]
                for f in edgefeatures:
                    if f.location.start>len(ct) and f.location.end>len(ct):                        
                        f.location+=(-len(ct))                    
                    elif f.location.end>len(ct):
                        f.location = _CompoundLocation((_FeatureLocation(f.location.start,_ExactPosition(len(ct))),_FeatureLocation(_ExactPosition(0), f.location.end-len(ct))))
                cps[cseguid]=_Contig( ct, features = edgefeatures, graph=sg, path=cp,nodemap=self.nodemap,circular=True)
        return sorted(cps.values(), key=len, reverse=True)

        
    def __repr__(self):
        # https://pyformat.info
        return _pretty_str( "Assembly (max_nodes={max_nodes})\n"
                            "fragments....: {sequences}\n"
                            "limit(bp)....: {limit}\n"
                            "G.nodes......: {nodes}\n"
                            "algorithm....: {al}".format(sequences = " ".join("{}bp".format(len(x)) for x in self.fragments),
                                                           limit     = self.limit,
                                                           nodes     = self.G.order(),
                                                           max_nodes = self.max_nodes,    
                                                           al        = self.algorithm.__name__))


example_fragments = ( _Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta",  name ="a"),
                      _Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc", name ="b"),
                      _Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg",name ="c") )


if __name__=="__main__":
    # asm = Assembly(example_fragments, limit=14)
    # lin = asm.assemble_linear()
    # print( lin )
    # for l in lin:
    #    print(l.figure())
 
    # crc = asm.assemble_circular()
    # for c in crc:
    #     print(c.figure())
          
                         
      
    # acgatgctatactgCCCCCtgtgctgtgctcta
    #                    tgtgctgtgctctaTTTTTtattctggctgtatc
    #                                       tattctggctgtatcGGGGGtacgatgctatactg
    
    # acgatgctatactgCCCCCtgtgctgtgctctaTTTTTtattctggctgtatcGGGGGtacgatgctatactg
    
    
    # acgatgctatactgCCCCCtgtgctgtgctctaTTTTTtattctggctgtatcGGGGGt    
    
    a = _Dseqrecord("acgatgctatactggCCCCCtgtgctgtgctctaGG",name="oneC")
    a.add_feature(1,33,label="first")
                                       # tgtgctgtgctcta 14 
    b =                     _Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatct", name="twoA")
    b.add_feature(1,34,label="scnd")
    b2=                     _Dseqrecord("tgtgctgtgctctaCCtattctggctgtatct",name="twoB")
                                                         # tattctggctgtatct 16 
    c =                                        _Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three")
    c.add_feature(1,36,label="third")
                                                          
                                                                               # acgatgctatactgg 15
                                                                               
    c2 = Assembly((a,b,b2,c), limit=14)
    assert c2.assemble_circular()[0].cseguid() == "t3mIjxv3Q5GK9SWpXD-UfyefANc"
    assert c2.assemble_circular()[1].cseguid() == "k9ztaDj9HsQYZvxzvkUWn6SY5Ks"
    assert str(c2.assemble_circular()[0].seq)=='acgatgctatactggCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGT'
    assert str(c2.assemble_circular()[1].seq)=='acgatgctatactggCCCCCtgtgctgtgctctaCCtattctggctgtatctGGGGGT'
    print("done")