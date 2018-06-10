#!/usr/bin/env python3
'''This module provides functions for assembly of sequences by homologous recombination and other
related techniques. Given a list of sequences (Dseqrecords), all sequences will be analyzed for

The assembly algorithm is based on graph theory where each overlapping region forms a node and
sequences separating the overlapping regions form edges.

'''

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

class _Fragment(_Dseqrecord):
    '''This class holds information about a DNA fragment in an assembly.
    This class is instantiated by the :class:`Assembly` class and is not
    meant to be instantiated directly.
    '''
    def __init__(self, record, *args, nodes=None ,**kwargs):
        self.nodes = nodes or []
        super().__init__(record, *args, **kwargs)
        self.seq = type(self.seq)(self.seq._data)

    def rc(self):
        answer             = _deepcopy(self)
        answer._seq         = answer._seq.rc()
        answer.name        = "{}_rc".format(self.name[:13])
        answer.description = self.description+"_rc"
        answer.id          = self.id+"_rc"   
        answer.nodes       = []
        return answer

    reverse_complement = rc


class _Memoize(type):
    @_memorize("Assembly")
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

    Examples
    --------

    >>> from pydna.assembly import Assembly
    >>> from pydna.dseqrecord import Dseqrecord
    >>> a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta")
    >>> b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc")
    >>> c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg")
    >>> x = Assembly((a,b,c), limit=14)
    >>> x
    ## Assembly object ##
    fragments....: 33bp 34bp 35bp
    limit(bp)....: 14
    G.nodes......: 10
    algorithm....: common_sub_strings
    linear(3)....: -73 -54 -14
    circular(1)..: o59
    >>> x.circular
    [Contig(o59)]
    >>> x.circular_products[0].seq.watson
    'acgatgctatactgCCCCCtgtgctgtgctctaTTTTTtattctggctgtatcGGGGGt'

    '''
    
    def __init__(self, fragments, limit = 25, algorithm=common_sub_strings, max_nodes=None):
        
        ''' Consider only terminal overlaps?'''
        self.limit = limit
        ''' The shortest common sub strings to be considered '''
        self.max_nodes = max_nodes or len(fragments)
        ''' The max number of nodes allowed. This can be reset to some other value'''

        # analyze_overlaps
        
        #print(123)
        
        fragments = [_Fragment(f) for f in fragments]
        rcfragments = { f.seguid():f.rc() for f in fragments}
        g=_nx.MultiDiGraph(selfloops=False)
        
#        for first, secnd in zip(fragments, rcfragments.values()):
#                    
#                matches = algorithm( str(first.seq).upper(),
#                                     str(secnd.seq).upper(),
#                                     self.limit)
#
#                for start_in_first, start_in_secnd, length in matches:
#                    node    = first[start_in_first:start_in_first+length]
#                    node_id = node.seguid()
#                    g.add_node(node_id, length = length, fragment=str(node.seq))
#                    first.nodes.append( (start_in_first, node_id))
#                    secnd.nodes.append( (start_in_secnd, node_id))

        for first, secnd in _itertools.combinations(fragments, 2):
            
                fseguid,rseguid = first.seguid(), secnd.seguid()
                
                if fseguid==rseguid: continue

                firrc = rcfragments[fseguid]
                secrc = rcfragments[rseguid]

                matches = algorithm( str(first.seq).upper(),
                                     str(secnd.seq).upper(),
                                     self.limit)
                
                for start_in_first, start_in_secnd, length in matches:
                    node    = first[start_in_first:start_in_first+length] 
                    node_id = node.seguid()
                    g.add_node(node_id, length = length, fragment=str(node.seq))
                    first.nodes.append((start_in_first, node_id))
                    secnd.nodes.append((start_in_secnd, node_id))
                    
                    start_in_firrc, start_in_secrc = len(first.seq) - start_in_first - length, len(secnd.seq) - start_in_secnd - length # 
                    node    = firrc[start_in_firrc:start_in_firrc+length] 
                    node_id = node.seguid()
                    g.add_node(node_id, length = length, fragment=str(node.seq))
                    firrc.nodes.append((start_in_firrc, node_id) )
                    secrc.nodes.append((start_in_secrc, node_id) )

                matches = algorithm( str(first.seq).upper(),
                                     str(secrc.seq).upper(),
                                     self.limit)
                
                for start_in_first, start_in_secrc, length in matches:
                    node    = first[start_in_first:start_in_first+length] 
                    node_id = node.seguid()
                    g.add_node(node_id, length = length, fragment=str(node.seq))
                    first.nodes.append( (start_in_first, node_id))
                    secrc.nodes.append( (start_in_secrc, node_id))

                    start_in_firrc, start_in_secnd = len(first.seq) - start_in_first - length, len(secnd.seq) - start_in_secrc - length
                    node    = firrc[start_in_firrc:start_in_firrc+length]
                    node_id = node.seguid()
                    g.add_node(node_id, length = length, fragment=str(node.seq))
                    firrc.nodes.append( (start_in_firrc, node_id) )
                    secnd.nodes.append( (start_in_secnd, node_id) )

        for f in _itertools.chain(fragments, rcfragments.values()):
            f.nodes.sort()
            #print(f.nodes)
            for (s1, n1),(s2, n2) in _itertools.combinations(f.nodes,2):
                if n1==n2:
                    continue
                if not g.has_edge(n1, n2) or (str(f._seq)[s1:s2].lower() not in (e["fragment"].lower() for e in g[n1][n2].values())):
                    #print(n1, n2, g.has_edge(n1, n2))
                    #print(str(f._seq)[:10]+"|"+str(f._seq)[-10:]+str(len(str(f._seq))))
#                    try: 
#                        print([e["fragment"][:10]+"|"+e["fragment"][-10:]+str(len(e["fragment"])) for e in g[n1][n2].values()])
#                    except:
#                        pass
                    #input(21)
                    feats = [f for f in f.features if s1<=int(f.location.start) and s2+g.node[n2]["length"]>=int(f.location.end)]
                    for feat in feats:
                        feat.location+=(-s1)
                    g.add_edge(n1, n2, fragment=str(f._seq)[s1:s2], 
                                       feats=feats, 
                                       length=s2-s1, 
                                       start=s1, 
                                       end=s2, 
                                       seq = f)
        #print(777)

        # add nodes "begin", "begin_rc", "end" and "end_rc" for linear assembles
        g.add_node("begin",    length=0, fragment="")
        g.add_node("begin_rc", length=0, fragment="")
        g.add_node("end",      length=0, fragment="")
        g.add_node("end_rc",   length=0, fragment="")
       
        # add edges from "begin" to nodes in the first sequence
        for start, node_id in fragments[0].nodes:
            feats = [f for f in fragments[0].features if start+g.node[node_id]["length"]>=int(f.location.end)]
            g.add_edge("begin", node_id, fragment=str(fragments[0]._seq)[0:start], feats=feats, length=start, start=0, end=start, seq = fragments[0])

        # add edges from "begin_rc" to nodes in the first reverse complement sequence
        for start, node_id in list(rcfragments.values())[0].nodes:
            feats = [f for f in list(rcfragments.values())[0].features if start+g.node[node_id]["length"]>=int(f.location.end)]
            g.add_edge("begin_rc", node_id, fragment=str(list(rcfragments.values())[0]._seq)[0:start], feats=feats, length =start, start=0, end=start, seq = list(rcfragments.values())[0])

        # add edges from nodes in last sequence to "end"
        for start_in_last, last_node_id in fragments[-1].nodes:
            feats = [f for f in fragments[-1].features if start_in_last<=int(f.location.start)]
            g.add_edge(last_node_id, "end", fragment=str(fragments[-1]._seq)[start_in_last:len(fragments[-1])],feats=feats, length =start, start=start_in_last, end=len(fragments[-1]), seq=fragments[-1])

        # add edges from nodes in last reverse complement sequence to "end_rc"
        for start_in_last, last_node_id in list(rcfragments.values())[-1].nodes:
            feats = [f for f in list(rcfragments.values())[-1].features if start_in_last<=int(f.location.start)]
            g.add_edge(last_node_id, "end_rc", fragment=str(list(rcfragments.values())[-1]._seq)[start_in_last:len(list(rcfragments.values())[-1])], feats=feats, length =start, start=start_in_last, end=len(fragments[-1]), seq=list(rcfragments.values())[-1])
        
        lps={}     #linear assembly
               
        lpths = _itertools.chain(_nx.all_simple_paths(_nx.DiGraph(g),"begin","end",      cutoff=self.max_nodes),
                                 _nx.all_simple_paths(_nx.DiGraph(g),"begin","end_rc",   cutoff=self.max_nodes),
                                 _nx.all_simple_paths(_nx.DiGraph(g),"begin_rc","end",   cutoff=self.max_nodes),
                                 _nx.all_simple_paths(_nx.DiGraph(g),"begin_rc","end_rc",cutoff=self.max_nodes))
        #print(g['h0FkKjiojpATbfkOM0C3u4zg9hs'][ 'gmfjuQLVSPP4ayjJMPuig1jxxmE'])
        
        #print(888)

        for lpath in lpths:
            #print(lpath)
            e1=[]
            for u,v in zip(lpath,lpath[1:]):
                e2=[]
                for d in g[u][v].values():
                    e2.append((u,v,d))
                e1.append(e2)
            #print(e1)
            #input(42)
            for edges in _itertools.product(*e1):
                #print("edges")
                sg=_nx.DiGraph(g.subgraph(lpath).copy())
                sg.add_edges_from(edges)
                ct = "".join(e[2]["fragment"] for e in edges)
                lseguid = _lseguid(ct)
                if lseguid in lps:
                    continue
                edgefeatures=[]
                offset=0
                for e in edges:
                    feats = _deepcopy(e[2]["feats"])
                    for f in feats:
                        f.location+=offset
                    edgefeatures.extend(feats)
                    offset+=e[2]["length"]
                lps[lseguid]= _Contig( ct, features=edgefeatures, graph=sg, path=lpath)
        #print(111)
        cps = {} # circular assembly
        nodes  = list(_itertools.chain.from_iterable([f.nodes for f in fragments]))
        nodes  = list(dict.fromkeys([n[1] for n in nodes]))

        cpaths = [list(x) for x in _nx.simple_cycles(g)]

        first_cpaths  = []
        second_cpaths = []
        #cpaths.sort(key=len)
        
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
                for d in g[u][v].values():
                    e2.append((u,v,d))
                e1.append(e2)
            for edges in _itertools.product(*e1):
                ct = "".join(e[2]["fragment"] for e in edges)
                if ct in [str(s.seq) for s in cps.values()]:
                    continue
                cseguid = _cseguid(ct)                
                if cseguid in cps:
                    continue           
                sg=_nx.DiGraph(g.subgraph(cp).copy())
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
                cps[cseguid]=_Contig( ct, features = edgefeatures, graph=sg, path=cp, circular=True)

        self.linear   = self.linear_products   = sorted(lps.values(), key=len, reverse=True)
        self.circular = self.circular_products = sorted(cps.values(), key=len, reverse=True)
        self.fragments = fragments
        self.G = g
        self.algorithm = algorithm
        #globals().update(locals());import sys;sys.exit(42)

    def list_circular(self):
        return _pretty_str("\n".join("{i} {repr(p)} {cs}".format(i=i,p=p,cs=p.cseguid()) for i,p in enumerate(self.circular)))
        
    list_circular_products = list_circular
    
    def list_linear(self):
        return _pretty_str("\n".join("{i} {repr(p)} {ls}".format(i=i,p=p,ls=p.lseguid()) for i,p in enumerate(self.linear)))
        
    list_linear_products = list_linear
        
    def __repr__(self):
        # https://pyformat.info
        return _pretty_str( "## Assembly object ##\n"
                            "fragments....: {sequences}\n"
                            "limit(bp)....: {limit}\n"
                            "G.nodes......: {nodes}\n"
                            "algorithm....: {al}\n"
                            "linear{lenlp:.<7}: {lp}\n"    
                            "circular{lencp:.<5}: {cp}".format(sequences       = " ".join("{}bp".format(len(x)) for x in self.fragments),
                                                                        limit           = self.limit,
                                                                        nodes           = self.G.order(),
                                                                        al              = self.algorithm.__name__,
                                                                        lencp           = "({lencp})".format(lencp=len(self.circular_products)),                                                          
                                                                        cp              = " ".join("o{}".format(len(x)) for x in self.circular_products[:8]),
                                                                        lenlp           = "({lenlp})".format(lenlp=len(self.linear_products)),
                                                                        lp              = " ".join("-{}".format(len(x)) for x in self.linear_products[:8])))

if __name__=="__main__":
    from pydna.readers import read
    
    pMEC1135 = read("pMEC1135.gb")
    hygromycin_product = read("hygromycin_product.gb")
    
    asm_hyg = Assembly((pMEC1135, hygromycin_product, pMEC1135))

    