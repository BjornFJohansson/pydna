#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides functions for assembly of sequences by homologous recombination and other
related techniques. Given a list of sequences (Dseqrecords), all sequences will be analyzed for
overlapping regions of DNA (common substrings).

The assembly algorithm is based on graph theory where each overlapping region forms a node and
sequences separating the overlapping regions form edges.

'''

import cPickle
import shelve

import logging
module_logger = logging.getLogger("pydna."+__name__)

import itertools
import networkx as nx
import operator
import random
import os

from copy import copy
from textwrap import dedent
from collections import defaultdict
from collections import namedtuple

from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature

from pydna.dsdna import Dseq
from pydna.dsdna import Dseqrecord
from pydna._simple_paths8 import all_simple_paths_edges, all_circular_paths_edges
from findsubstrings_suffix_arrays_python import common_sub_strings
from findsubstrings_suffix_arrays_python import terminal_overlap

from orderedset import OrderedSet
from pydna.pretty                   import pretty_str



class Fragment(Dseqrecord):
    '''This class holds information about a DNA fragment in an assembly.
    This class is instantiated by the :class:`Assembly` class and is not
    meant to be instantiated directly.

    '''

    def __init__(self, record, start1    = 0,
                               end1      = 0,
                               start2    = 0,
                               end2      = 0,
                               alignment = 0,
                               i         = 0, *args, **kwargs):

        super(Fragment, self).__init__(record, *args, **kwargs)

        self.start1             = start1
        self.end1               = end1
        self.left_overlap_size  = end1-start1
        self.start2             = start2
        self.end2               = end2
        self.right_overlap_size = end2-start2
        self.alignment          = alignment
        self.i                  = i

    def __str__(self):
        return ("Fragment alignment {}\n").format(self.alignment)+super(Fragment, self).__str__()

class Contig(Dseqrecord):
    '''This class holds information about a DNA assembly. This class is instantiated by
    the :class:`Assembly` class and is not meant to be instantiated directly.

    '''

    def __init__(self,
                 record,
                 source_fragments=[],
                 *args, **kwargs):

        super(Contig, self).__init__(record, *args, **kwargs)
        self.source_fragments = source_fragments
        self.number_of_fragments = len(self.source_fragments)

    def __repr__(self):
        return "Contig({}{})".format({True:"-", False:"o"}[self.linear],len(self))

    def detailed_figure(self):
        '''Synonym of :func:`detailed_fig`'''
        return self.detailed_fig()

    def detailed_fig(self):
        fig=""
        for s in self.source_fragments:
            fig +="{}{}\n".format(" "*s.alignment, str(s.seq))
        return fig

    def figure(self):
        '''Synonym of :func:`small_fig`'''
        return self.small_fig()

    def small_figure(self):
        '''Synonym of :func:`small_fig`'''
        return self.small_fig()

    def small_fig(self):
        '''
        Returns a small ascii representation of the assembled fragments. Each fragment is
        represented by:

        ::

         Size of common 5' substring|Name and size of DNA fragment| Size of common 5' substring

        Linear:

        ::

          frag20| 6
                 \\/
                 /\\
                  6|frag23| 6
                           \\/
                           /\\
                            6|frag14


        Circular:

        ::

          -|2577|61
         |       \\/
         |       /\\
         |       61|5681|98
         |               \\/
         |               /\\
         |               98|2389|557
         |                       \\/
         |                       /\\
         |                       557-
         |                          |
          --------------------------


        '''

        if self.linear:
            '''
            frag20| 6
                   \/
                   /\
                    6|frag23| 6
                             \/
                             /\
                              6|frag14
            '''
            f = self.source_fragments[0]
            space2 = len(f.name)


            fig = ("{name}|{o2:>2}\n"
                   "{space2} \/\n"
                   "{space2} /\\\n").format(name = f.name,
                                            o2 = f.right_overlap_size,
                                            space2 = " "*space2)
            space = len(f.name)

            for f in self.source_fragments[1:-1]:
                name= "{o1:>2}|{name}|".format(o1   = f.left_overlap_size,
                                               name = f.name)
                space2 = len(name)
                fig +=("{space} {name}{o2:>2}\n"
                       "{space} {space2}\/\n"
                       "{space} {space2}/\\\n").format( name = name,
                                                        o2 = f.right_overlap_size,
                                                        space = " "*space,
                                                        space2 = " "*space2)
                space +=space2
            f = self.source_fragments[-1]
            fig += ("{space} {o1:>2}|{name}").format(name = f.name,
                                                    o1 = f.left_overlap_size,
                                                    space = " "*(space))



        else:
            '''
             -|2577|61
            |       \/
            |       /\
            |       61|5681|98
            |               \/
            |               /\
            |               98|2389|557
            |                       \/
            |                       /\
            |                       557-
            |                          |
             --------------------------
            '''
            f = self.source_fragments[0]
            space = len(f.name)+3
            fig =(" -|{name}|{o2:>2}\n"
                  "|{space}\/\n"
                  "|{space}/\\\n").format(name = f.name,
                                           o2 = f.right_overlap_size,
                                           space = " "*space)
            for f in self.source_fragments[1:]:
                name= "{o1:>2}|{name}|".format(o1 = f.left_overlap_size,
                                                      name = f.name)
                space2 = len(name)
                fig +=("|{space}{name}{o2:>2}\n"
                       "|{space}{space2}\/\n"
                       "|{space}{space2}/\\\n").format(o2 = f.right_overlap_size,
                                                       name = name,
                                                       space = " "*space,
                                                       space2 = " "*space2)
                space +=space2

            fig +="|{space}{o1:>2}-\n".format(space=" "*(space), o1=self.source_fragments[0].left_overlap_size)
            fig +="|{space}   |\n".format(space=" "*(space))
            fig +=" {space}".format(space="-"*(space+3))
        return pretty_str(dedent(fig))

class Assembly(object):
    '''Assembly of a list of linear DNA fragments into linear or circular constructs.
    The Assembly is meant to replace the Assembly method as it is easier to use.
    Accepts a list of Dseqrecords (source fragments) to initiate an Assembly object.
    Several methods are available for analysis of overlapping sequences, graph construction
    and assembly.

    Parameters
    ----------

    dsrecs : list
        a list of Dseqrecord objects.

    Examples
    --------

    >>> from pydna import Assembly, Dseqrecord
    >>> a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta")
    >>> b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc")
    >>> c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg")
    >>> x = Assembly((a,b,c), limit=14)
    >>> x
    Assembly:
    Sequences........................: [33] [34] [35]
    Sequences with shared homologies.: [33] [34] [35]
    Homology limit (bp)..............: 14
    Number of overlaps...............: 3
    Nodes in graph(incl. 5' & 3')....: 5
    Only terminal overlaps...........: No
    Circular products................: [59]
    Linear products..................: [74] [73] [73] [54] [54] [53] [15] [14] [14]
    >>> x.circular_products
    [Contig(o59)]
    >>> x.circular_products[0].seq.watson
    'CCCCCtgtgctgtgctctaTTTTTtattctggctgtatcGGGGGtacgatgctatactg'

    '''

    def __init__(self, dsrecs, limit = 25, only_terminal_overlaps=False, max_nodes=None):

        refresh = False
        cached  = None

        key = "|".join(sorted([str(r.seguid()) for r in dsrecs]))+str(limit)+str(only_terminal_overlaps)+str(max_nodes)

        if os.environ["pydna_cache"] in ("compare", "cached"):

            module_logger.info('open shelf file {}'.format(os.path.join(os.environ["pydna_data_dir"],"assembly")))

            cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "assembly"), protocol=cPickle.HIGHEST_PROTOCOL, writeback=False)

            module_logger.info('created key = {}'.format(key))
            module_logger.info( "pydna_cache = {}".format(os.environ["pydna_cache"]) )

            try:
                cached = cache[key]
            except:
                if os.environ["pydna_cache"] == "compare":
                    raise Exception("no result for this key!")
                else:
                    refresh = True

            cache.close()

        if refresh or os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):
            self.dsrecs    = dsrecs
            ''' Sequences fed to this class is stored in this property'''
            self.only_terminal_overlaps = only_terminal_overlaps
            ''' Consider only terminal overlaps?'''
            self.limit     = limit
            ''' The shortest common sub strings to be considered '''
            self.max_nodes = max_nodes or len(self.dsrecs)
            ''' The max number of nodes allowed. This can be reset to some other value'''
            self.key = key
            self.only_terminal_overlaps = only_terminal_overlaps
            self._assemble()

        if os.environ["pydna_cache"] == "compare":
            self._compare(cached)

        if refresh or os.environ["pydna_cache"] == "refresh":
            self._save()

        elif cached and os.environ["pydna_cache"] not in ("nocache", "refresh"):
            for key, value in cached.__dict__.items():
                setattr(self, key, value )
            cache.close()

    def _compare(self, cached):
        if str(self) != str(cached):
            module_logger.warning('Assembly error')

    def _save(self):
        cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "assembly"), protocol=cPickle.HIGHEST_PROTOCOL, writeback=False)
        cache[self.key] = self
        cache.close()

    def _assemble(self):

        for dr in self.dsrecs:
            if dr.name in ("",".", "<unknown name>", None):
                dr.name = "frag{}".format(len(dr))

        if self.only_terminal_overlaps:
            algorithm = terminal_overlap
        else:
            algorithm = common_sub_strings

        # analyze_overlaps
        cols = {}
        for dsrec in self.dsrecs:
            dsrec.features = [f for f in dsrec.features if f.type!="overlap"]
            dsrec.seq = Dseq(dsrec.seq.todata)
        rcs = {dsrec:dsrec.rc() for dsrec in self.dsrecs}
        matches=[]
        dsset=OrderedSet()

        for a, b in itertools.combinations(self.dsrecs, 2):
            match = algorithm( str(a.seq).upper(),
                               str(b.seq).upper(),
                               self.limit)
            if match:
                matches.append((a, b, match))
                dsset.add(a)
                dsset.add(b)
            match = algorithm( str(a.seq).upper(),
                               str(rcs[b].seq).upper(),
                               self.limit)
            if match:
                matches.append((a, rcs[b], match))
                dsset.add(a)
                dsset.add(rcs[b])
                matches.append((rcs[a], b, [(len(a)-sa-le,len(b)-sb-le,le) for sa,sb,le in match]))
                dsset.add(b)
                dsset.add(rcs[a])

        self.no_of_olaps=0

        for a, b, match in matches:
            for start_in_a, start_in_b, length in match:
                self.no_of_olaps+=1
                chksum = a[start_in_a:start_in_a+length].seguid()
                #assert chksum == b[start_in_b:start_in_b+length].seguid()

                try:
                    fcol, revcol = cols[chksum]
                except KeyError:
                    fcol = '#%02X%02X%02X' % (random.randint(175,255),random.randint(175,255),random.randint(175,255))
                    rcol = '#%02X%02X%02X' % (random.randint(175,255),random.randint(175,255),random.randint(175,255))
                    cols[chksum] = fcol,rcol

                qual      = {"note"             : ["olp_{}".format(chksum)],
                             "chksum"           : [chksum],
                             "ApEinfo_fwdcolor" : [fcol],
                             "ApEinfo_revcolor" : [rcol]}

                if not chksum in [f.qualifiers["chksum"][0] for f in a.features if f.type == "overlap"]:
                    a.features.append( SeqFeature( FeatureLocation(start_in_a,
                                                                   start_in_a + length),
                                                                   type = "overlap",
                                                                   qualifiers = qual))
                if not chksum in [f.qualifiers["chksum"][0] for f in b.features if f.type == "overlap"]:
                    b.features.append( SeqFeature( FeatureLocation(start_in_b,
                                                                   start_in_b + length),
                                                                   type = "overlap",
                                                                   qualifiers = qual))
        for ds in dsset:
            ds.features = sorted([f for f in ds.features], key = operator.attrgetter("location.start"))

        self.analyzed_dsrecs = list(dsset)


        # Create graph

        self.G=nx.MultiDiGraph(multiedges=True, name ="original graph" , selfloops=False)
        self.G.add_node( '5' )
        self.G.add_node( '3' )

        for i, dsrec in enumerate(self.analyzed_dsrecs):

            overlaps = sorted( {f.qualifiers['chksum'][0]:f for f in dsrec.features
                                if f.type=='overlap'}.values(),
                               key = operator.attrgetter('location.start'))

            if overlaps:
                overlaps = ([SeqFeature(FeatureLocation(0, 0),
                             type = 'overlap',
                             qualifiers = {'chksum':['5']})]+
                             overlaps+
                            [SeqFeature(FeatureLocation(len(dsrec),len(dsrec)),
                                        type = 'overlap',
                                        qualifiers = {'chksum':['3']})])

                for olp1, olp2 in itertools.combinations(overlaps, 2):

                    n1 = olp1.qualifiers['chksum'][0]
                    n2 = olp2.qualifiers['chksum'][0]

                    if n1 == '5' and n2=='3':
                        continue

                    s1,e1,s2,e2 = (olp1.location.start.position,
                                   olp1.location.end.position,
                                   olp2.location.start.position,
                                   olp2.location.end.position,)

                    source_fragment = Fragment(dsrec,s1,e1,s2,e2,i)

                    self.G.add_edge( n1, n2,
                                     frag=source_fragment,
                                     weight = s1-e1,
                                     i = i)

        #linear assembly

        linear_products=defaultdict(list)

        for path in all_simple_paths_edges(self.G, '5', '3', data=True, cutoff=self.max_nodes):

            pred_frag = copy(path[0][2].values().pop()['frag'])
            source_fragments = [pred_frag, ]

            if pred_frag.start2<pred_frag.end1:
                result=pred_frag[pred_frag.start2+(pred_frag.end1-pred_frag.start2):pred_frag.end2]
            else:
                result=pred_frag[pred_frag.end1:pred_frag.end2]

            for first_node, second_node, edgedict in path[1:]:

                edgedict = edgedict.values().pop()

                f  = copy(edgedict['frag'])

                f.alignment =  pred_frag.alignment + pred_frag.start2- f.start1
                source_fragments.append(f)

                if f.start2>f.end1:
                    result+=f[f.end1:f.end2]
                else:
                    result+=f[f.start2+(f.end1-f.start2):f.end2]

                pred_frag = f

            add=True
            for lp in linear_products[len(result)]:
                if (str(result.seq).lower() == str(lp.seq).lower()
                    or
                    str(result.seq).lower() == str(lp.seq.reverse_complement()).lower()):
                    add=False
            for dsrec in self.dsrecs:
                if (str(result.seq).lower() == str(dsrec.seq).lower()
                    or
                    str(result.seq).lower() == str(dsrec.seq.reverse_complement()).lower()):
                    add=False
            if add:
                linear_products[len(result)].append(Contig( result, source_fragments))

        self.linear_products = list(itertools.chain.from_iterable(linear_products[size] for size in sorted(linear_products, reverse=True)))


        # circular assembly

        self.cG = self.G.copy()
        self.cG.remove_nodes_from(('5','3'))
        #circular_products=defaultdict(list)
        circular_products={}

        for pth in all_circular_paths_edges(self.cG):

            ns = min(enumerate(pth), key = lambda x:x[1][2]['i'])[0]

            path = pth[ns:]+pth[:ns]

            pred_frag = copy(path[0][2]['frag'])

            source_fragments = [pred_frag, ]

            if pred_frag.start2<pred_frag.end1:
                result=pred_frag[pred_frag.start2+(pred_frag.end1-pred_frag.start2):pred_frag.end2]
            else:
                result=pred_frag[pred_frag.end1:pred_frag.end2]

            result.seq = Dseq(str(result.seq))

            for first_node, second_node, edgedict in path[1:]:

                f  = copy(edgedict['frag'])

                f.alignment =  pred_frag.alignment + pred_frag.start2- f.start1
                source_fragments.append(f)

                if f.start2>f.end1:
                    nxt = f[f.end1:f.end2]
                else:
                    nxt =f[f.start2+(f.end1-f.start2):f.end2]
                nxt.seq = Dseq(str(nxt.seq))
                result+=nxt

                pred_frag = f

            #add=True
            #for cp in circular_products[len(result)]:
            #    if (str(result.seq).lower() in str(cp.seq).lower()*2
            #        or
            #        str(result.seq).lower() == str(cp.seq.reverse_complement()).lower()*2):
            #        pass
            #        add=False
            #        print "##--"
            #if add:
            #    circular_products[len(result)].append( Contig( Dseqrecord(result, circular=True), source_fragments))

            r = Dseqrecord(result, circular=True)
            circular_products[r.cseguid()] = Contig(r, source_fragments )


        #self.circular_products = list(itertools.chain.from_iterable(circular_products[size] for size in sorted(circular_products, reverse=True)))
        self.circular_products = sorted(circular_products.values(), key=len, reverse=True)


    def __repr__(self):
        return   ( "Assembly:\n"
                   "Sequences........................: {sequences}\n"
                   "Sequences with shared homologies.: {analyzed_dsrecs}\n"
                   "Homology limit (bp)..............: {limit}\n"
                   "Number of overlaps...............: {no_of_olaps}\n"
                   "Nodes in graph(incl. 5' & 3')....: {nodes}\n"
                   "Only terminal overlaps...........: {pr}\n"
                   "Circular products................: {cp}\n"
                   "Linear products..................: {lp}"    ).format(sequences       = " ".join("[{}]".format(len(x)) for x in self.dsrecs),
                                                                         analyzed_dsrecs = " ".join("[{}]".format(len(x)) for x in self.analyzed_dsrecs),
                                                                         limit           = self.limit,
                                                                         no_of_olaps     = self.no_of_olaps,
                                                                         nodes           = self.G.order(),
                                                                         pr              = {True:"Yes",False:"No"}[self.only_terminal_overlaps],
                                                                         cp              = " ".join("[{}]".format(len(x)) for x in self.circular_products),
                                                                         lp              = " ".join("[{}]".format(len(x)) for x in self.linear_products))

if __name__=="__main__":
    import doctest
    doctest.testmod()
