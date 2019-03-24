import textwrap as _textwrap
import networkx as _nx
from pydna._pretty import pretty_str    as _pretty_str
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.utils import rc

class Contig(_Dseqrecord):
    '''This class holds information about a DNA assembly. This class is instantiated by
    the :class:`Assembly` class and is not meant to be used directly.

    '''

    def __init__(self,
                 record,
                 *args,
                 graph = None,
                 path  = None,
                 nodemap = None,
                 **kwargs):

        super().__init__(record, *args, **kwargs)
        self.graph = graph
        self.path = path or []
        self.nodemap = nodemap or {}

    def __repr__(self):
        return "Contig({}{})".format({True:"-", False:"o"}[self.linear],len(self))

        
    def _repr_pretty_(self, p, cycle):
        '''returns a short string representation of the object'''
        p.text("Contig({}{})".format({True:"-", False:"o"}[self.linear],len(self)))

           
    def _repr_html_(self):
        return "<pre>"+self.small_fig()+"</pre>"
    
    
    def reverse_complement(self):
        answer = type(self)(super().reverse_complement())
        answer.graph = _nx.relabel_nodes(self.graph.reverse(), self.nodemap)
        for edge in answer.graph.edges():
            answer.graph.edges[edge]["seq"] = rc(answer.graph.edges[edge]["seq"])
            n = answer.graph.edges[edge]["name"]
            n = n[:-3] if n.endswith("_rc") else "{}_rc".format(n)[:13]
            answer.graph.edges[edge]["name"] = n
            #n = "{}_rc".format(n) if not n.endswith("_rc") else answer.graph.edges[edge]["name"][:-3]
            l=len(answer.graph.edges[edge]["seq"])
            answer.graph.edges[edge]["start"] = l-answer.graph.edges[edge]["end"]   - answer.graph.node[edge[0]]["length"]  
            answer.graph.edges[edge]["end"]   = l-answer.graph.edges[edge]["start"] - answer.graph.node[edge[1]]["length"]
            answer.graph.edges[edge]["length"] = answer.graph.edges[edge]["end"] - answer.graph.edges[edge]["start"]
        answer.path  = [self.nodemap[n] for n in self.path[::-1]]
        nm = self.nodemap
        answer.nodemap = {nm[i]:i for i in nm}
        return answer
    
    
    rc = reverse_complement


    def detailed_fig(self):
        """Returns a text representation of the assembled fragments.

        Linear:

        ::

            acgatgctatactgCCCCCtgtgctgtgctcta
                               TGTGCTGTGCTCTA
                               tgtgctgtgctctaTTTTTtattctggctgtatc



        Circular:

        ::

            ||||||||||||||
            acgatgctatactgCCCCCtgtgctgtgctcta
                               TGTGCTGTGCTCTA
                               tgtgctgtgctctaTTTTTtattctggctgtatc
                                                  TATTCTGGCTGTATC
                                                  tattctggctgtatcGGGGGtacgatgctatactg
                                                                       ACGATGCTATACTG


        """


        fig=""
        fragmentposition=0
        nodeposition = self.graph[self.path[0]][self.path[1]]["length"]
        nodeposition = 0
        mylist = []
        for u, v, e in [(u, v, self.graph[u][v]) for u, v in zip(self.path, self.path[1:])]:
            nodeposition += e["length"]
            fragmentposition-=e["start"]
            mylist.append([fragmentposition, str(e["seq"])])
            mylist.append([nodeposition, v.upper()])
            fragmentposition+=e["end"]
                
        if self.circular:
            nodeposition = self.graph[self.path[0]][self.path[1]]["start"]
            nodelength   = len(v) #len(self.path[-1])
            mylist= [[nodeposition,"|"*nodelength]] + mylist
        else:
            mylist = mylist[:-1]
        
        firstpos = -1 * min(0, min(mylist)[0])
        
        for p,s in mylist:
            fig+="{}{}\n".format(" "*(p+firstpos), s)
            
        return _pretty_str(fig)
    

    detailed_figure = detailed_fig
    

    def small_fig(self):
        '''Returns a small ascii representation of the assembled fragments. Each fragment is
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

        edges = [self.graph[u][v] for u,v in zip(self.path, self.path[1:])]
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

            f = edges[0]
            space2 = len(f["name"])


            fig = ("{name}|{o2:>2}\n"
                   "{space2} \\/\n"
                   "{space2} /\\\n").format(name   = f["name"],
                                            o2     = self.graph.node[self.path[1]]["length"],
                                            space2 = " "*space2)
            space = space2 #len(f.name)

            for i,f in enumerate( edges[1:-1] ):
                name = "{o1:>2}|{name}|".format(o1   = self.graph.node[self.path[i+1]]["length"],
                                                name = f["name"])
                space2 = len(name)
                
                fig +=("{space} {name}{o2:>2}\n"
                       "{space} {space2}\\/\n"
                       "{space} {space2}/\\\n").format( name   = name,
                                                        o2     = self.graph.node[self.path[i+2]]["length"],
                                                        space  = " "*space,
                                                        space2 = " "*space2)
                space +=space2
            f = edges[-1]
            fig += ("{space} {o1:>2}|{name}").format(name  = f["name"],
                                                     o1    = self.graph.node[self.path[-2]]["length"],
                                                     space = " "*(space))



        else:
            r'''
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

            f = edges[0]
            
            space = len(f["name"])+3
            
            fig =(" -|{name}|{o2:>2}\n"
                  "|{space}\\/\n"
                  "|{space}/\\\n").format( name = f["name"],
                                           o2 = self.graph.node[self.path[1]]["length"],
                                           space = " "*space )
            
            for i, f in enumerate( edges[1:] ):
                name= "{o1:>2}|{name}|".format(o1   = self.graph.node[self.path[i+1]]["length"],
                                               name = f["name"])
                space2 = len(name)
                fig +=("|{space}{name}{o2:>2}\n"
                       "|{space}{space2}\\/\n"
                       "|{space}{space2}/\\\n").format( o2     = self.graph.node[self.path[i+2]]["length"],
                                                        name   = name,
                                                        space  = " "*space,
                                                        space2 = " "*space2 )
                space +=space2

            fig +="|{space}{o1:>2}-\n".format(space=" "*(space), o1=self.graph.node[self.path[0]]["length"])
            fig +="|{space}   |\n".format(space=" "*(space))
            fig +=" {space}".format(space="-"*(space+3))
        return _pretty_str(_textwrap.dedent(fig))


    figure = small_figure = small_fig
    
    

if __name__=="__main__":
    import os as _os
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"]=""
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"]=cached
