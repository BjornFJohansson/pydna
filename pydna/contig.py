import textwrap as _textwrap
from pydna._pretty import pretty_str    as _pretty_str
from pydna.dseqrecord import Dseqrecord as _Dseqrecord

class Contig(_Dseqrecord):
    '''This class holds information about a DNA assembly. This class is instantiated by
    the :class:`Assembly` class and is not meant to be instantiated directly.

    '''

    def __init__(self,
                 record,
                 *args,
                 source_fragments=[],
                 **kwargs):

        super().__init__(record, *args, **kwargs)
        self.source_fragments = source_fragments
        self.number_of_fragments = len(self.source_fragments)

    def __repr__(self):
        return "Contig({}{})".format({True:"-", False:"o"}[self.linear],len(self))
        
    def _repr_pretty_(self, p, cycle):
        '''returns a short string representation of the object'''
        p.text("Contig({}{})".format({True:"-", False:"o"}[self.linear],len(self)))
            
    def _repr_html_(self):
        return "<pre>"+self.small_fig()+"</pre>"
        #"Contig({}{})".format({True:"-", False:"o"}[self.linear],len(self))

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
        return _pretty_str(_textwrap.dedent(fig))


if __name__=="__main__":
    import os        as _os
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True)
    _os.environ["pydna_cache"]=cache
