#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''The pretty_str class is same as str but has a _repr_pretty_ method for
   for nicer string output in the IPython shell'''

class pretty_str(str):
    ''' Thanks to Min RK, UC Berkeley for this'''
    def _repr_pretty_(self, p, cycle):
        p.text(self)
    def _repr_html_(self):
        return "123"

class pretty_unicode(str):
    def _repr_pretty_(self, p, cycle):
        p.text(self)

class pretty_string(str):
    def _repr_pretty_(self, p, cycle):
        p.text(self)
        
        
from IPython.lib.latextools import latex_to_png

class Circle(object):

    def __init__(self, radius):
        self.radius = radius

    def _repr_pretty_(self, p, cycle):
        p.text(u"\u25CB")

    def _repr_html_(self):
        return "<h1>Cirle: radius=%s</h1>" % self.radius

    def _repr_svg_(self):
        return """<svg>
<circle cx="100" cy="50" r="40" stroke="black" stroke-width="2" fill="red"/>
</svg>"""

    def _repr_png_(self):
        return latex_to_png('$\circle$')

    def _repr_javascript_():
        
        
        
if __name__=="__main__":
    import doctest
    doctest.testmod()