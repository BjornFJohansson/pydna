#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''The pretty_str class is same as str but has a _repr_pretty_ method for
   for nicer string output in the IPython shell'''

class pretty_str(str):
    ''' Thanks to Min RK, UC Berkeley for this'''
    def _repr_pretty_(self, p, cycle):
        p.text(self)

class pretty_unicode(str):
    def _repr_pretty_(self, p, cycle):
        p.text(self)

class pretty_string(str):
    def _repr_pretty_(self, p, cycle):
        p.text(self)
        
if __name__=="__main__":
    import doctest
    doctest.testmod()