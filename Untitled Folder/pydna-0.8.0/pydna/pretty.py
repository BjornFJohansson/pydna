#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''The pretty_str class is same as unicode but has a _repr_pretty_ method for
   for nicer string output in the IPython shell'''

class pretty_str(unicode):
    ''' Thanks to Min RK, UC Berkeley for this'''
    def _repr_pretty_(self, p, cycle):
        p.text(self)

if __name__=="__main__":
    a = ''' a
            b'''
    print a
    b = pretty_str(a)
    print b

    c = ''' ä
            ö'''

    print c
    d = u'''
                                 Taq (rate {rate} nt/s) 35 cycles             |{size}bp
                                 95.0°C    |95.0°C                 |      |SantaLucia 1998
                                 |_________|_____          72.0°C  |72.0°C|SaltC {saltc:2}mM
                                 | 03min00s|30s  \         ________|______|
                                 |         |      \ {ta}°C/{0:2}min{1:2}s| 5min |
                                 |         |       \_____/         |      |
                                 |         |         30s           |      |4-12°C'''

    print d