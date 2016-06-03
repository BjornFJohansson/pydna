#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''The pretty_str class is same as str but has a _repr_pretty_ method for
   for nicer string output in the IPython shell'''

class pretty_str(str):
    ''' Thanks to Min RK, UC Berkeley for this'''
    def _repr_pretty_(self, p, cycle):
        p.text(self)

class pretty_unicode(unicode):
    def _repr_pretty_(self, p, cycle):
        p.text(self)

class pretty_string(str):
    def _repr_pretty_(self, p, cycle):
        p.text(self)

if __name__=="__main__":
    import pydna

    print pydna.read("/home/bjorn/Desktop/python_packages/pydna/pydna/pydna_read_test.txt").format()
    print pydna.read("/home/bjorn/Desktop/python_packages/pydna/pydna/pydna_read_test2.txt").format()[3270:3281]
    import sys;sys.exit(42)
    import StringIO
    from Bio import SeqIO
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

    import textwrap, re

    raw = open("pydna_read_test.txt", 'rU').read()

    pattern =  r"(?:>.+\n^(?:^[^>]+?)(?=\n\n|>|LOCUS|ID))|(?:(?:LOCUS|ID)(?:(?:.|\n)+?)^//)"

    rawseq = re.findall(pattern, textwrap.dedent(raw + "\n\n"), flags=re.MULTILINE).pop(0)

    handle = StringIO.StringIO(raw)

    sr = SeqIO.read(handle, "genbank", alphabet=IUPACAmbiguousDNA())

    s = sr.format("gb").strip()

    print pretty_string(s[:55]+"circular"+s[63:])[3200:3300]

    from pydna import Dseqrecord

    v = Dseqrecord(sr)

    print v.format("gb")

