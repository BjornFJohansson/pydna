#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# import pydna
# template = pydna.read(">t\ntacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
# p1 = pydna.read(">p1\ntacactcaccgtctatcattatc", ds = False)
# p2 = pydna.read(">p2\ngtgctatcagatgatacagtcg", ds = False)
# ann = pydna.Anneal((p1, p2), template)

import pydna

a = pydna.Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
b = pydna.Dseqrecord(
    "ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc"
)
c = pydna.Dseqrecord(
    "tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt"
)
primer_pairs = pydna.assembly_primers([a, b, c], circular=True)
p = []
for t, (f, r) in zip([a, b, c], primer_pairs):
    p.append(pydna.pcr(f, r, t))
obj = pydna.Assembly(p)
# print(assemblyobj)

# @memorize("function")
# def times_two(n):
#    return n*2
#
# @memorize("method")
# class Doubler(object):
#    def __init__(self, n):
#        self.n=n
#        self.key=str(n)
#    def result(self):
#        return self.n*2
#
# def peek(fn):
#    import shelve
#    ca = shelve.open("/home/bjorn/.local/share/pydna/{}".format(fn))
#    for key in ca.keys():
#        print("key         :", key)
#        print("cached value:", ca[key])
#        print()
#    ca.close()

# if __name__=="__main__":
#    print(_os.environ["pydna_cache"])
#
#    print("times_two(2) = ",times_two(2))
#    peek("function")
#
#    x=Doubler(3)
#    print( "x.result() = ", x.result() )
#    peek("method")

#    import pydna
#    a=pydna.Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
#    b=pydna.Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc")
#    c=pydna.Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt")
#    primer_pairs = pydna.assembly_primers([a,b,c], circular = True)
#    p=[]
#    for t, (f,r) in zip([a,b,c], primer_pairs):
#        p.append(pydna.pcr(f,r,t))
#    obj = pydna.Assembly(p)
# print(assemblyobj)


# from mementos import MementoMetaclass
#
## use the whatever memoize implementation you like
#
# class Memoize(type):
#    @memorize("hej")
#    def __call__(cls, *args, **kwargs):
#        return super().__call__(*args)
#
# class Mather(object, metaclass = Memoize):
#    def __init__(self, *args, **kwargs):
#        self.values = args
#    def sum(self):
#        return sum(self.values)
#
#
# o1 = Mather(1,2,3, type="float")


# class MyMetaClass(type):
#    @classmethod
#    def __prepare__(metacls, name, bases, fn="", **kargs):
#        print(name, bases, kargs)
#        return super().__prepare__(name, bases, **kargs)
#    def __new__(metacls, name, bases, namespace, fn="", **kargs):
#        return super().__new__(metacls, name, bases, namespace)
#    def __init__(cls, name, bases, namespace, fn="", **kargs):
#        cls.fn=fn
#        super().__init__(name, bases, namespace)
#    @memorize("hejhejhej")
#    def __call__(cls, *args, **kwargs):
#        print(cls)
#        return super().__call__(*args)
#
# class Mather2(object, metaclass=MyMetaClass, fn="hejhej"):
#    def __init__(self, *args, **kwargs):
#        self.values = args
#    def sum(self):
#        return sum(self.values)
#
# x=Mather2(1,2,3, type="float")
