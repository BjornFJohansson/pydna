#!/usr/bin/env python
# -*- coding: utf-8 -*-

def str SmallestRotation(str s):
    cdef str prev="",w="",ds=""
    cdef int rep=0,old=0,k=0,i=0,j=0
    ds=s+s
    lends=len(ds)
    while k < lends:
        i,j = k,k+1
        while j < lends and ds[i] <= ds[j]:
            i = (ds[i] == ds[j]) and i+1 or k
            j += 1
        while k < i+1:
            k += j-i
            prev=w
            w=ds[old:k]
            old = k
            if w == prev:
                rep += 1
            else:
                prev,rep = w,1
            if len(w)*rep == len(s):
                return w*rep