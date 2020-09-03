#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 08:03:27 2017

@author: bjorn
"""

from pydna.readers import read

template = read(">t\\ntacactcaccgtctatcattatctactatcgactgtatcatctgatagcac")
p1 = read(">p1\\ntacactcaccgtctatcattatc", ds=False)
p2 = read(">p2\\ngtgctatcagatgatacagtcg", ds=False)
# ann = pydna.Anneal((p1, p2), template)
