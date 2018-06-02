#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 31 10:58:29 2017

@author: bjorn
"""
import io
from unittest.mock import patch
from urllib.request import urlopen
from Bio import Entrez

@patch('Bio.Entrez._urlopen')
def openPatch(urlopenMock):
    urlopenMock.return_value = open("X60065.gb", "rb")
    Entrez._urlopen = urlopenMock
    Entrez.email = "bjornjobb@gmail.com"
    Entrez.tool  = "pydna"
    #https://www.ncbi.nlm.nih.gov/nuccore/5
    handle = Entrez.efetch(db="nuccore",
                           id="AJ515744.1",
                           rettype="gb",
                           retmode="text")
    from Bio import SeqIO
    result = SeqIO.read(handle, "genbank")
    print(result)

openPatch()