#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nose
from pydna import parse, read

from Bio.SeqIO import read as BPread
from Bio.SeqIO import parse as BPparse

def test_pth1():

    q = BPread("read1.gb", "gb")
    w = BPread("read2.gb", "gb")
    e = BPread("read3.fasta", "fasta")
    r = BPread("read4.fasta", "fasta")

    q.format("gb")
    w.format("gb")

    a, b = BPparse("pth1.txt", "gb")

    x, y = parse("pth1.txt")

    assert "".join(a.format("gb").splitlines()[1:]) == "".join(x.format("gb").splitlines()[1:])
    assert "".join(b.format("gb").strip().splitlines()[4:]) == "".join(y.format("gb").splitlines()[4:])


if __name__ == '__main__':
    nose.runmodule()

