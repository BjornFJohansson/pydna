#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_myprimers(monkeypatch):
    monkeypatch.setenv("pydna_primers", "primers_linux_line_endings.txt")
    from pydna import myprimers
    from pydna.parsers import parse_primers
    from importlib import reload

    reload(myprimers)
    newlist = parse_primers("primers_linux_line_endings.txt")[::-1]
    primerdict = myprimers.primerdict()

    assert len(primerdict) == 4
    primer_list = myprimers.primerlist()
    assert primer_list == newlist


def test_prepend_primerlist(monkeypatch):
    monkeypatch.setenv("pydna_primers", "primers_linux_line_endings.txt")
    from pydna import myprimers
    from pydna.parsers import parse_primers
    from importlib import reload

    # >3_primer
    # aaaaaaaa
    # >2_primer
    # cccccccc
    # >1_primer
    # gggggggg
    # >0_primer
    # tttttttt

    oldlist = myprimers.primerlist()

    reload(myprimers)

    newlist = parse_primers("""
                            >abc
                            aaa
                            >efg
                            ttt
                            """)

    np = myprimers.prepend_primerlist(newlist, oldlist)

    assert [s.name for s in parse_primers(np)] == ["5_abc",
                                                   "4_efg"]


def test_check_primer_list(monkeypatch):
    monkeypatch.setenv("pydna_primers",
                       "primers_linux_line_endings_not_unique.txt")
    from pydna import myprimers
    from pydna.parsers import parse_primers
    from importlib import reload

    reload(myprimers)

    from textwrap import dedent

    m = dedent("""\
    6 primers, 5 unique primer sequences
    1 primer(s) without sequence (N)
    Wrong number: 5 6_primer	atatatat
    1_primer 3_primer gggggggg""")

    pl = myprimers.primerlist()

    print(myprimers.check_primer_list(pl) == m)





if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
