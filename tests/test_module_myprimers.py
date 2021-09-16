#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

monkeypatch = pytest.MonkeyPatch()


def test_PrimerList_init(monkeypatch):

    monkeypatch.setenv("pydna_primers", "primers_linux_line_endings.txt")

    from pydna.parsers import parse_primers

    primer_source = parse_primers("primers_linux_line_endings.txt")[::-1]

    from pydna import myprimers

    from importlib import reload

    reload(myprimers)

    pl1 = myprimers.PrimerList()

    assert len(pl1) == 4

    assert pl1 == primer_source

    pl2 = myprimers.PrimerList(path="primers_linux_line_endings.txt")

    assert pl1 == pl2

    newlist = parse_primers("""
                            >abc
                            aaa
                            >efg
                            ttt
                            """)

    np = pl1.assign_numbers_to_new_primers(newlist)

    assert [s.name for s in parse_primers(np)] == ["5_abc",
                                                   "4_efg"]

    newlist = parse_primers("""
                            >abc
                            aaa
                            >efg
                            tttttttt
                            """)

    np = pl1.assign_numbers_to_new_primers(newlist)

    assert [s.name for s in parse_primers(np)] == ["4_abc"]

    with pytest.raises(ValueError):
        pl1.pydna_code_from_list(newlist)

    import textwrap

    code = textwrap.dedent("""\
    p = {}

    p[0], p[1], p[2], p[3] = parse_primers('''

    >0_primer
    tttttttt

    >1_primer
    gggggggg

    >2_primer
    cccccccc

    >3_primer
    aaaaaaaa

    ''')""")

    assert pl1.pydna_code_from_list(pl1) == code

    assert pl1.pydna_code_from_indices([0, 1, 2, 3]) == code

    assert pl1.pydna_code_from_accessed() == code

    pl3 = myprimers.PrimerList(path="primers_linux_line_endings.txt")

    assert pl3.accessed == []
    pl3[1] = primer_source[1]
    assert pl3.accessed == [1]

    pl3.accessed == []
    assert pl3[2:4] == primer_source[2:4]
    assert pl3.accessed == [1, 2, 3]

    with pytest.raises(ValueError):
        myprimers.PrimerList(identifier="/")

    with pytest.raises(ValueError):
        pl3[2] = primer_source[1]

    with pytest.raises(IndexError):
        pl3[999] = primer_source[1]

    from unittest import mock
    
    subp = mock.MagicMock()
    
    monkeypatch.setattr("sys.platform", "linux")
    monkeypatch.setattr("subprocess.run", subp)
    
    pl3.open_folder()

    subp.assert_called_with(["xdg-open", pl3.path.parent])


def test_check_primer_numbers(monkeypatch):
    from pydna import myprimers
    pl = myprimers.PrimerList(path="primers_linux_line_endings_not_unique.txt")
    assert myprimers.check_primer_numbers(pl) == [pl[5]]


def test_undefined_sequence(monkeypatch):
    from pydna import myprimers
    pl = myprimers.PrimerList(path="primers_linux_line_endings_not_unique.txt")
    assert myprimers.undefined_sequence(pl) == [pl[2]]


def test_find_duplicate_primers(monkeypatch):
    from pydna import myprimers
    pl = myprimers.PrimerList(path="primers_linux_line_endings_not_unique.txt")
    assert myprimers.find_duplicate_primers(pl) == [[pl[1], pl[3]]]


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
