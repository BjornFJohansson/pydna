import pytest

def test_contig(monkeypatch):
    monkeypatch.setenv("pydna_cached_funcs", "")

    from pydna import contig    
    from pydna.assembly import Assembly
    from pydna.dseqrecord import Dseqrecord
    
    a = Dseqrecord("acgatgctatactgCCCCCtgtgctgtgctcta",   name="one")
    b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatc",  name="two")
    c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactg", name="three")
    asm = Assembly((a,b,c), limit=14)

    cnt = asm.circular_products[0]
    
    assert repr(cnt) == "Contig(o59)"
    
    assert cnt.detailed_figure() == str(   "||||||||||||||\n"
                                           "acgatgctatactgCCCCCtgtgctgtgctcta\n"
                                           "                   TGTGCTGTGCTCTA\n"
                                           "                   tgtgctgtgctctaTTTTTtattctggctgtatc\n"
                                           "                                      TATTCTGGCTGTATC\n"
                                           "                                      tattctggctgtatcGGGGGtacgatgctatactg\n"
                                           "                                                           ACGATGCTATACTG\n")



    from textwrap import indent

    fig = """ -|one|14
|      \\/
|      /\\
|      14|two|15
|             \\/
|             /\\
|             15|three|14
|                      \\/
|                      /\\
|                      14-
|                         |
 -------------------------"""
    assert fig == cnt.small_fig()
    
    cnt2 = asm.linear_products[0]

    
    fig = ('one|14\n'
           '    \\/\n'
           '    /\\\n'
           '    14|two|15\n'
           '           \\/\n'
           '           /\\\n'
           '           15|three')

    assert fig == cnt2.small_fig() == cnt2.figure() == cnt2.small_figure()
    
    assert repr(cnt2) == 'Contig(-73)'

    #print(repr(cnt2._repr_html_()))
    
    assert cnt2._repr_html_() == '<pre>one|14\n    \\/\n    /\\\n    14|two|15\n           \\/\n           /\\\n           15|three</pre>'

    from unittest.mock import MagicMock
    
    pp = MagicMock()
    
    cnt2._repr_pretty_(pp, None)
    
    pp.text.assert_called_with('Contig(-73)')
    
def test_reverse_complement(monkeypatch):
    from pydna._pretty import pretty_str
    from pydna.assembly import Assembly
    from pydna.dseqrecord import Dseqrecord
    a = Dseqrecord("acgatgctatactgtgCCNCCtgtgctgtgctcta")
    b = Dseqrecord("tgtgctgtgctctaTTTTTTTtattctggctgtatc")
    c = Dseqrecord("tattctggctgtatcGGGGGtacgatgctatactgtg")
    a.name="aaa"
    b.name="bbb"
    c.name="ccc"
    x = Assembly((a,b,c), limit=14)
    y = x.circular[0]
    
    yfig = '''
 -|aaa|14
|      \\/
|      /\\
|      14|bbb|15
|             \\/
|             /\\
|             15|ccc|16
|                    \\/
|                    /\\
|                    16-
|                       |
 -----------------------
     '''[1:].rstrip()
     
     
     
     
    ydfig= pretty_str('''
||||||||||||||||
acgatgctatactgtgCCNCCtgtgctgtgctcta
                     TGTGCTGTGCTCTA
                     tgtgctgtgctctaTTTTTTTtattctggctgtatc
                                          TATTCTGGCTGTATC
                                          tattctggctgtatcGGGGGtacgatgctatactgtg
                                                               ACGATGCTATACTGTG
    '''[1:].rstrip()+"\n")
    
    

    assert y.figure() == yfig
    assert y.detailed_figure() == ydfig
    
    z=y.rc()
    
    zfig = '''
 -|ccc_rc|15
|         \\/
|         /\\
|         15|bbb_rc|14
|                   \\/
|                   /\\
|                   14|aaa_rc|16
|                             \\/
|                             /\\
|                             16-
|                                |
 --------------------------------
     '''[1:].rstrip()
     
     
    zdfig= '''
||||||||||||||||
cacagtatagcatcgtaCCCCCgatacagccagaata
                      GATACAGCCAGAATA
                      gatacagccagaataAAAAAAAtagagcacagcaca
                                            TAGAGCACAGCACA
                                            tagagcacagcacaGGNGGcacagtatagcatcgt
                                                               CACAGTATAGCATCGT
    '''[1:].rstrip()+"\n"
    
    assert z.figure() == zfig
    assert z.detailed_figure() == zdfig
    

if __name__ == '__main__':
    pytest.main([__file__, "-vv", "-s","--cov=pydna","--cov-report=html"])

