
import pytest


def test_amplicon():
    
    from pydna.amplify    import Anneal
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer     import Primer
    
    template = Dseqrecord("AAAtacactcaccgtctatcattatctactatcgactgtatcatctgatagcacTTT")
    
    p1 = Primer("CCCtacactcaccgtctatcattatc")
    p2 = Primer("GGGgtgctatcagatgatacagtcg")
    
    ann = Anneal((p1,p2),template)
    
    prod = ann.products[0]
    
    assert repr(prod) == 'Amplicon(57)'
    
    assert prod._repr_html_() == 'Amplicon(57)'
    
    from unittest.mock import MagicMock
    
    pp = MagicMock()
    
    prod._repr_pretty_(pp, None)
    
    #assert pp.text.assert_called_with('Amplicon(57)')
    
    
    fig='''    5tacactcaccgtctatcattatc...cgactgtatcatctgatagcac3
                               |||||||||||||||||||||| tm 55.9 (dbd) 60.5
                              3gctgacatagtagactatcgtgGGG5
 5CCCtacactcaccgtctatcattatc3
     ||||||||||||||||||||||| tm 54.6 (dbd) 58.8
    3atgtgagtggcagatagtaatag...gctgacatagtagactatcgtg5'''
    
    import textwrap
    
    assert prod.figure() == textwrap.dedent(fig)
    
    assert prod.program() == prod.taq_program()
    
    assert prod.pfu_sso7d_program() == prod.dbd_program()

    from pydna.amplicon import Amplicon

    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
    from pydna.seqrecord import SeqRecord
    
    arg = SeqRecord(Seq("aaa", IUPACAmbiguousDNA()))
    
    x = Amplicon.from_SeqRecord(arg)
    



def test_amplicon_dbd():
    
    from pydna.amplify    import Anneal
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer     import Primer
    
    template = Dseqrecord("GCGTCCAGCGGCTGCCCGAGGCGCCAAGTGCCCGGGCCGAGCCCGCATCTGAGGCCGCCGCGGGC")
    
    p1 = Primer("GCGTCCAGCGGCTGCCCGAGG")
    p2 = Primer("GCCCGCGGCGGCCTCAGATGCGG")
    
    ann = Anneal((p1,p2),template)
    
    prod = ann.products[0]
    
    assert repr(prod) == 'Amplicon(65)'

    fig =( r'\nPfu-Sso7d (rate 15s/kb)\n'
            'Two-step|    30 cycles |      |65bp\n'
            '98.0°C  |98.0C         |      |Tm formula: Pydna tmbresluc\n'
            '_____ __|_____         |      |SaltC 50mM\n'
            '00min30s|10s  \        |      |Primer1C 1.0µM\n'
            '        |      \ 72.0°C|72.0°C|Primer2C 1.0µM\n'
            '        |       \______|______|GC 81%\n'
            '        |       0min 0s|10min |4-12°C\n')
     
    assert str(prod.pfu_sso7d_program()) == fig
    
def test_amplicon_dbd_low_gc():
    
    from pydna.amplify    import Anneal
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer     import Primer
    
    template = Dseqrecord("AAAATATTTTTATACATAATACAATTGTATATTCTTAAATAAAAAATACGTCATC")
    
    p1 = Primer("AAAATATTTTTATACAT")
    p2 = Primer("GATGACGTATTTTTTAT")
    
    ann = Anneal((p1,p2),template)
    
    prod = ann.products[0]
    
    assert repr(prod) == 'Amplicon(55)'

    fig =(r'\nPfu-Sso7d (rate 15s/kb)\n'
             'Two-step|    30 cycles |      |65bp\n'
             '98.0°C  |98.0C         |      |Tm formula: Pydna tmbresluc\n'
             '_____ __|_____         |      |SaltC 50mM\n'
             '00min30s|10s  \        |      |Primer1C 1.0µM\n'
             '        |      \ 72.0°C|72.0°C|Primer2C 1.0µM\n'
             '        |       \______|______|GC 81%\n'
             '        |       0min 0s|10min |4-12°C\n')
     
    assert str(prod.pfu_sso7d_program())

if __name__ == '__main__':
    pytest.main([__file__, "-v", "-s","--cov=pydna","--cov-report=html"])
