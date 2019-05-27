#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

def test_genbankfile():
    
    from pydna import genbankrecord
    
    gbr = genbankrecord.GenbankRecord("aaa")
    
    assert gbr.hyperlink == "<a href='https://www.ncbi.nlm.nih.gov/nuccore/accession?from=&to=&strand=1' target='_blank'>accession</a>"
    
    assert repr(gbr) == 'Gbank(accession)(-3)'
    
    assert gbr._repr_html_()=="<a href='https://www.ncbi.nlm.nih.gov/nuccore/accession?from=&to=&strand=1' target='_blank'>accession</a>"
    
    from unittest.mock import MagicMock
    
    pp = MagicMock()
    
    gbr._repr_pretty_(pp, None)
    
    pp.text.assert_called_with('Gbank(accession)(-3)')
    
    gbr = genbankrecord.GenbankRecord("aaa", start=1, stop=2)

    assert gbr.hyperlink == "<a href='https://www.ncbi.nlm.nih.gov/nuccore/accession?from=1&to=2&strand=1' target='_blank'>accession 1-2</a>"
    
    gbr_rc = gbr.rc()
    
    assert gbr_rc.strand==2
    
    from Bio.Seq import Seq
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
    from pydna.seqrecord import SeqRecord
    
    arg = SeqRecord(Seq("aaa", IUPACAmbiguousDNA()))
    
    genbankrecord.GenbankRecord.from_SeqRecord(arg)

if __name__ == '__main__':
    pytest.main([__file__, "-vv", "-s","--cov=pydna","--cov-report=html"])



