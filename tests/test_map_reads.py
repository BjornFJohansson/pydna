#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import sys

from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydna.readers import read

from Bio.SeqIO import read as abiread


def test_map():

    traces = []

    import glob
    for name in glob.glob('*.ab1'):
        traces.append( abiread( name, "abi") )

    for t in traces:

        d = Dseqrecord(t.seq)

        if "ITVFFKEYPYDVPDYAIEGIFHAT" in d:

            tag = "tat cca tat gac gtt cca gac tat gca"
            trc = "ata ggt ata ctg caa ggt ctg ata cgt"[::-1]

            s = Dseqrecord(Dseq(tag,trc))
            sl = s.find_aa("YPYDVPDYA")
            assert str( s[sl].seq.translate() ) == "YPYDVPDYA"
            assert "YPYDVPDYA" in s

            tag = "AAA tat cca tat gac gtt cca gac tat gca"
            trc = "    ata ggt ata ctg caa ggt ctg ata cgt"[::-1]

            s = Dseqrecord(Dseq(tag,trc))
            sl = s.find_aa("YPYDVPDYA")
            assert str( s[sl].seq.translate() ) == "YPYDVPDYA" 
            assert "YPYDVPDYA" in s 

            tag = "    tat cca tat gac gtt cca gac tat gca"
            trc = "AAA ata ggt ata ctg caa ggt ctg ata cgt"[::-1]

            s = Dseqrecord(Dseq(tag,trc))
            sl = s.find_aa("YPYDVPDYA")
            assert str( s[sl].seq.translate() ) == "YPYDVPDYA" 
            assert "YPYDVPDYA" in s 

            tag = "    tat cca tat gac gtt cca gac tat gca"
            trc = "AAA ata ggt ata ctg caa ggt ctg ata cgt"[::-1]

            s = Dseqrecord(Dseq(tag,trc))
            sl = s.find_aa("YPYDVPDYA")
            assert str( s[sl].seq.translate() ) == "YPYDVPDYA" 

            tag = "tat cca tat gac gtt cca gac tat gca"
            trc = "ata ggt ata ctg caa ggt ctg ata cgt"[::-1]

            tag, trc = trc, tag

            s = Dseqrecord(Dseq(tag,trc))
            sl = s.rc().find_aa("YPYDVPDYA")

            assert str( s.rc()[sl].seq.translate()) == "YPYDVPDYA" 
            assert "YPYDVPDYA" in s.rc() 

            tag = "aaa tat cca tat gac gtt cca gac tat gca"
            trc = "ttt ata ggt ata ctg caa ggt ctg ata cgt"[::-1]

            s = Dseqrecord(Dseq(tag,trc, circular=True))
            sl = s.find_aa("YPYDVPDYA")
            assert str( s[sl].seq.translate() ) == "YPYDVPDYA" 
            assert "YPYDVPDYA" in s 


def test_map2():
    pCR_MCT1_HA46 = read("pCR_MCT1_HA46.gb")

    slc = pCR_MCT1_HA46.find_aa("VFFKE YPYDVPDYA IEG".replace(" ", ""))

    pCR_MCT1_HA46.map_target = slc

    map_ = pCR_MCT1_HA46.map_trace_files("*.ab1")

    assert set(map_)==set(['28-1rev_D04_026.ab1', '32-3rev_H04_018.ab1', '36-5rev_D05_041.ab1'])

    assert set([x.fname for x in pCR_MCT1_HA46.matching_reads]) == set(['28-1rev_D04_026.ab1', '32-3rev_H04_018.ab1', '36-5rev_D05_041.ab1'])
                                                                    
    assert set([x.fname for x in pCR_MCT1_HA46.not_matching_reads])==set(['02-G1_B01_013.ab1'])

    assert pCR_MCT1_HA46.find_aa("YPYDVPDYA".replace(" ", "")) == slice(1088, 1115, None)

    assert pCR_MCT1_HA46.find_aa("VFFKE YPYDVPDYA IEG".replace(" ", "")) == slice(1073, 1124, None)


if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])

