#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test empty
'''

import nose, sys

from nose.tools import assert_true, assert_false, assert_equal, assert_raises, assert_raises_regexp

from pydna import Dseqrecord, Dseq, parse, read,  eq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord as Srec


def test_linear_circular():
    ''' test Dseqrecord linear & circular property'''
    a=Dseqrecord("attt")
    a.stamp()
    assert_true( a.verify_stamp() )
    a=Dseqrecord("attt", linear = True)
    assert_true(a.linear     == True)
    assert_true(a.circular   == False)
    assert_true(a.rc().linear     == True)
    assert_true(a.rc().circular   == False)
    assert_true(a.seq.linear     == True)
    assert_true(a.seq.circular   == False)

    a=Dseqrecord("attt", linear = False)
    assert_true(a.linear     == False)
    assert_true(a.circular   == True)
    assert_true(a.rc().linear     == False)
    assert_true(a.rc().circular   == True)
    assert_true(a.seq.linear     == False)
    assert_true(a.seq.circular   == True)

    a=Dseqrecord("attt", circular = True)
    assert_true(a.linear     == False)
    assert_true(a.circular   == True)
    assert_true(a.rc().linear     == False)
    assert_true(a.rc().circular   == True)
    assert_true(a.seq.linear     == False)
    assert_true(a.seq.circular   == True)

    a=Dseqrecord("attt", circular = False)
    assert_true(a.linear     == True)
    assert_true(a.circular   == False)
    assert_true(a.rc().linear     == True)
    assert_true(a.rc().circular   == False)
    assert_true(a.seq.linear     == True)
    assert_true(a.seq.circular   == False)

def test_stamp():

    a=Dseqrecord("attt")
    a.stamp()
    assert_true( a.verify_stamp() )

def test_revcomp():

    a=Dseqrecord("attt")

    rc = a.rc()

    assert_true( str(rc.seq) == "aaat")

def test_initialization():
    a=[]

    a.append(       Dseqrecord("attt")          )
    a.append(       Dseqrecord(Dseq("attt"))    )
    a.append(       Dseqrecord(Seq("attt"))     )
    a.append(       Dseqrecord(Srec(Seq("attt"))))
    a.append(       Dseqrecord(Dseqrecord("attt")) )

    for b in a:
       assert_true( type(b.seq) == Dseq          )
       assert_true( str(b.seq.watson) == "attt" )
       assert_true( str(b.seq.crick)  == "aaat" )
       assert_true( str(b.seq) == "attt"        )
       assert_true( str(b.seq) == "attt"        )
       assert_true(b.linear     == b.seq.linear  )
       assert_true(b.linear     == True         )
       assert_true(b.circular   == False        )
       assert_true(b.seq.linear     == True     )
       assert_true(b.seq.circular   == False    )

    a=[]
    a.append(       Dseqrecord("attt", circular=True)           )
    a.append(       Dseqrecord(Dseq("attt"), circular=True)     )
    a.append(       Dseqrecord(Seq("attt"), circular=True)      )
    a.append(       Dseqrecord(Srec(Seq("attt")), circular=True))
    a.append(       Dseqrecord(Dseqrecord("attt"), circular=True  ))

    for b in a:
       assert_true( type(b.seq) == Dseq          )
       assert_true( str(b.seq.watson) == "attt" )
       assert_true( str(b.seq.crick)  == "aaat" )
       assert_true( str(b.seq) == "attt"        )
       assert_true( str(b.seq) == "attt"        )
       assert_true(b.linear     == b.seq.linear  )
       assert_true(b.linear     == False        )
       assert_true(b.circular   == True         )
       assert_true(b.seq.linear     == False    )
       assert_true(b.seq.circular   == True     )

    a=[]
    a.append(Dseqrecord(Dseq("attt",circular=True), circular=True))
    a.append(Dseqrecord(Dseq("attt",circular=False), circular=True))
    a.append(Dseqrecord(Dseq("attt",circular=True), circular=False))
    a.append(Dseqrecord(Dseq("attt",circular=False), circular=False))

    circular = [True,True,False,False]
    linear   = [False,False,True,True]

    for b,ci,li in zip(a,circular,linear):
       assert_true( type(b.seq) == Dseq          )
       assert_true( str(b.seq.watson) == "attt" )
       assert_true( str(b.seq.crick)  == "aaat" )
       assert_true( str(b.seq) == "attt"        )
       assert_true( str(b.seq) == "attt"        )
       assert_true(b.linear     == b.seq.linear  )
       assert_true(b.linear     == li  )
       assert_true(b.circular   == ci         )
       assert_true(b.seq.linear     == li    )
       assert_true(b.seq.circular   == ci     )

    a=[]
    ds = Dseq("attt", "taaa")
    assert_true(ds.linear == True)
    assert_true(ds.ovhg == -1)
    assert_true( str(ds.watson) == "attt" )
    assert_true( str(ds.crick) == "taaa" )

    #   attt
    #    aaat

    a.append(Dseqrecord(ds, circular = False))
    assert_true(ds.linear == True)
    a.append(Dseqrecord(ds, linear  = True))
    assert_true(ds.linear == True)

    a.append(Dseqrecord(ds, circular=True))
    assert_true(ds.linear == True)
    a.append(Dseqrecord(ds, linear=False))
    assert_true(ds.linear == True)

    circular = [False,False,True,True]
    linear   = [True,True,False,False]
    crick    = ["taaa","taaa","aaat","aaat"]
    sek = ["attta","attta","attt", "attt"]
    for b,ci,li,s,cri in zip(a,circular,linear, sek, crick):

       assert_true( type(b.seq) == Dseq          )
       assert_true( str(b.seq.watson) == "attt" )
       assert_true( str(b.seq.crick)  == cri     )
       assert_true( str(b.seq) == s              )

       assert_true(b.linear     == b.seq.linear  )
       assert_true(b.linear     == li  )
       assert_true(b.circular   == ci         )
       assert_true(b.seq.linear     == li    )
       assert_true(b.seq.circular   == ci     )

    a=[]
    ds = Dseq("attt", "caaa")
    assert_true(ds.linear == True)
    assert_true(ds.ovhg == -1)

    a.append(Dseqrecord(ds, circular=False))
    assert_true(ds.linear == True)
    a.append(Dseqrecord(ds, linear=True))
    assert_true(ds.linear == True)

    with assert_raises(TypeError):
        Dseqrecord(ds, circular=True)

    assert_true(ds.linear == True)

    with assert_raises(TypeError):
        Dseqrecord(ds, linear=False)

    assert_true(ds.linear == True)

    with assert_raises(TypeError):
        b = Dseqrecord([])

    with assert_raises(TypeError):
        b = Dseqrecord(("a",))

    with assert_raises(TypeError):
        b = Dseqrecord(0)

    from pydna import read

    input =   '''
            LOCUS       New_DNA                    4 bp ds-DNA     linear       30-MAR-2013
            DEFINITION  .
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            COMMENT
            COMMENT     ApEinfo:methylated:1
            FEATURES             Location/Qualifiers
                 misc_feature    2..3
                                 /label=NewFeature
                                 /ApEinfo_fwdcolor=cyan
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
            ORIGIN
                    1 acgt
            //
            '''
    a = read(input)

    assert_equal( a.features[0].extract(a).seq.watson, "CG")


    b = a+a

    for f in b.features:
        assert_equal( b.features[0].extract(a).seq.watson, "CG")

    feature = a.features[0]

    s = Dseq("agctt","agcta")
    #print s.fig()
    #Dseq(-6)
    # agctt
    #atcga
    b = Dseqrecord(s)
    b.features.append(feature)
    cb = Dseqrecord(b,circular=True)
    assert_equal(b.features[0].extract(b).seq.watson.lower(), cb.features[0].extract(b).seq.watson.lower() )
    assert_equal(b.features[0].extract(b).seq.crick.lower(),  cb.features[0].extract(b).seq.crick.lower() )

    s = Dseq("aagct","aagct")
    #print s.fig()
    #Dseq(-6)
    #aagct
    # tcgaa
    b = Dseqrecord(s)
    with assert_raises(TypeError):
        cb = Dseqrecord(b, circular=True)

    s = Dseq("agctt","agcta")
    #print s.fig()
    #Dseq(-6)
    # agcta
    #ttcga

    b = Dseqrecord(s)
    b.features.append(feature)
    cb = Dseqrecord(b,circular=True)
    assert_equal(b.features[0].extract(b).seq.watson.lower(), cb.features[0].extract(b).seq.watson.lower() )
    assert_equal(b.features[0].extract(b).seq.crick.lower(),  cb.features[0].extract(b).seq.crick.lower() )

def test_Dseqrecord_cutting_circular():

    from Bio.Restriction import BsaI, KpnI, Acc65I

    test = "aaaaaaGGTACCggtctcaaaa"

    for i in range(len(test)):
        nt = test[i:]+test[:i]

        d = Dseqrecord(nt, circular = True).cut(Acc65I)[0]
        assert_equal(d.seq.watson.upper(), "GTACCGGTCTCAAAAAAAAAAG")
        assert_equal(d.seq.crick.upper(),  "GTACCTTTTTTTTTTGAGACCG")
        assert_equal(d.seq.ovhg, -4)

        d = Dseqrecord(nt, circular = True).cut(KpnI)[0]
        assert_equal(d.seq.watson.upper(), "CGGTCTCAAAAAAAAAAGGTAC")
        assert_equal(d.seq.crick.upper(), "CTTTTTTTTTTGAGACCGGTAC")
        assert_equal(d.seq.ovhg, 4)

        d = Dseqrecord(nt, circular = True).cut(BsaI)[0]
        assert_equal(d.seq.watson.upper(), "AAAAAAAAAGGTACCGGTCTCA")
        assert_equal(d.seq.crick.upper(), "TTTTTGAGACCGGTACCTTTTT")
        assert_equal(d.seq.ovhg, -4)


def test_Dseq_cutting_adding():

    from Bio.Seq import Seq
    from Bio.Restriction import BamHI,EcoRI, PstI, EcoRV, SmaI
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
    from Bio.SeqUtils.CheckSum import seguid
    from pydna import Dseq


    a = Dseq('GGATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGGATCC',
             'CCTAGGagtagatgatagtagcatcgcatgactagataagacgacgagtagtagccatgagagatattaatatatatatacgcgcaCCTAGG'[::-1],
             linear=True,
             ovhg=0)

    b = a.cut(BamHI)[1]


    assert_equal( b.watson , "GATCCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtG")
    assert_equal( b.crick  , "GATCCacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaG")

    c = Dseq('nCTGCAGtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtGAATTCn',
             'nGACGTCagtagatgatagtagcatcgcatgactagataagacgacgagtagtagccatgagagatattaatatatatatacgcgcaCTTAAGn'[::-1],
             linear=True,
             ovhg=0)

    f,d,l = c.cut((EcoRI, PstI))

    assert_equal( d.watson  , "GtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtG")
    assert_equal( d.crick   , "AATTCacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaCTGCA")


    e =    Dseq("nGAATTCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtCTGCAGn",
                "nCTTAAGagtagatgatagtagcatcgcatgactagataagacgacgagtagtagccatgagagatattaatatatatatacgcgcaGACGTCn"[::-1],
                linear=True,
                ovhg=0)

    f = e.cut((EcoRI,PstI))[1]

    assert_equal( f.watson ,"AATTCtcatctactatcatcgtagcgtactgatctattctgctgctcatcatcggtactctctataattatatatatatgcgcgtCTGCA")
    assert_equal( f.crick  , "GacgcgcatatatatataattatagagagtaccgatgatgagcagcagaatagatcagtacgctacgatgatagtagatgaG")



    ''' blunt cloning '''


    pUC19 = read("tests/pUC19.gb")

    assert_false( pUC19.linear )

    assert_true( len(pUC19) == 2686 )
    assert_true( len(pUC19.seq.watson) == 2686 )
    assert_true( len(pUC19.seq.crick) == 2686 )

    assert_true( pUC19.seq.circular == True)
    assert_true( pUC19.seq.linear   == False)

    pUC19_SmaI = pUC19.cut(SmaI)
    assert_true( len(pUC19_SmaI) == 1)
    pUC19_SmaI = pUC19_SmaI.pop()


    assert_true( pUC19_SmaI.linear )
    assert_true( len(pUC19_SmaI) == 2686 )
    assert_true( pUC19_SmaI.linear )

    pUC19_SmaI_a = pUC19_SmaI.seq + a

    assert_true(  pUC19_SmaI_a.linear   )
    assert_false( pUC19_SmaI_a.circular )

    pUC19_SmaI_a=pUC19_SmaI_a.looped()
    assert_true( len(pUC19_SmaI_a) == 2778 )

    assert_true(  pUC19_SmaI_a.circular )
    assert_false( pUC19_SmaI_a.linear   )
    assert_true( eq(pUC19_SmaI_a, read("tests/pUC19-SmaI-a.gb")   ))

    ''' sticky end cloning '''

    pUC19_BamHI = pUC19.cut(BamHI)

    assert_true( len(pUC19_BamHI) == 1)

    pUC19_BamHI = pUC19_BamHI.pop().seq

    assert_true( len(pUC19_BamHI.watson) == len(pUC19_BamHI.crick) == 2686 )

    pUC19_BamHI_a = pUC19_BamHI+b

    assert_true( len(pUC19_BamHI_a.watson) == len(pUC19_BamHI_a.crick) == 2772 )

    assert_true( pUC19_BamHI_a.circular == False)
    assert_true( pUC19_BamHI_a.linear   == True)

    pUC19_BamHI_a = pUC19_BamHI_a.looped()

    assert_true( pUC19_BamHI_a.circular == True)
    assert_true( pUC19_BamHI_a.linear   == False)

    assert_true( eq(pUC19_BamHI_a, read("tests/pUC19-BamHI-a.gb")))

    pUC19_BamHI_a_rc = pUC19_BamHI+b.rc()

    pUC19_BamHI_a_rc = pUC19_BamHI_a_rc.looped()


    assert_true( pUC19_BamHI_a.circular == True)
    assert_true( pUC19_BamHI_a.linear   == False)
    assert_true( eq(pUC19_BamHI_a_rc, read("tests/pUC19-BamHI-a-rc.gb")))

    ''' adding (ligating) dsDNA objects '''
    with assert_raises_regexp(TypeError, "circular"):
        pUC19+a
    with assert_raises_regexp(TypeError, "circular"):
        a+pUC19
    with assert_raises_regexp(TypeError, "compatible"):
        a+b
    with assert_raises_regexp(TypeError, "compatible"):
        b+a
    with assert_raises_regexp(TypeError, "compatible"):
        d+d

    ''' directional cloning '''

    pUC19_EcoRI_PstI = pUC19.cut(EcoRI, PstI).pop(0)

    with assert_raises_regexp(TypeError, "compatible"):
        pUC19_EcoRI_PstI + d

    pUC19_EcoRI_PstI_d = pUC19_EcoRI_PstI + d.rc()

    pUC19_EcoRI_PstI_d =  pUC19_EcoRI_PstI_d.looped()

    assert_true( eq(pUC19_EcoRI_PstI_d,      read("tests/pUC19-EcoRI_PstI-d-rc.gb")))
    assert_true( eq(pUC19_EcoRI_PstI_d.rc(), read("tests/pUC19-EcoRI_PstI-d-rc.gb")))


def test_Dseqrecord_cutting_adding():
    from Bio.Restriction import Bsu36I, BstAPI
    pCAPs = read("tests/pCAPs.gb")
    a,b = pCAPs.cut(Bsu36I, BstAPI)
    c=(a+b).looped()
    assert_true( eq(c, pCAPs) )


    a = (Dseqrecord( Dseq(  'AATTCACANGGTACCNGGTACCNGCGGATATC',
                             'GTGTNCCATGGNCCATGGNCGCCTATAG'[::-1], -4)),

         Dseqrecord( Dseq(      'CACANGGTACCNGGTACCNGCGGATATC',
                             'GTGTNCCATGGNCCATGGNCGCCTATAG'[::-1], 0)),

         Dseqrecord( Dseq(    'CACANGGTACCNGGTACCNGCGGATATC',
                       'AATTGTGTNCCATGGNCCATGGNCGCCTATAG'[::-1], 4)),)

    from Bio.Restriction import KpnI, Acc65I, NlaIV

    enzymes = [Acc65I, NlaIV, KpnI]

    for enz in enzymes:
        for f in a:
            b,c,d = f.cut(enz)
            e=b+c+d
            assert str(e.seq).lower() == str(f.seq).lower()



    #from pydna import *
    #from pydna_helper import gb, ape
    from Bio.Restriction import KpnI, BamHI, Acc65I, NlaIV, EcoRI, EcoRV

    a=read('''

LOCUS       New_DNA                   10 bp ds-DNA     linear       02-APR-2013
DEFINITION
ACCESSION   New_DNA
VERSION     New_DNA
KEYWORDS    .
SOURCE
  ORGANISM  . .
COMMENT
COMMENT     ApEinfo:methylated:1
FEATURES             Location/Qualifiers
     misc_feature    1..1
                     /label=1
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    2..2
                     /label=2
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    3..3
                     /label=3
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    4..4
                     /label=4
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    5..5
                     /label=5
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    6..6
                     /label=6
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    7..7
                     /label=7
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    8..8
                     /label=8
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    9..9
                     /label=9
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    10..10
                     /label=10
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
ORIGIN
        1 ttGGTACCgg
//''')

    b,c = a.cut(Acc65I)

    assert_equal( [f.qualifiers["label"] for f in b.features], [['1'], ['2'], ['3'], ['4'], ['5'], ['6'], ['7']])
    assert_equal( [f.qualifiers["label"] for f in c.features], [['4'], ['5'], ['6'], ['7'], ['8'], ['9'], ['10']])




    a=read('''

LOCUS       New_DNA                   33 bp ds-DNA     linear       08-NOV-2012
DEFINITION  .
ACCESSION
VERSION
SOURCE      .
  ORGANISM  .
COMMENT
COMMENT     ApEinfo:methylated:1
FEATURES             Location/Qualifiers
     misc_feature    1..11
                     /label=Acc65I-1
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    12..18
                     /label=Acc65I-2
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    19..33
                     /label=Acc65I-3
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    1..15
                     /label=KpnI-1
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    16..22
                     /label=KpnI-2
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    23..33
                     /label=KpnI-3
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    1..13
                     /label=NlaIV-1
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    14..20
                     /label=NlaIV-2
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
     misc_feature    21..33
                     /label=NlaIV-3
                     /ApEinfo_fwdcolor=cyan
                     /ApEinfo_revcolor=green
                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                     width 5 offset 0
ORIGIN
        1 GAATTCacan ggtaccnGGT ACCngcgGAT ATC
//

    ''')

    assert_true( a.seguid()=="di3hL8t2G4iQQsxlm_CtvnUMBz8" )

    assert_true( ([x.qualifiers["label"][0] for x in a.features] ==
    ['Acc65I-1', 'Acc65I-2', 'Acc65I-3', 'KpnI-1', 'KpnI-2',
     'KpnI-3', 'NlaIV-1', 'NlaIV-2', 'NlaIV-3']))

    b,c,d = a.cut(Acc65I)

    assert_true( [x.qualifiers["label"][0] for x in b.features] == ['Acc65I-1', 'KpnI-1', 'NlaIV-1'])
    assert_true( [x.qualifiers["label"][0] for x in c.features] == ['Acc65I-2', 'KpnI-2', 'NlaIV-2'])
    assert_true( [x.qualifiers["label"][0] for x in d.features] == ['Acc65I-3', 'KpnI-3', 'NlaIV-3'])
    e = b+c+d
    assert_true( sorted([x.qualifiers["label"][0] for x in e.features])  == [x.qualifiers["label"][0] for x in a.features])
    assert_true( str(a.seq)==str(e.seq))

    b,c,d = a.cut(KpnI)
    assert_true( [x.qualifiers["label"][0] for x in b.features] == ['Acc65I-1', 'KpnI-1', 'NlaIV-1'])
    assert_true( [x.qualifiers["label"][0] for x in c.features] == ['Acc65I-2', 'KpnI-2', 'NlaIV-2'])
    assert_true( [x.qualifiers["label"][0] for x in d.features] == ['Acc65I-3', 'KpnI-3', 'NlaIV-3'])
    e = b+c+d
    assert_true( sorted([x.qualifiers["label"][0] for x in e.features])  == [x.qualifiers["label"][0] for x in a.features])

    b,c,d = a.cut(NlaIV)
    assert_true( [x.qualifiers["label"][0] for x in b.features] == ['Acc65I-1', 'NlaIV-1'])
    assert_true( [x.qualifiers["label"][0] for x in c.features] == ['NlaIV-2'])
    assert_true( [x.qualifiers["label"][0] for x in d.features] == [ 'KpnI-3', 'NlaIV-3'])
    e = b+c+d
    assert_true( str(a.seq)==str(e.seq))

    b,c = a.cut(EcoRI)
    e = b+c
    assert_true( str(a.seq)==str(e.seq))

    b,c = a.cut(EcoRV)
    e = b+c
    assert_true( str(a.seq)==str(e.seq))

    b,c,d = a.cut(EcoRI,EcoRV)
    e = b+c+d

    assert_true( str(a.seq)==str(e.seq))

    b,c,d, f = a.cut(Acc65I,EcoRI)
    e = b+c+d+f
    assert_true( str(a.seq)==str(e.seq))

    b,c,d, f = a.cut(EcoRI,Acc65I)
    e = b+c+d+f
    assert_true( str(a.seq)==str(e.seq))




def test_Dseq_slicing():
    from Bio.Restriction import BamHI
    a=Dseq("ggatcc","ggatcc",0)

    assert_equal( a[:].watson, a.watson )
    assert_equal( a[:].crick, a.crick )
    assert_equal( a.ovhg, a[:].ovhg )

    b,c = a.cut(BamHI)
    d = b[1:5]
    e = d.rc()
    assert_equal( d+e, Dseq("gatc","gatc",0))



def test_Dseq_slicing2():

    from Bio.Restriction import BamHI, EcoRI, KpnI

    a = Dseq("aaGGATCCnnnnnnnnnGAATTCccc", circular=True)

    assert_equal( a.cut(EcoRI, BamHI, KpnI,), a.cut(BamHI, EcoRI, KpnI,))

    a = Dseqrecord("aaGGATCCnnnnnnnnnGAATTCccc", circular=True)

    assert_equal( a.cut( EcoRI, BamHI,  KpnI,)[0].seq, a.cut( BamHI, EcoRI,  KpnI,)[0].seq)



def test_rogerstager():

    from Bio.Seq import Seq
    from Bio.Restriction import BsaI

    answ = []
    answ.append( Dseq('aaaaaaaaaaaaggtctca', 'ttttttttccagagttttt'[::-1]) )
    answ.append( Dseq('aaaaaaaaaggtctca', 'tttttccagagttttt'[::-1]) )

    tests = [Seq("aaaaaaggtctcaaaaaaa"),
             Seq("aaaaaaggtctcaaaa")  ]
    for s in tests:
        d = Dseqrecord(s).looped()
        for f in d.cut(BsaI):
            a = answ.pop(0)
            assert_true( f.seq.watson == a.watson )
            assert_true( f.seq.crick  == a.crick )
            assert_true( f.seq.ovhg == a.ovhg )
            assert_true( eq(f.seq, a) )



def test_cut_around_and_religate():

    from Bio.Restriction import KpnI, BamHI, Acc65I

    def cut_and_religate_Dseq(seq_string, enz, top):
        ds = Dseq(seq_string, linear=top)
        frags = ds.cut(enz)
        if not frags:
            return
        a = frags.pop(0)
        for f in frags:
            a+=f
        if not top:
            a=a.looped()
        assert_true( eq(a,ds) )

    seqs = [
            ("aaaGGTACCcccGGTACCCgggGGTACCttt", BamHI,  False),
            ("aaaGGTACCcccGGTACCCgggGGTACCttt", Acc65I, False),
            ("aaaGGTACCcccGGTACCCgggGGTACCttt", KpnI,   False),

            ("aaaGGTACCcccGGTACCCgggGGTACCttt", Acc65I, True),
            ("aaaGGTACCcccGGTACCCgggGGTACCttt", KpnI,   True),
            ("aaaGGTACCcccGGTACCCgggGGTACCttt", BamHI,  True),

            ("aaaGGTACCcccGGATCCCgggGGTACCttt", [Acc65I, BamHI], False),
            ("aaaGGTACCcccGGATCCCgggGGTACCttt", [KpnI, BamHI],   False),
            ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI,Acc65I],  False),
            ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI,KpnI],    False),

            ("aaaGGTACCcccGGATCCCgggGGTACCttt", [Acc65I, BamHI], True),
            ("aaaGGTACCcccGGATCCCgggGGTACCttt", [KpnI, BamHI],   True),
            ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI,Acc65I],  True),
            ("aaaGGTACCcccGGATCCCgggGGTACCttt", [BamHI,KpnI],    True),

            ]

    for s in seqs:
        sek, enz, lin = s
        for i in range(len(sek)):
            zek = sek[i:]+sek[:i]
            cut_and_religate_Dseq(zek, enz, lin)


def test_features_change_ori():

    s = read('''
                LOCUS       New_DNA                   13 bp ds-DNA     circular     12-NOV-2013
                DEFINITION  .
                ACCESSION
                VERSION
                SOURCE      .
                ORGANISM  .
                COMMENT
                COMMENT     ApEinfo:methylated:1
                FEATURES             Location/Qualifiers
                     misc_feature    join(9..10,12..13,1..1,3..6)
                                     /label=hej
                                     /ApEinfo_fwdcolor=cyan
                                     /ApEinfo_revcolor=green
                                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                     width 5 offset 0
                ORIGIN
                        1 gattttaatc acc
                //''')

    #from pydna_helper import ape

    for i in range(1, len(s)):
        b=s.shifted(i)
        assert_true(  str(b.features[0].extract(b).seq).lower()=="tcccgtttt")

    s = read('''
            LOCUS       New_DNA                   21 bp ds-DNA     circular     03-APR-2013
            DEFINITION  a
            ACCESSION
            VERSION
            SOURCE      .
              ORGANISM  .
            COMMENT
            COMMENT     ApEinfo:methylated:1
            FEATURES             Location/Qualifiers
                 misc_feature    join(18..21,1..4)
                                 /label=bb
                                 /ApEinfo_fwdcolor=cyan
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
                 misc_feature    5..17
                                 /label=ins
                                 /ApEinfo_fwdcolor=#e03c2b
                                 /ApEinfo_revcolor=green
                                 /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}
                                 width 5 offset 0
            ORIGIN
                    1 aaaGGTACCt ttGGATCCggg
            //
        ''')

    assert_true( str(s.features[0].extract(s).seq)  == "CGGGAAAG" )
    assert_true( str(s.features[1].extract(s).seq)  == "GTACCTTTGGATC" )

    for i in range(1, len(s)):

        b = s.shifted(i)
        assert_true( [str(f.extract(b).seq) for f in b.features if f.qualifiers["label"][0]=='ins'][0] == "GTACCTTTGGATC" )
        assert_true( [str(f.extract(b).seq) for f in b.features if f.qualifiers["label"][0]=='bb'][0] == "CGGGAAAG" )

    from Bio.Restriction import Acc65I,KpnI, BamHI

    bb1, ins1 = sorted(s.cut(Acc65I, BamHI), key=len, reverse=True)

    for i in range(1, len(s)):
        b = s.shifted(i)

        bb, ins = sorted(b.cut(Acc65I, BamHI), key=len, reverse=True)

        assert_true( eq(bb1, bb) )
        assert_true( eq(ins1,ins) )

        assert_true( bb.features[0].extract(bb).seq.watson == "CGGGAAAG" )
        assert_true( bb.features[0].extract(bb).seq.crick  == "CTTTCCCG" )

        assert_true( eq(bb.features[0].extract(bb), s.features[0].extract(s) ) )

        assert_true( ins.features[0].extract(ins).seq.watson == "GTACCTTTG" )
        assert_true( ins.features[0].extract(ins).seq.crick  == "GATCCAAAG" )

        assert_true( str(ins.features[0].extract(ins).seq) == str(s.features[1].extract(s).seq) )


def test_synced():
    pUC19        = read("tests/pUC19.gb")
    pUC19_LAC4   = read("tests/pUC_LAC4.gb")
    pUC19_LAC4_c = read("tests/pUC_LAC4_correct_rotation.gb")

    correct = str(pUC19_LAC4_c.seq).upper()

    for i in range(1, len(pUC19_LAC4), 500):
        cand = pUC19_LAC4.shifted(i)
        assert_equal(str(cand.synced("tcgcgcgtttcggtgatgacggtga").seq).upper(),
                         correct,
                         str(pUC19_LAC4.synced(pUC19).seq).upper())
        print(i, end=' ')
    print()



def test_synced3():
    pGUP1 = read("tests/pGUP1_correct.gb")
    pGREG505 = read("tests/pGREG505.gb")
    pGUP1_not_synced =  read("tests/pGUP1_not_synced.gb")
    assert_equal(pGUP1_not_synced.synced(pGREG505).seguid(), '42wIByERn2kSe_Exn405RYwhffU')


def test_synced2():
    pUC19            = read("tests/pUC19.gb")
    pUC19_small_gene = read("tests/pUC19_small_gene.gb")

    correct = str(pUC19_small_gene.seq).upper()

    for i in range(1, len(pUC19_small_gene), 500):
        cand = pUC19_small_gene.shifted(i)
        assert_equal(str(cand.synced("tcgcgcgtttcggtgatgacggtga").seq).upper(),
                         correct,
                         str(cand.synced(pUC19).seq).upper())
        print(i, end=' ')
    print()

if __name__ == '__main__':
    nose.runmodule(argv=[sys.argv[0], '--nocapture'])

