#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
utils tests
'''

import unittest

from pydna import eq, shift_origin, read, Dseqrecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class test_utils(unittest.TestCase):

    def test_eq(self):
        ''' test eq'''
        self.assertTrue(  eq( "AAA" ,"TTT", linear   = True ) )
        self.assertTrue(  eq( "AAA" ,"TTT", linear   = False) )

        self.assertTrue(  eq( "aAA" ,"TtT", linear   = True ) )
        self.assertTrue(  eq( "AAa" ,"TtT", linear   = False) )


        self.assertTrue(  eq( "ATA" ,"AAT", circular = True ) )
        self.assertFalse( eq( "ATA" ,"AAT", circular = False) )
        self.assertTrue(  eq( "AAA" ,"AAA", linear   = True ) )
        self.assertTrue(  eq( "AAA" ,"AAA", linear   = False) )

        self.assertTrue(  eq( "ATA" ,Seq("AAT"), circular = True ) )
        self.assertFalse( eq( "ATA" ,Seq("AAT"), circular = False) )
        self.assertTrue(  eq( "AAA" ,Seq("AAA"), linear   = True ) )
        self.assertTrue(  eq( "AAA" ,Seq("AAA"), linear   = False) )

        self.assertTrue(  eq( "ATA" ,SeqRecord("AAT"), circular = True ) )
        self.assertFalse( eq( "ATA" ,SeqRecord("AAT"), circular = False) )
        self.assertTrue(  eq( "AAA" ,SeqRecord("AAA"), linear   = True ) )
        self.assertTrue(  eq( "AAA" ,SeqRecord("AAA"), linear   = False) )

        self.assertTrue(  eq( "ATA" ,Dseqrecord("AAT"), circular = True ) )
        self.assertFalse( eq( "ATA" ,Dseqrecord("AAT"), circular = False) )
        self.assertTrue(  eq( "AAA" ,Dseqrecord("AAA"), linear   = True ) )
        self.assertTrue(  eq( "AAA" ,Dseqrecord("AAA"), linear   = False) )

        self.assertTrue(  eq( Seq("ATA") ,SeqRecord("AAT"), circular = True ) )
        self.assertFalse( eq( Seq("ATA") ,SeqRecord("AAT"), circular = False) )
        self.assertTrue(  eq( Seq("AAA") ,SeqRecord("AAA"), linear   = True ) )
        self.assertTrue(  eq( Seq("AAA") ,SeqRecord("AAA"), linear   = False) )

        self.assertTrue(  eq( Seq("ATA") ,Dseqrecord("AAT"), circular = True ) )
        self.assertFalse( eq( Seq("ATA") ,Dseqrecord("AAT"), circular = False) )
        self.assertTrue(  eq( Seq("AAA") ,Dseqrecord("AAA"), linear   = True ) )
        self.assertTrue(  eq( Seq("AAA") ,Dseqrecord("AAA"), linear   = False) )

        self.assertTrue(  eq( Dseqrecord("AAA",circular=False) ,Dseqrecord("AAA",circular=False)) )
        self.assertTrue(  eq( Dseqrecord("AAA",circular=True)  ,Dseqrecord("AAA",circular=True))   )
        self.assertFalse( eq( Dseqrecord("ATA",circular=False) ,Dseqrecord("AAT",circular=False)) )
        self.assertTrue(  eq( Dseqrecord("ATA",circular=True)  ,Dseqrecord("AAT",circular=True)) )

        self.assertEqual( True , True )
        self.assertTrue( True )
        self.assertFalse( False )

        with self.assertRaises(TypeError):
            1+"1"

        with self.assertRaisesRegexp(TypeError, "unsupported"):
            1+"1"

    def test_shift_origin(self):

        pCAPs   = read("./pCAPs.gb")
        self.assertTrue( pCAPs.circular )
        pCAPs_b = shift_origin(pCAPs, 200)
        self.assertEqual( len(pCAPs), len(pCAPs_b) )
        self.assertTrue( pCAPs_b.circular )
        self.assertTrue( eq(pCAPs, pCAPs_b) )
        pCAPs_b_linear = pCAPs_b.tolinear()
        self.assertTrue( eq(pCAPs, pCAPs_b_linear, circular=True) )
        pCAPs_c = pCAPs[200:]+pCAPs[:200]
        self.assertTrue( eq(pCAPs, pCAPs_c, circular=True) )
        with self.assertRaisesRegexp(ValueError, "shift"):
            pCAPs_b = shift_origin(pCAPs, 20000)


    def test_copy_features(self):

        from pydna.utils import seguid
        from pydna import read,copy_features
        a=read("./pCAPs.gb")
        b=read("./pCAPs_fasta.txt")

        for sh in [1,2,3,3127,3128,3129]:
            newb = (b[sh:]+b[:sh]).looped()
            copy_features(a, newb)
            #print "a",[len(str(f.extract(a).seq.lower()) for f in a.features if len(f)>10]
            #print "b",[len(str(f.extract(newb).seq).lower()) for f in newb.features]

            self.assertTrue( sorted([str(f.extract(a).seq).lower() for f in a.features if len(f)>10],key=len)
                            == sorted([str(f.extract(newb).seq).lower() for f in newb.features],key=len))

        b=b.rc()

        for sh in [1,2,3,3127,3128,3129]:
            newb = b[sh:]+b[:sh]
            copy_features(a, newb)
            self.assertTrue( sorted([str(f.extract(a).seq).lower() for f in a.features if len(f)>10],key=len) == sorted([str(f.extract(newb).seq).lower() for f in newb.features],key=len))

        seguid_bla = "riT98j2v4NxVS8sbw_Q8epCwQwo"
        seguid_cre = "xLZ2xs2O8CUMmWh2OrhmNFp5ZLg"

        copy_features(a, b)
        assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre, seguid_cre, seguid_bla, seguid_bla]

        b=read("./pCAPs_fasta.txt").looped()

        b=b.synced("attaacgagtgccgtaaacgacgatggttttacc")

        copy_features(a, b)
        assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]

        b=read("./pCAPs_fasta.txt").looped()
        b=b.synced("ttaacgagtgccgtaaacgacgatggttttacc")

        copy_features(a, b)
        assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]

        b=read("./pCAPs_fasta.txt").looped()
        b=b.synced("taacgagtgccgtaaacgacgatggttttacc")

        copy_features(a, b)
        assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_bla,seguid_bla]

        b=read("./pCAPs_fasta.txt").looped()
        b=b.synced("gttaccaatgcttaatcagtgaggcacctatctcagc")

        copy_features(a, b)
        assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]

        b=read("./pCAPs_fasta.txt").looped()
        b=b.synced("ttaccaatgcttaatcagtgaggcacctatctcagc")

        copy_features(a, b)
        assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,seguid_bla,seguid_bla]

        b=read("./pCAPs_fasta.txt").looped()
        b=b.synced("taccaatgcttaatcagtgaggcacctatctcagc")

        copy_features(a, b)
        assert [seguid(f.extract(b).seq) for f in b.features] == [seguid_cre,seguid_cre,]


    def test_cseguid(self):
        from pydna.utils import cseguid
        x="tcgcgcgtttcggtgatgacggtgaaaacctctgacacatgcagctcccggagacggtcacagcttgtctgtaagcggatgccgggagcagacaagcccgtcagggcgcgtcagcgggtgttggcgggtgtcggggctggcttaactatgcggcatcagagcagattgtactgagagtgcaccatatgcggtgtgaaataccgcacagatgcgtaaggagaaaataccgcatcaggcgccattcgccattcaggctgcgcaactgttgggaagggcgatcggtgcgggcctcttcgctattacgccagctggcgaaagggggatgtgctgcaaggcgattaagttgggtaacgccagggttttcccagtcacgacgttgtaaaacgacggccagtgaattcgagctcggtacccgggGATCTATGAATATGGATCCGACTTACTGCAGGAATTCAAGCTACTGTTAGAgatcctctagagtcgacctgcaggcatgcaagcttggcgtaatcatggtcatagctgtttcctgtgtgaaattgttatccgctcacaattccacacaacatacgagccggaagcataaagtgtaaagcctggggtgcctaatgagtgagctaactcacattaattgcgttgcgctcactgcccgctttccagtcgggaaacctgtcgtgccagctgcattaatgaatcggccaacgcgcggggagaggcggtttgcgtattgggcgctcttccgcttcctcgctcactgactcgctgcgctcggtcgttcggctgcggcgagcggtatcagctcactcaaaggcggtaatacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgagattatcaaaaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttaccaatgcttaatcagtgaggcacctatctcagcgatctgtctatttcgttcatccatagttgcctgactccccgtcgtgtagataactacgatacgggagggcttaccatctggccccagtgctgcaatgataccgcgagacccacgctcaccggctccagatttatcagcaataaaccagccagccggaagggccgagcgcagaagtggtcctgcaactttatccgcctccatccagtctattaattgttgccgggaagctagagtaagtagttcgccagttaatagtttgcgcaacgttgttgccattgctacaggcatcgtggtgtcacgctcgtcgtttggtatggcttcattcagctccggttcccaacgatcaaggcgagttacatgatcccccatgttgtgcaaaaaagcggttagctccttcggtcctccgatcgttgtcagaagtaagttggccgcagtgttatcactcatggttatggcagcactgcataattctcttactgtcatgccatccgtaagatgcttttctgtgactggtgagtactcaaccaagtcattctgagaatagtgtatgcggcgaccgagttgctcttgcccggcgtcaatacgggataataccgcgccacatagcagaactttaaaagtgctcatcattggaaaacgttcttcggggcgaaaactctcaaggatcttaccgctgttgagatccagttcgatgtaacccactcgtgcacccaactgatcttcagcatcttttactttcaccagcgtttctgggtgagcaaaaacaggaaggcaaaatgccgcaaaaaagggaataagggcgacacggaaatgttgaatactcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccacctgacgtctaagaaaccattattatcatgacattaacctataaaaataggcgtatcacgaggccctttcgtc"
        assert cseguid(x) == cseguid(x.upper()) == cseguid(x.lower()) == 'JgiKgBksd2v3q99NStKLepoQCm8'
        from Bio.SeqUtils.CheckSum  import seguid as base64_seguid



if __name__ == '__main__':
    unittest.main()
