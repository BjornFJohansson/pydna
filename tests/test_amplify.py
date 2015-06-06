#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test parse
'''

import unittest
from pydna import parse
from pydna import pcr
from Bio.SeqUtils.CheckSum import seguid

class test_pcr(unittest.TestCase):

    def test_circ_pcr(self):
        ''' test circ pcr'''
        s = parse('''
        >MCT4_flaghis_rv
        gccgcaagcttgtcgtcatcgtctttgtagtcCATggccaggagggttggttccgcc
        >MCT4_flaghis_fw
        ccaccaccaccaccaccaccaccaccaccacggaggggccgtggtggacgagggcc''', ds=False)

        t = parse('''
        >hej circular
        TCAGGTGAGGCGGAACCAACCCTCCTGGCCATGggaggggccgtggtggacgagggccccacaggcg
        ''')

        p = pcr(s,t)

        self.assertTrue( p.seguid() == 'wYEpbWxO8tfgA_AL67wqGzYmc5k')



    def test_pcr(self):
        ''' test pcr'''

        raw=[]

        raw.append(

        ('7JOV1MJBZJp2Smja/7KFGhS2SWY',

        parse('''
        >524_pFA6aF (29-mer)
        cacatacgatttaggtgacactatagaac

        >523_AgTEF1tpR (21-mer)
        ggttgtttatgttcggatgtg
        '''),

        parse("./pAG25.gb"),)
        )

        raw.append(

        ('7pPxy/bQvs4+7CaOgiywQVzUFDc',

        parse('''
        >lowgc_f
        TTTCACTAGTTACTTGTAGTCGacgtgccatctgtgcagacaaacgcatcaggatat

        >lowgc_r
        AAGTTGGAAATCTAGCTTTTCTTgacgtcagcggccgcattgcaca
        '''),
        parse("./pCAPs.gb"),)
        )

        raw.append(

        ('7JOV1MJBZJp2Smja/7KFGhS2SWY',

        parse('''
        >524_pFA6aF (29-mer)
        cacatacgatttaggtgacactatagaac

        >523_AgTEF1tpR (21-mer)
        ggttgtttatgttcggatgtg
        '''),

        parse("../tests/pAG25.gb"),

        ))

        raw.append(

        ('yshvYTXr9iXCnh3YytWQRDBNQzI',

        parse('''
        >ForwardPrimer1
        gctactacacacgtactgactg

        >ReversePrimer
        tgtggttactgactctatcttg

        LOCUS       sequence_50_bp            46 bp    DNA     circular UNK 08-FEB-2013
        DEFINITION  sequence_50_bp circular
        ACCESSION   sequence_50_bp
        VERSION     sequence_50_bp
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 ccaagataga gtcagtaacc acagctacta cacacgtact gactgt
        //

        '''),))

        raw.append(

        ('yshvYTXr9iXCnh3YytWQRDBNQzI',

        parse('''
        >ForwardPrimer2
        gctactacacacgtactgactg

        >ReversePrimer
        tgtggttactgactctatcttg

        LOCUS       template                  46 bp    DNA     circular UNK 15-OCT-2012
        DEFINITION  template circular
        ACCESSION   template
        VERSION     template
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 ccaagataga gtcagtaacc acagctacta cacacgtact gactgt
        //
        '''),))
        raw.append(

        ('yshvYTXr9iXCnh3YytWQRDBNQzI',

        parse('''
        >ForwardPrimer3
        gctactacacacgtactgactg

        >ReversePrimer
        tgtggttactgactctatcttg

        LOCUS       template                  46 bp    DNA     circular UNK 08-FEB-2013
        DEFINITION  template circular
        ACCESSION   template
        VERSION     template
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 tccaagatag agtcagtaac cacagctact acacacgtac tgactg
        //
        '''),))

        raw.append(
        ('yshvYTXr9iXCnh3YytWQRDBNQzI',

        parse('''
        >ForwardPrimer4
        gctactacacacgtactgactg

        >ReversePrimer
        tgtggttactgactctatcttg

        LOCUS       template                  46 bp    DNA     circular UNK 15-OCT-2012
        DEFINITION  template circular
        ACCESSION   template
        VERSION     template
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 gtccaagata gagtcagtaa ccacagctac tacacacgta ctgact
        //
        '''),))

        raw.append(
        ('60meNXeGKO7ahZwcIl5yXHFC3Yg',
        parse('''
        >fw1
        cacatacgatttaggtgacactatagaac
        >rv
        ggttgtttatgttcggatgtg

        LOCUS       tm                        50 bp    DNA     circular UNK 15-OCT-2012
        DEFINITION  tm circular
        ACCESSION   tm
        VERSION     tm
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 cacatccgaa cataaacaac ccacatacga tttaggtgac actatagaac
        //
        '''),)
        )

        raw.append(

        ('60meNXeGKO7ahZwcIl5yXHFC3Yg',
        parse('''
        >fw2
        cacatacgatttaggtgacactatagaac
        >rv
        ggttgtttatgttcggatgtg

        LOCUS       tm                        50 bp    DNA     circular UNK 15-OCT-2012
        DEFINITION  tm circular
        ACCESSION   tm
        VERSION     tm
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 acatccgaac ataaacaacc cacatacgat ttaggtgaca ctatagaacc
        //

        '''),)
        )

        raw.append(

        ('60meNXeGKO7ahZwcIl5yXHFC3Yg',
        parse('''
        >fw3
        cacatacgatttaggtgacactatagaac
        >rv
        ggttgtttatgttcggatgtg
        LOCUS       tm                        50 bp    DNA     circular UNK 15-OCT-2012
        DEFINITION  tm circular
        ACCESSION   tm
        VERSION     tm
        KEYWORDS    .
        SOURCE      .
          ORGANISM  .
                    .
        FEATURES             Location/Qualifiers
        ORIGIN
                1 ccacatccga acataaacaa cccacatacg atttaggtga cactatagaa
        //

        '''),)
        )


        raw.append(
        ('y6ohCJ4O+8Is012DItz4F4saxNo',
        parse('''
        >f_Eric_Ma
        ARATGAGTCTTCTRACCGAGGTCG
        >r_Eric_Ma
        TGAAAAGACATCYTCAAGYYTCTG
        >templ
        AGCAAAAGCAGGTAGATATTGAAAAATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCGCAGAGACTTGAAGATGTCTCTGCAGGGAAGAACACTGATCTCGAGGCTCTCATGGAATGGCTAAAGACAAGACCAATCCTGTCACCTCTGACTAAGGGGATTTTAGGGTTTGTGTTCACGCTCACCGTGCCCAGTGAGCGAGGACTGCAGCGTAGACGCTTTGTCCAGAATGCCTTAAATGGGAATGGAGACCCAAACAACATGGACAGGGCAGTCAAACTATACAGGAAGCTGAAAAGAGAGATAACATTCCATGGGGCTAAAGAGGTTGCACTCAGCTATTCAACCGGTGCACTTGCCAGTTGCATGGGTCTCATATACAACAGGATGGGAACGGTAACCACAGAAGTAGCTTTTGGCCTGGTGTGTGCCACTTGTGAGCAGATTGCTGACTCACAGCATCGATCTCACAGACAGATGGTGACTACCACCAACCCACTAATCAGGCATGAAAACAGAATGGTGCTGGCCAGCACTACAGCTAAGGCTATGGAGCAGATGGCTGGATCGAGTGAACAGGCAGCGGAAGCCATGGAGGTTGCTAGTCAGGCTAGGCAGATGGTGCAGGCAATGAGGACAATTGGGACTCACCCTAGCTCCAGTGCCGGTCTGAAAGATGATCTTCTTGAAAATTTGCAGGCCTACCAGAAGCGGATGGGAGTGCAAATGCAGCGATTCAAGTGATCCTCTCGTTATTGCCGCAAGTATCATTGGGATCTTGCACTTGATATTGTGGATTCTTGATCGTCCTTTCTTCAAATGTATTTATCGTCGCCTTAAATACGGTTTGAAAAGAGGGCCTTCTACGGAAGGAGTGCCTGAGTCTATGAGGGAAGAGTATCGGCAGGAACAGCAGAGTGCTGTGGATGTTGACGATGGTCATTTTGTCAACATAGAGCTGGAGTAAAAAACTACCTTGTTTCTACT
        '''),)
        )

        for key, tst in enumerate(raw):
            print key
            print tst[1][0].name
            print "pcr test", key
            print
            self.assertEqual(tst[0], seguid(pcr(tst[1:]).seq))

if __name__ == '__main__':
    unittest.main()






'''
ARATGAGTCTTCTRACCGAGGTCG
||||||||||||||||||||||||
AAATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCGCAGAGACTTGAAGATGTCT
                                                                             |||||||||||||||||||
                                                                             CAGARRCTTGARGATGTCTTTTCA


'''


