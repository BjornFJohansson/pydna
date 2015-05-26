#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test parse
'''

import unittest

from pydna import parse, read

class test_parse(unittest.TestCase):

    def test_parse1(self):
        ''' test parsing fasta sequences from a text'''

        text   =  '''
                points....: 1

                The sequence seq below represents a double stranded linear DNA molecule.

                >seq
                CTCCCCTATCACCAGGGTACCGATAGCCACGAATCT

                Give the sequence(s) of the fragment(s) formed after digesting seq
                with the restriction enzyme Acc65I in the order that they appear in seq.

                Use FASTA format and give the Watson strand(s) in 5'-3' direction below.
                Give the sequences the names frag1,frag2,... etc.
                >frag1
                CTCCCCTATCACCAGG

                >frag2
                GTACCGATAGCCACGAATCT

                *********** Question 4 ***********

                QuestionID:
                '''
        result = parse(text)

        correct = ['CTCCCCTATCACCAGGGTACCGATAGCCACGAATCT',
                   'CTCCCCTATCACCAGG',
                   'GTACCGATAGCCACGAATCT']

        self.assertEqual( [str(s.seq) for s in result], correct )

        self.assertEqual( [s.linear for s in result], [True,True,True] )


        input =   '''
                LOCUS       ScCYC1                   330 bp    DNA              UNK 01-JAN-1980
                DEFINITION  ScCYC1
                ACCESSION   ScCYC1
                VERSION     ScCYC1
                KEYWORDS    .
                SOURCE      .
                  ORGANISM  .
                            .
                FEATURES             Location/Qualifiers
                ORIGIN
                        1 ATGACTGAAT TCAAGGCCGG TTCTGCTAAG AAAGGTGCTA CACTTTTCAA GACTAGATGT
                       61 CTACAATGCC ACACCGTGGA AAAGGGTGGC CCACATAAGG TTGGTCCAAA CTTGCATGGT
                      121 ATCTTTGGCA GACACTCTGG TCAAGCTGAA GGGTATTCGT ACACAGATGC CAATATCAAG
                      181 AAAAACGTGT TGTGGGACGA AAATAACATG TCAGAGTACT TGACTAACCC AAAGAAATAT
                      241 ATTCCTGGTA CCAAGATGGC CTTTGGTGGG TTGAAGAAGG AAAAAGACAG AAACGACTTA
                      301 ATTACCTACT TGAAAAAAGC CTGTGAGTAA
                //
                '''
        result = parse(input).pop()

        self.assertEqual( str(result.seq) , str(read(input).seq) )

        correct = '''ATGACTGAATTCAAGGCCGGTTCTGCTAAGAAAGGTGCTACACTTTTCAAGACTAGATGTCTACAATGCCACACCGTGGAAAAGGGTGGCCCACATAAGGTTGGTCCAAACTTGCATGGTATCTTTGGCAGACACTCTGGTCAAGCTGAAGGGTATTCGTACACAGATGCCAATATCAAGAAAAACGTGTTGTGGGACGAAAATAACATGTCAGAGTACTTGACTAACCCAAAGAAATATATTCCTGGTACCAAGATGGCCTTTGGTGGGTTGAAGAAGGAAAAAGACAGAAACGACTTAATTACCTACTTGAAAAAAGCCTGTGAGTAA'''

        self.assertEqual( str(result.seq) , correct )

        self.assertTrue( result.linear   == True )
        self.assertTrue( result.circular == False )

        seqs = parse('./RefDataBjorn.fas')

        self.assertEqual( len(seqs) , 771 )
        self.assertEqual( list(set([len (a) for a in seqs])) ,[901])

        pAG25 = read("./pAG25.gb")

        self.assertTrue( pAG25.circular == True )
        self.assertTrue( pAG25.linear   == False)

        pCAPs = read("./pCAPs.gb")

        self.assertTrue( pCAPs.circular == True )
        self.assertTrue( pCAPs.linear   == False)

        pUC19 = read("./pUC19.gb")

        self.assertTrue( pUC19.circular == True )
        self.assertTrue( pUC19.linear   == False)



    def test_parse2(self):
        from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

        seqs = parse('./RefDataBjorn.fas')

        self.assertTrue( len(seqs) == 771 )
        self.assertTrue( list(set([len (a) for a in seqs])) == [901] )

        for i,s in enumerate(seqs):
            a = s.description
            b = a.split()
            c =  "|".join([b[0],b[1],b[3]])
            s.id = b[2].replace(" ","_")+"_"+str(i)
            s.description = ""
            if b[3]=="Zenion hololepis":
                s.id = b[3].replace(" ","_")+"_"+str(i)
            s.seq.alphabet = IUPACAmbiguousDNA()


if __name__ == '__main__':
    unittest.main()









