#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test parse
'''

import nose, sys

from nose.tools import assert_equal, assert_true

from pydna import parse, read


def test_parse1():
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

    assert_equal( [str(s.seq) for s in result], correct )

    assert_equal( [s.linear for s in result], [True,True,True] )


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

    assert_equal( str(result.seq) , str(read(input).seq) )

    correct = '''ATGACTGAATTCAAGGCCGGTTCTGCTAAGAAAGGTGCTACACTTTTCAAGACTAGATGTCTACAATGCCACACCGTGGAAAAGGGTGGCCCACATAAGGTTGGTCCAAACTTGCATGGTATCTTTGGCAGACACTCTGGTCAAGCTGAAGGGTATTCGTACACAGATGCCAATATCAAGAAAAACGTGTTGTGGGACGAAAATAACATGTCAGAGTACTTGACTAACCCAAAGAAATATATTCCTGGTACCAAGATGGCCTTTGGTGGGTTGAAGAAGGAAAAAGACAGAAACGACTTAATTACCTACTTGAAAAAAGCCTGTGAGTAA'''

    assert_equal( str(result.seq) , correct )

    assert_true( result.linear   == True )
    assert_true( result.circular == False )

    seqs = parse('tests/RefDataBjorn.fas')

    assert_equal( len(seqs) , 771 )
    assert_equal( list(set([len (a) for a in seqs])) ,[901])

    pAG25 = read("tests/pAG25.gb")

    assert_true( pAG25.circular == True )
    assert_true( pAG25.linear   == False)

    pCAPs = read("tests/pCAPs.gb")

    assert_true( pCAPs.circular == True )
    assert_true( pCAPs.linear   == False)

    pUC19 = read("tests/pUC19.gb")

    assert_true( pUC19.circular == True )
    assert_true( pUC19.linear   == False)



def test_parse2():
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

    seqs = parse('tests/RefDataBjorn.fas')

    assert_true( len(seqs) == 771 )
    assert_true( list(set([len (a) for a in seqs])) == [901] )

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
    nose.runmodule(argv=[sys.argv[0], '--nocapture'])








