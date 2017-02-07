#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test parse
'''

import pytest
import sys

from pydna.parsers import parse
from pydna.readers import read

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

    assert [str(s.seq) for s in result] == correct
    assert [s.linear for s in result] == [True,True,True]


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

    assert str(result.seq) == str(read(input).seq)
    correct = '''ATGACTGAATTCAAGGCCGGTTCTGCTAAGAAAGGTGCTACACTTTTCAAGACTAGATGTCTACAATGCCACACCGTGGAAAAGGGTGGCCCACATAAGGTTGGTCCAAACTTGCATGGTATCTTTGGCAGACACTCTGGTCAAGCTGAAGGGTATTCGTACACAGATGCCAATATCAAGAAAAACGTGTTGTGGGACGAAAATAACATGTCAGAGTACTTGACTAACCCAAAGAAATATATTCCTGGTACCAAGATGGCCTTTGGTGGGTTGAAGAAGGAAAAAGACAGAAACGACTTAATTACCTACTTGAAAAAAGCCTGTGAGTAA'''

    assert str(result.seq) == correct
    assert result.linear   == True 
    assert result.circular == False 

    seqs = parse('RefDataBjorn.fas')

    assert len(seqs) == 771    
    assert list(set([len (a) for a in seqs])) == [901]
    pAG25 = read("pAG25.gb")

    assert pAG25.circular == True 
    assert pAG25.linear   == False

    pCAPs = read("pCAPs.gb")

    assert pCAPs.circular == True 
    assert pCAPs.linear   == False

    pUC19 = read("pUC19.gb")

    assert pUC19.circular == True 
    assert pUC19.linear   == False



def test_parse2():
    from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

    seqs = parse('RefDataBjorn.fas')

    assert len(seqs) == 771
    assert list(set([len (a) for a in seqs])) == [901]

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
    pytest.cmdline.main([__file__, "-v", "-s"])








