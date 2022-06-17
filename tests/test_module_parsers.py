#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test parse
"""

import pytest


def test_parse1():
    from pydna.parsers import parse
    from pydna.readers import read

    """ test parsing fasta sequences from a text"""

    text = """
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
            """
    result = parse(text)

    correct = [
        "CTCCCCTATCACCAGGGTACCGATAGCCACGAATCT",
        "CTCCCCTATCACCAGG",
        "GTACCGATAGCCACGAATCT",
    ]

    assert [str(s.seq) for s in result] == correct
    assert [s.linear for s in result] == [True, True, True]

    input = """
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
            """
    result = parse(input).pop()

    assert str(result.seq) == str(read(input).seq)
    correct = """ATGACTGAATTCAAGGCCGGTTCTGCTAAGAAAGGTGCTACACTTTTCAAGACTAGATGTCTACAATGCCACACCGTGGAAAAGGGTGGCCCACATAAGGTTGGTCCAAACTTGCATGGTATCTTTGGCAGACACTCTGGTCAAGCTGAAGGGTATTCGTACACAGATGCCAATATCAAGAAAAACGTGTTGTGGGACGAAAATAACATGTCAGAGTACTTGACTAACCCAAAGAAATATATTCCTGGTACCAAGATGGCCTTTGGTGGGTTGAAGAAGGAAAAAGACAGAAACGACTTAATTACCTACTTGAAAAAAGCCTGTGAGTAA"""

    assert str(result.seq) == correct
    assert result.linear == True
    assert result.circular == False

    seqs = parse("RefDataBjorn.fas")

    assert len(seqs) == 771
    assert list(set([len(a) for a in seqs])) == [901]
    pAG25 = read("pAG25.gb")

    assert pAG25.circular == True
    assert pAG25.linear == False

    pCAPs = read("pCAPs.gb")

    assert pCAPs.circular == True
    assert pCAPs.linear == False

    pUC19 = read("pUC19.gb")

    assert pUC19.circular == True
    assert pUC19.linear == False

    input = """
    ID   example    standard; DNA; UNC; 3 BP.
    SQ   Sequence 3 BP;
         aaa                                                                       3
    //
    """
    result = parse(input).pop()
    input = """
    ID   name?      standard; circular DNA; UNK; 100 BP.
    XX
    DT   25-DEC-2017
    XX
    DE   description?.
    XX
    AC   id?;
    XX
    SV   id?
    XX
    KW   .
    XX
    OS   .
    OC   .
    OC   .
    XX
    FH   Key             Location/Qualifiers
    SQ   Sequence 100 BP;
         aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa        60
         aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa                             100
    //
    """
    result = parse(input).pop()


def test_parse2():

    from pydna.parsers import parse
    from pydna.readers import read

    seqs = parse("RefDataBjorn.fas")

    assert len(seqs) == 771
    assert list(set([len(a) for a in seqs])) == [901]

    for i, s in enumerate(seqs):
        a = s.description
        b = a.split()
        c = "|".join([b[0], b[1], b[3]])
        s.id = b[2].replace(" ", "_") + "_" + str(i)
        s.description = ""
        if b[3] == "Zenion hololepis":
            s.id = b[3].replace(" ", "_") + "_" + str(i)


def test_parse_primers():
    from pydna.parsers import parse_primers

    data = str(">1\n" "aaaa\n" ">2\n" "cccc\n")
    parse_primers(data)


def test_parse_error():
    from pydna.parsers import parse

    data = """
LOCUS
DATA_IS_NOT_A_SEQUENCE
//"""
    parse(data)
    assert parse(data) == []


def test_parse_list():
    from pydna.parsers import parse_primers

    data = str(">1\n" "aaaa\n" ">2\n" "cccc\n")

    parse_primers([data, data])


def test_misc_parse():
    from pydna.parsers import parse

    from Bio.SeqIO import read as BPread
    from Bio.SeqIO import parse as BPparse

    q = BPread("read1.gb", "gb")
    w = BPread("read2.gb", "gb")
    e = BPread("read3.fasta", "fasta")
    r = BPread("read4.fasta", "fasta")

    with open("pth1.txt", "r", encoding="utf-8") as f:
        a, b = BPparse(f, "gb")

    assert "|" + a.features[13].qualifiers["label"][0] + "|" == "|2micron 2µ|"
    assert "|" + a.format("gb")[3314:3324] + "|" == '|olor="gree|'

    assert a.features[13].qualifiers["label"][0] == "2micron 2µ"
    assert a.format("gb")[3268:3278] == "2micron 2µ"

    x, y = parse("pth1.txt")

    assert "".join(a.format("gb").splitlines()[1:]) == "".join(
        x.format("gb").splitlines()[1:]
    )
    assert "".join(b.format("gb").strip().splitlines()[4:]) == "".join(
        y.format("gb").splitlines()[4:]
    )


def test_dna2949():
    from pydna.parsers import parse
    with open("dna2943.gb") as f:
        f.read()
    seqlist = parse("dna2943.gb", ds=True)
    assert len(seqlist) == 1
    assert seqlist[0].useguid() == "jkB2Ky9pW-hR7XCocz40PR_RKl4"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
