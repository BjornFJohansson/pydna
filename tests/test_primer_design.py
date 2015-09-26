#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test parse
'''

import nose, sys, os

from pydna import pcr, Assembly, Dseqrecord, assembly_primers, cloning_primers, parse

from textwrap import dedent

def test_primer_Design():
    ''' test_primer_design'''

    a=Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")
    b=Dseqrecord("ccaaacccaccaggtaccttatgtaagtacttcaagtcgccagaagacttcttggtcaagttgcc")
    c=Dseqrecord("tgtactggtgctgaaccttgtatcaagttgggtgttgacgccattgccccaggtggtcgtttcgtt")

    primer_pairs = assembly_primers([a,b,c])

    frags=[]

    for (f,r),t in zip(primer_pairs,[a,b,c]):
        frags.append(pcr(f,r,t))

    asm=Assembly(frags)

    assert asm.linear_products[0].seguid() == "1eNv3d_1PqDPP8qJZIVoA45Www8"

    frags=[]

    primer_pairs = assembly_primers([a,b,c], circular=True)

    for (f,r),t in zip(primer_pairs,[a,b,c]):
        frags.append(pcr(f,r,t))

    #print frags

    asm=Assembly(frags)

    assert asm.circular_products[0].cseguid() == "V3Mi8zilejgyoH833UbjJOtDMbc"



def test_primer_Design_with_linker():
    ''' test_primer_design'''


    from pydna import Dseqrecord, Assembly, pcr, assembly_primers

    b  = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    l  = Dseqrecord("AAATTTCCCGGG")
    c  = Dseqrecord("tctgatctactatctgactgtactgatctattgacactgtgatcattctagtgtattactc")

    ((bf,br),(cf,cr)) = assembly_primers((b,l,c))

    nb = pcr((bf,br),b)
    nc = pcr((cf,cr),c)

    asm1 = Assembly((nb,nc))

    assert asm1.linear_products[0].seguid(),(b+l+c).seguid() == 'l95igKB8iKAKrvvqE9CYksyNx40'


    b  = Dseqrecord("agctactgactattaggggttattctgatcatctgatctactatctgactgtactgatcta")
    l  = Dseqrecord("AAATTTCCCGGG")
    c  = Dseqrecord("tctgatctactatctgactgtactgatctattgacactgtgatcattctagtgtattactc")

    ((bf,br),(cf,cr)) = assembly_primers((b,l,c), circular = True)

    nb = pcr((bf,br),b)
    nc = pcr((cf,cr),c)

    asm = Assembly((nb,nc))

    #print (b+l+c).looped().seq

    assert (b+l+c).looped().cseguid() == asm.circular_products[0].cseguid()
    #print (b+l+c).looped().cseguid() == 'jdHXfQI5k4Sk2ESiZYfKv4oP2FI'

    assert (b+l+c).looped().cseguid() == 'jdHXfQI5k4Sk2ESiZYfKv4oP2FI'


def test_primer_Design_saving_to_text_file():

    files_ = [dedent(x) for x in
    ('''\
    >fw64 t-sequence
    atgactgctaacccttc
    >rv64 t-sequence
    catcgtaagtttcgaac''',

    '''\
    >fw64 t-sequence
    atgactgctaacccttcN
    >rv64 t-sequence
    catcgtaagtttcgaacN''',

    '''\
    >fw64 t-sequence
    atgactgctaacccttcN''',

    '''\
    >fw64 t-sequence
    atgactgctaacccttcN
    >rv64 t-sequence
    catcgtaagtttcgaac''',

    '''\
    >rv64 t-sequence
    catcgtaagtttcgaacN''',

    '''\
    >rv64 t-sequence
    catcgtaagtttcgaacN
    >fw64 t-sequence
    atgactgctaacccttc''',

    '''\
    >fw64 t-sequenceZ
    atgactgctaacccttc
    >rv64 t-sequenceZ
    catcgtaagtttcgaac''',

    '''\
    >fw64 t-sequenceZ
    atgactgctaacccttc
    >rv64 t-sequenceZ
    catcgtaagtttcgaac
    >fw64 t-sequence
    atgactgctaacccttc
    >rv64 t-sequence
    catcgtaagtttcgaac''',

    '''\
    >fw64x t-sequence
    atgactgctaacccttc
    >rv64x t-sequence
    catcgtaagtttcgaac''',

    '''\
    >fw64x t-sequence
    atgactgctaacccttc
    >rv64x t-sequence
    catcgtaagtttcgaac
    >fw64 t-sequence
    atgactgctaacccttc
    >rv64 t-sequence
    catcgtaagtttcgaac''',
    )]

    templ=Dseqrecord("atgactgctaacccttccttggtgttgaacaagatcgacgacatttcgttcgaaacttacgatg")

    templ.accession = "t-sequence"

    try:
        os.remove("PRIMERS.TXT")
    except OSError:
        pass

    pf, pr = cloning_primers(templ, path="PRIMERS.TXT")

    assert os.path.isfile("PRIMERS.TXT")

    with open("PRIMERS.TXT", "r") as f: text=f.read().strip()

    assert text == files_[0]

    with open("PRIMERS.TXT", "w") as f: text=f.write(files_[1])

    pf, pr = cloning_primers(templ, path="PRIMERS.TXT")

    with open("PRIMERS.TXT", "r") as f: text=f.read().strip()

    assert text == files_[1]

    with open("PRIMERS.TXT", "w") as f: text=f.write(files_[2])

    pf, pr = cloning_primers(templ, path="PRIMERS.TXT")

    with open("PRIMERS.TXT", "r") as f: text=f.read().strip()

    assert text == files_[3]

    with open("PRIMERS.TXT", "w") as f: text=f.write(files_[4])

    pf, pr = cloning_primers(templ, path="PRIMERS.TXT")

    with open("PRIMERS.TXT", "r") as f: text=f.read().strip()

    assert text == files_[5]

    with open("PRIMERS.TXT", "w") as f: text=f.write(files_[6])

    pf, pr = cloning_primers(templ, path="PRIMERS.TXT")

    with open("PRIMERS.TXT", "r") as f: text=f.read().strip()

    assert text == files_[7]

    with open("PRIMERS.TXT", "w") as f: text=f.write(files_[8])

    pf, pr = cloning_primers(templ, path="PRIMERS.TXT")

    with open("PRIMERS.TXT", "r") as f: text=f.read().strip()

    assert text == files_[9]

    with open("PRIMERS.TXT", "w") as f: text=f.write('')

    pf, pr = cloning_primers(templ, path="PRIMERS.TXT")

    f, r = parse("PRIMERS.TXT", ds=False)

    assert f.description == "fw64 t-sequence"
    assert r.description == "rv64 t-sequence"

if __name__ == '__main__':
    nose.runmodule(argv=[sys.argv[0], '--nocapture'])









