#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydna.all import *

MySeq2 = read("/home/bjorn/Desktop/pydna_test/MySeq2_copied.gb")
with open("/home/bjorn/Desktop/pydna_test/MySeq2_copied.gb") as f:
    MySeq2_copied = read(f.read())

# reading looses text uppercase/lowercase

assert MySeq2_copied.format() == MySeq2.format()

# print(MySeq2.list_features())
# +-----+--------------+-----+-----+-----+-----+-----------+------+
# | Ft# | Label/Note   | Dir | Sta | End | Len | type      | orf? |
# +-----+--------------+-----+-----+-----+-----+-----------+------+
# |   0 | L:MyGene     | --> | 7   | 73  |  66 | CDS       | yes  |
# |   1 | L:MyMutation | --> | 31  | 32  |   1 | variation |  no  |
# +-----+--------------+-----+-----+-----+-----+-----------+------+

mg = MySeq2.extract_feature(0)
mg2 = MySeq2.features[0].extract(MySeq2)

assert mg.format() == mg2.format(), mg == mg2

mgc = MySeq2_copied.extract_feature(0)
mgc2 = MySeq2_copied.features[0].extract(MySeq2_copied)

assert mgc.format() == mgc2.format(), mgc == mgc2

assert mg == mg2 != mgc == mgc2  # a path property is added to mg and mg2 from MySeq2

assert mg2.format() == mgc.format()

mp = primer_design(mg)

## default Dseqrecord field content should be globally available somehow


slc = MySeq2[2:-2]

assert slc.name == slc.id == "MySeq2"


MySeq = Dseqrecord(
    "atgtcgtATGaaaccgttatcgatcatatgtGcgaaatgtcgcgcgtcatctacgtatcatcgatctactTAAacgtgta"
)

# LOCUS       name                      80 bp    DNA     linear   UNK 18-APR-2018
# DEFINITION  description.
# ACCESSION   id
# VERSION     id
# KEYWORDS    .
# SOURCE      .
#  ORGANISM  .
#            .
# FEATURES             Location/Qualifiers
# ORIGIN
#        1 atgtcgtatg aaaccgttat cgatcatatg tgcgaaatgt cgcgcgtcat ctacgtatca
#       61 tcgatctact taaacgtgta
# //
