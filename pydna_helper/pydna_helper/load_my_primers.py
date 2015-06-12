#!/usr/bin/env python
# -*- coding: utf-8 -*-

global primer
global old_primer
global new_primer_dict
global primer_dict
global old_primer_dict

from settings import page_with_primers
from pydna import parse

primer     = parse(page_with_primers, ds=False)

primer     = primer[::-1]
old_primer = primer[:37+18]
primer     = primer[37+18:]

primer_dict             = dict((p.id, p) for p in primer)
old_primer_dict         = dict((p.id, p) for p in old_primer)

for i, p in enumerate(primer):
    globals()["p{:03d}".format(i)] = p

if __name__=="__main__":
    print "primers loaded into memory"
    print "{:3d} old primers    -> old_primer [list]".format(len(old_primer))
    print "{:3d} primers        -> primer [list]".format(len(primer))
    print "first primer ", primer[1].id
    print "last ({}) primer is {}".format(len(primer)-1, primer[-1].id)


    assert str(primer_dict["509_mycGFPr"].seq) == "CTACTTGTACAGCTCGTCCA"
    assert primer[0].id == "0_S1"
    assert primer[580].id == "580_GXF1_YPK_fwd"
    assert primer[803].id == "803_pYPK0_del_r"

    assert primer[100].id == "100_Jen1dom1rev"
    assert primer[200].id == "200_SRO7_1756_rev"
    assert primer[300].id == "300_C-STE2"
    assert primer[400].id == "400_RPS10_rv"
    assert primer[500].id == "500_pCAPsEcoRV5Pr"
    assert primer[600].id == "600_Sc_ade2_rev"
    assert primer[700].id == "700_sc_fas2-B1:"
    assert primer[800].id == "800_ELO1_rv"
    assert primer[900].id == "900_AIDxrv"
