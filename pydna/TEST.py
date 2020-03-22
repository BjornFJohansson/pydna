#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sun Feb  9 07:31:39 2020 @author: bjorn

from pydna.dseqrecord import Dseqrecord

from pydna.amplify import pcr
from pydna.primer import Primer
template  = Dseqrecord("ATGGCAGTTGAGAAGActaTCTTCTCAACTGCCAT")

p1 = Primer("ATGGCAGTTGAGAAGA")
#p2 = Primer("cgactgtatcatctgatagcac").reverse_complement()

x = pcr(p1, template)

print(x.forward_primer.tm())

print(x.forward_primer.concentration)
print(x.concentration)
print(x.program())