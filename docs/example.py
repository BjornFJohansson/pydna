from pydna.dseqrecord import Dseqrecord
dsr = Dseqrecord("ccccGGATCCatgccctaaGGATCCaaaa")
dsr.add_feature(10,19) # a small gene
dsr.figure()
# In[]
from Bio.Restriction import BamHI
a,b,c  = dsr.cut(BamHI)

print(a.figure())

print(b.figure())

print(c.figure())
# In[]
vector = Dseqrecord("aatgtttttcgctgacaatcataATAGATCTtgctatgcatcatcgatct", circular=True)
from Bio.Restriction import BglII

linear_vector = vector.linearize(BglII)

rec_vector = (linear_vector + b).looped().synced(vector)

rec_vector.figure()
