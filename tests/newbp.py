from Bio.Seq import Seq

from Bio.SeqRecord import SeqRecord

sr = SeqRecord(Seq("GATC"))

from Bio.SeqUtils import CheckSum
'dYoa6Z6T6WtPxu5msRTifUobfm0'
CheckSum.seguid(sr.seq)

sr.__str__()
sr.seq.__str__()
sr.seq.__repr__()

from pydna.dseqrecord import Dseqrecord

s = Dseqrecord("GATC")






from pydna.dseq import Dseq

fdseq = Dseq("GAT", "ATC")
sdseq = Dseq("GATC", "GATC")

nw = Dseq.quick("GATC", "GATC")

from pydna.dseqrecord import Dseqrecord

fdseqr = Dseqrecord(fdseq)
rdseqr = Dseqrecord(sdseq)

fdseqr+rdseqr

fdseqr.seguid()
from Bio.Restriction import BamHI
fdseq = Dseq.from_string("GGATCC")
fdseq = Dseq("GGATCC")
fdseq = Dseq("GGATCC", linear=False)
fdseq.cut(BamHI)



from pydna.readers import read
from pydna.amplify import Anneal
from pydna.dseqrecord import Dseqrecord

t = Dseqrecord("tacactcaccgtctatcattatcta" +
               "gatc"*240 +
               "ctatcgactgtatcatctgatagcac")

t.__str__()




from pydna.readers import read
from Bio.SeqRecord import SeqRecord
p1 = read(">p1\ntacactcaccgtctatcattatc", ds = False)
p2 = read(">p2\ngtgctatcagatgatacagtcg", ds = False)

p1.__str__()

ann = Anneal((p1, p2), t)
print(ann.report())
ann.products
amplicon_list = ann.products
amplicon = amplicon_list.pop()
print(amplicon.figure())

amplicon.__str__()

print(amplicon)




from pydna.dseq  import Dseq 
s=Dseq("atgtacgatcgtatgctggttatattttag")
s.translate()


