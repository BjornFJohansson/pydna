from pydna.dseqrecord import Dseqrecord
from pydna          import assembly
from pydna.parsers  import parse
from pydna.utils    import eq
    
a = Dseqrecord("ACTACGGCCTTCTCTCCCCCtgtgctgtgctcta",name="one34")
a.add_feature(1,33,label="first")  

b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatct",name="two35")
b.add_feature(1,34,label="scnd")                     

c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg",name="three37")
c.add_feature(1,36,label="third")

ln0 = assembly.Assembly((a,b,c), limit=14)
l = ln0.assemble_linear()[0]
assert str(l.seq)=='ACTACGGCCTTCTCTCCCCCtgtgctgtgctctaTTTTTtattctggctgtatctGGGGGTacgatgctatactgg'

a = Dseqrecord("acgatgctatactggCCCCCtgtgctgtgctct", name="one")
b = Dseqrecord("tgtgctgtgctctTTTTTtattctggctgtat", name="two")
c = Dseqrecord("tattctggctgtatGGGGGtacgatgctatactgg", name="three")

c0 = assembly.Assembly((a,b,c), limit=13)

z=c0.assemble_circular()[0]

y=z.rc()