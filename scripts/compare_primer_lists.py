# coding: utf-8
from pydna import parse, eq

new_primer = parse("/home/bjorn/Dropbox/wikidata/PrimersToBuy.wiki", ds=False)
primer = parse("/home/bjorn/Dropbox/wikidata/Primers.wiki", ds=False)

primer = primer[::-1]
old_primer = primer[: 37 + 18]
primer = primer[37 + 18 :]

new_primer_dict = dict((p.id, p) for p in new_primer)
primer_dict = dict((p.id, p) for p in primer)
old_primer_dict = dict((p.id, p) for p in old_primer)

assert str(primer_dict["509_mycGFPr"].seq) == "CTACTTGTACAGCTCGTCCA"
assert primer[0].id == "0_S1"
assert primer[1].id == "1_5CYC1clone"
assert primer[580].id == "580_GXF1_YPK_fwd"

fprimer = parse("/home/bjorn/Dropbox/wikidata/New10.wiki", ds=False)
fprimer = fprimer[::-1]
fold_primer = fprimer[: 37 + 18]
fprimer = fprimer[37 + 18 :]

print len(primer), len(fprimer)

for i in range(1, len(primer)):
    if not eq(primer[i], fprimer[i]):
        print i
        print primer[i].id, fprimer[i].id
        raw_input()
