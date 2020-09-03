# coding: utf-8

from pydna import read, parse, pcr, sync, Genbank
from Bio.Restriction import BglII, BamHI

print """
YEp24PGK
========

This script carries out the cloning project described in the first
example in the cookbook.

This script simulates the construction of the vector YEp24PGK_XK by
PCR amplification of the Saccharomyces cerevisiae XKS1 (Genbank
Z72979) gene from chromosomal DNA using primer1 and primer3 and the
subsequent sub cloning of the PCR fragment in vector YEP24PGK
(Genbank KC562906) by ligation of vector digested with BglII and the
PCR product digested with BamHI as described in the Materials and
Methods section on page 4250 in:

Johansson B et al. (2001), Xylulokinase Overexpression in Two Strains
of Saccharomyces Cerevisiae Also Expressing Xylose Reductase and
Xylitol Dehydrogenase and Its Effect on Fermentation of Xylose and
Lignocellulosic Hydrolysate, Applied and Environmental Microbiology
67 4249–4255.

"""

raw_input("press return!\n")


gb = Genbank("me@home.org")

if gb.test():
    xks1_gene = gb.nucleotide("Z72979")
    print "Genbank record Z72979 downloaded from NCBI"
    YEp24PGK = gb.nucleotide("KC562906")
    print "Genbank record KC562906 downloaded from NCBI\n"
else:
    xks1_gene = read("Z72979.gb")
    print "A local copy of Genbank record Z72979 is used"
    YEp24PGK = read("KC562906.gb")
    print "A local copy of Genbank record KC562906 is used\n"

raw_input("press return!\n")

primers = """
>primer1
GCGGATCCTCTAGAATGGTTTGTTCAGTAATTCAG
>primer3
AGATCTGGATCCTTAGATGAGAGTCTTTTCCAG
"""
primer1, primer2 = parse(primers, ds=False)
xks1_pcr_product = pcr(primer1, primer2, xks1_gene)

YEp24PGK_bgl = YEp24PGK.cut(BglII).pop()
stuffer1, xks1_bam, stuffer2 = xks1_pcr_product.cut(BamHI)

YEp24PGK_XK = (YEp24PGK_bgl + xks1_bam.rc()).looped()

YEp24PGK_XK = YEp24PGK_XK.synced(YEp24PGK)

print "The sequence of YEp24PGK_XK was generated"

print "Seguid of YEp24PGK_XK is correct", YEp24PGK_XK.seguid() == "HRVpCEKWcFsKhw/W+25ednUfldI"

YEp24PGK_XK.write("YEp24PGK_XK.gb")

print """
done! The file YEp24PGK_XK.gb was written to the current
working directory
"""
