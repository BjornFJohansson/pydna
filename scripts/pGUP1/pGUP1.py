#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydna import Dseqrecord, read, pcr, circular_assembly, sync
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid


info = """
pGUP1
=====
This script carries out the cloning project described in the second
example in the cookbook.

This file shows an example construction of the vector pGUP1, described
in reference (1) below.

Quote from the publication (1):

"The expression vectors harboring GUP1 or GUP1H447A were obtained as
follows: the open reading frame of GUP1 was amplified by PCR using plasmid
pBH2178 (kind gift from Morten Kielland-Brandt) as a template and using
primers GUP1rec1sens (5'-gaattcgatatcaagcttatcgataccgatgtcgctgatcagcatcctgtc-
tcc-3') and GUP1rec2AS (5'-gacataactaattacatgactcgaggtcgactcagcattttaggtaaatt-
ccg-3'), underlined sequences being homologous to the target vector pGREG505
(Jansen et al., 2005). The PCR fragment was purified by a PCR purification kit
(QIAGEN, Chatsworth, CA) and introduced into pGREG505 by cotransfection
into yeast cells thus generating pGUP1 (Jansen et al., 2005)."

This is a cloning in three steps:

A. PCR of the GUP1 locus using GUP1rec1sens GUP1rec2AS, resulting in
   a linear insert.

B. Digestion of the plasmid pGREG505 with SalI, This step is not
   mentioned above, but evident from (2). This digestion removes a
   DNA fragment containing the HIS3 marker gene from the final
   construct.

C. Recombination between the linear insert and the linear vector.

This cloning procedure is replicated using Python-dna. Execution of
this file will write the sequence of the pGUP1 plasmid in Genbank
format to a file called "pGUP1.gb". The sequence will be in Genbank
format.

References
----------

1) Régine Bosson, Malika Jaquenoud, and Andreas Conzelmann,
“GUP1 of Saccharomyces Cerevisiae Encodes an O-acyltransferase
Involved in Remodeling of the GPI Anchor,” Molecular Biology of
the Cell 17, no. 6 (June 2006): 2636–2645.
http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1474799/

2) Jansen G, Wu C, Schade B, Thomas DY, Whiteway M. 2005. Drag&Drop
cloning in yeast. Gene, 344: 43–51.
http://www.ncbi.nlm.nih.gov/pubmed/15656971
"""

print info

raw_input("Press any key and wait for the script to finish!")

# Establish the two primers. These sequences can be found in (1)
GUP1rec1sens = SeqRecord(Seq("gaattcgatatcaagcttatcgataccgatgtcgctgatcagcatcctgtctcc"))
GUP1rec2AS = SeqRecord(Seq("gacataactaattacatgactcgaggtcgactcagcattttaggtaaattccg"))

# Read the GUP1 locus sequence into a Dseqrecord object
# This sequence was taken from the Saccharomyces genome Database:
# http://www.yeastgenome.org/cgi-bin/getSeq?query=YGL084C&flankl=1000&flankr=1000&format=fasta
GUP1 = read("GUP1_locus.gb")

# The insert is formed by PCR using the two primers and the template sequence
insert = pcr(GUP1rec1sens, GUP1rec2AS, GUP1)

# The sequence for the plasmid is read into a Dseqrecord object called pGREG505
# this sequence was found at
# http://www.euroscarf.de/plasmid_details.php?accno=P30350
# This sequence is circular, this information is parsed from the Genbank file.
pGREG505 = read("pGREG505.gb")

# Import the SalI restriction enzyme from Biopython
from Bio.Restriction import SalI

# Cut the circular pGREG505 plasmid with SalI
# this enzyme cuts twice, so two fragments are formed
linear_vector, his3 = pGREG505.cut(SalI)

# Circular recombination products are formed
# limit is the length of the necessary regions of
# homology between the seqences
formatted_sequences, circular_recombination_products = circular_assembly(
    (insert, linear_vector), limit=28
)

# The circular recombination products are returned
# by order of size. In this case, the largest one is the
# correct one
pGUP1 = circular_recombination_products[0]

# The circular recombination products are returned
# by order of size. In this case, the largest one is the
# correct one
pGUP1 = sync(pGUP1, pGREG505)

pGUP1.write("pGUP1.gb")

print "done! look at the file pGUP1.gb"
