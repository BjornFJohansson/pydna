from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import EcoRI

# Create a Dseqrecord with your DNA sequence
sequence = "GAATTCGCGGATCCGTCGACGAATTC"
record = Dseqrecord(sequence)

# Cut with a single enzyme
cut_records = record.cut(EcoRI)

# Cut with multiple enzymes
cut_records_multiple = record.cut([EcoRI, BamHI])

# Display the resulting fragments
for frag in cut_records:
    print(frag)

for frag in cut_records_multiple:
    print(frag)