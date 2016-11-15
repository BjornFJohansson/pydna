import pydna

from Bio.SeqIO import read as BPread
from Bio.SeqIO import parse as BPparse

q = BPread("read1.gb", "gb")
w = BPread("read2.gb", "gb")
e = BPread("read3.fasta", "fasta")
r = BPread("read4.fasta", "fasta")


with open("pth1.txt", "r", encoding="utf-8") as f:
    a, b = BPparse(f, "gb")

print("|"+a.features[13].qualifiers['label'][0]+"|")
print("|"+a.format("gb")[3314:3324]+"|")

assert a.features[13].qualifiers['label'][0] == '2micron 2µ'
assert a.format("gb")[3314:3324] == '2micron 2µ'

x, y = pydna.parse("pth1.txt")

assert "".join(a.format("gb").splitlines()[1:]) == "".join(x.format("gb").splitlines()[1:])
assert "".join(b.format("gb").strip().splitlines()[4:]) == "".join(y.format("gb").splitlines()[4:])
