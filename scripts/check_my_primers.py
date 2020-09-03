#!/home/bjorn/anaconda3/envs/bjorn36/bin/python
# -*- coding: utf8 -*-
"""
@author: bjorn
"""

from pydna.myprimers import list_primers

print("number of primers", len(list_primers))
print("unique primers", len(set(str(p.seq).lower() for p in list_primers)))
print("no seq (n)", len([p for p in list_primers if set(p.seq.lower()) == set("n")]))

for i, p in enumerate(list_primers):
    if not p.name.startswith(str(i)):
        print(i, p.format("tab"))

print("names checked for primer number for ", i + 1, "primers")

duplicates = []
for i1, p1 in enumerate(list_primers):
    for i2, p2 in enumerate(list_primers):
        if set(p1.seq.lower()) == set("n") or set(p2.seq.lower()) == set("n"):
            continue
        if i1 != i2:
            if str(p1.seq).lower() == str(p2.seq).lower():
                duplicates.extend((p1, p2))


duplicates = list(dict.fromkeys(duplicates))

for d in duplicates:
    print(d.format("fasta"))


if __name__ == "__main__":
    input("press enter to close")
