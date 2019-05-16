#!/home/bjorn/anaconda3/envs/bjorn36/bin/python
# -*- coding: utf8 -*-
"""
@author: bjorn
"""

from pydna.myprimers import list_primers

for i, p in enumerate(list_primers):
    print(i, end=" ")
    assert p.name.startswith(str(i))

if __name__ is "__main__":
    pass
    input("press enter to close")
