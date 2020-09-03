#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from lxml import etree

doc = etree.parse("labhost_Bacillus_subtilis.xml")

root = doc.getroot()

for comp in root.iter("{http://sbols.org/v1#}component"):
    for c in next(comp.iterchildren()):
        print(c.tag)
        print(c)

        break
