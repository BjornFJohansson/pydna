#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat Mar 30 06:59:58 2019 @author: bjorn
      
import colorsys

from pydna.dseqrecord import Dseqrecord

from Bio.SeqFeature import SeqFeature, FeatureLocation

seq = Dseqrecord("")

newcodes =[]

for hue in range(0,360,36):
    for s in [0.2, 0.4, 0.6]:
        r,g,b = colorsys.hsv_to_rgb(hue/360, s, 255)
        r,g,b = int(r),int(g),int(b)
        hex="#{0:02x}{1:02x}{2:02x}".format(r,g,b)
        sf = SeqFeature(FeatureLocation(1,9, strand=1), type='misc_feature')
        sf.qualifiers["label"]=[hex,]
        sf.qualifiers["ApEinfo_fwdcolor"]=[hex,]
        se = Dseqrecord("agtagtcgta")
        se.features.append(sf)
        seq+=se
        
        print(hex)

from pydna.editor import ape

ape(seq)

#ffcccc
#ff9999
#ff6666
#ffeacc
o#ffd699
#ffc166
#f4ffcc
#eaff99
#e0ff66
#d6ffcc
#adff99
#84ff66
#ccffe0
#99ffc1
#66ffa3
#ccffff
#99ffff
#66ffff
#cce0ff
#99c1ff
#66a3ff
#d6ccff
#ad99ff
#8466ff
#f4ccff
#ea99ff
#e066ff
#ffccea
#ff99d6
#ff66c1
    
    
    