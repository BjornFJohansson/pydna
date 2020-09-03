#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Sat Mar 30 06:59:58 2019 @author: bjorn

import colorsys

from pydna.dseqrecord import Dseqrecord

from Bio.SeqFeature import SeqFeature, FeatureLocation

seq = Dseqrecord("")

colors = []

for hue in range(0, 360, 36):
    for s in [0.2, 0.4, 0.6]:
        r, g, b = colorsys.hsv_to_rgb(hue / 360, s, 255)
        r, g, b = int(r), int(g), int(b)
        hex = "#{0:02x}{1:02x}{2:02x}".format(r, g, b)
        colors.append(f"{hex}")

colors = colors[::3] + colors[1::3] + colors[2::3]
colors = colors[::-1]
for hex in colors:
    sf = SeqFeature(FeatureLocation(1, 9, strand=1), type="misc_feature")
    sf.qualifiers["label"] = [
        hex,
    ]
    sf.qualifiers["ApEinfo_fwdcolor"] = [
        hex,
    ]
    se = Dseqrecord("agtagtcgta")
    se.features.append(sf)
    seq += se

from pydna.editor import ape

ape(seq)

print(colors)
