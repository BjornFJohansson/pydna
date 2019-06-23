#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 09:46:13 2019

@author: bjorn
"""


from bandwagon import BandsPattern, BandsPatternsSet, LADDER_100_to_4k

ladder = LADDER_100_to_4k.modified(label="Ladder", background_color="#ffffaf")

patterns = [
    BandsPattern([100, 500, 3500], ladder, label="C1"),
    BandsPattern([300, 400, 1500], ladder, label="C2"),
    BandsPattern([100, 1200, 1400, 3000], ladder, label="C3"),
    BandsPattern([100, 700], ladder, label="C4"),
]
patterns_set = BandsPatternsSet(patterns=[ladder] + patterns, ladder=ladder,
                                label="Test pattern", ladder_ticks=3)
ax = patterns_set.plot()
#ax.figure.savefig("simple_band_patterns.png", bbox_inches="tight", dpi=200)