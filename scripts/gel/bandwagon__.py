#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from bandwagon import BandsPattern
from bandwagon import BandsPatternsSet
from bandwagon import custom_ladder

LADDER_100_to_4k = custom_ladder(
    "100-4k",
    {
        100: 205,
        200: 186,
        300: 171,
        400: 158,
        500: 149,
        650: 139,
        850: 128,
        1000: 121,
        1650: 100,
        2000: 90,
        3000: 73,
        4000: 65,
    },
)

ladder = LADDER_100_to_4k.modified(label="Ladder", background_color="#ffffaf")


patterns = [
    BandsPattern(
        [10000, 7750, 4000, 3000, 2000, 1500, 1000, 750, 500], ladder, label="C1"
    ),
    BandsPattern([100, 500, 3500], ladder, label="C1"),
    BandsPattern([300, 400, 1500], ladder, label="C2"),
    BandsPattern([100, 1200, 1400, 3000], ladder, label="C3"),
    BandsPattern([100, 700], ladder, label="C4"),
]
patterns_set = BandsPatternsSet(
    patterns=[ladder] + patterns, ladder=ladder, label="Test pattern", ladder_ticks=3
)
ax = patterns_set.plot()
# ax.figure.savefig("simple_band_patterns.png", bbox_inches="tight", dpi=200)


from scipy.interpolate import CubicSpline

ipolatr = CubicSpline(
    [100, 200, 300, 400, 500, 650, 850, 1000, 1650, 2000, 3000, 4000],
    [205, 186, 171, 158, 149, 139, 128, 121, 100, 90, 73, 65],
    bc_type="natural",
    extrapolate=True,
)

for b in [10000, 7750, 4000, 3000, 2000, 1500, 1000, 750, 500]:
    print(ipolatr(b))


from math import log

x = 100
fx = 395 - 39.8 * log(x)
print(fx)
