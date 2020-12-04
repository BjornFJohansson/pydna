#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# doctest: +NORMALIZE_WHITESPACE
# doctest: +SKIP
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""docstring."""

from PIL import Image, ImageDraw
import numpy as np
import math
from scipy.interpolate import CubicSpline
from pydna.ladders import PennStateLadder


interpolator = CubicSpline(
    [int(bp) for bp in ("500 750 1000 1500 2000"
                        " 3000 4000 5000 7750 10000").split()],
    [int(px) for px in "366 296 246 183 146 104 84 70 50 41".split()],
    bc_type="natural",
    extrapolate=False,)


def gel(samples=[PennStateLadder, ],
        gel_length=600,
        margin=50,
        interpolator=interpolator):
    """docstring."""
    max_intensity = 256
    lane_width = 50
    lanesep = 10
    start = 10
    width = int(60 + (lane_width + lanesep) * len(samples))
    lanes = np.zeros((len(samples), gel_length), dtype=np.int)
    image = Image.new("RGB", (width, gel_length), "#ddd")
    draw = ImageDraw.Draw(image)
    draw.rectangle((0, 0, (width, gel_length)), fill=(0, 0, 0))
    scale = (gel_length-margin)/interpolator(min(interpolator.x))

    for labelsource in samples[0]:
        peak_centre = (interpolator(len(labelsource))) * scale - 5 + start
        label = f"{len(labelsource):<5} -"
        draw.text((2, peak_centre), label, fill=(255, 255, 255))

    for lane_number, lane in enumerate(samples):
        for band in lane:
            log = math.log(len(band), 10)
            height = (band.m() / (240 * log)) * 1e10
            peak_centre = interpolator(len(band)) * scale + start
            max_spread = 10
            if len(band) < 50:
                peak_centre += 50
                max_spread *= 4
                max_intensity /= 10
            band_spread = max_spread / log
            for i in range(max_spread, 0, -1):
                y1 = peak_centre - i
                y2 = peak_centre + i
                intensity = height * math.exp(-float(
                    ((y1 - peak_centre) ** 2)) / (2 * (band_spread**2))
                    ) * max_intensity
                for y in range(int(y1), int(y2)):
                    lanes[lane_number][y] += intensity

    for i, lane in enumerate(lanes):
        max_intensity = np.amax(lanes[i])
        if max_intensity > 256:
            lanes[i] = np.multiply(lanes[i], 256)
            lanes[i] = np.divide(lanes[i], max_intensity)

    for i, lane in enumerate(lanes):
        x1 = 50 + i * (lane_width+lanesep)
        x2 = x1 + lane_width
        for y, intensity in enumerate(lane):
            y1 = y
            y2 = y + 1
            draw.rectangle((x1, y1, x2, y2),
                           fill=(intensity,
                                 intensity,
                                 intensity))
        draw.rectangle((x1, 5, x2, start),
                       fill=(0, 0, 0),
                       outline=(256, 256, 256),
                       width=1)

    return image

    #  im.rotate(90, expand=1)
    #  im_invert = ImageOps.invert(im)


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
