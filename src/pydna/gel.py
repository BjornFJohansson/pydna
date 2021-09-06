#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# doctest: +NORMALIZE_WHITESPACE
# doctest: +SKIP
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""docstring."""

from PIL import Image as _Image
from PIL import ImageDraw as _ImageDraw
import numpy as _np
import math as _math
from scipy.interpolate import CubicSpline as _CubicSpline
from pydna.ladders import GeneRuler_1kb_plus as _mwstd


def interpolator(mwstd):
    """docstring."""
    interpolator = _CubicSpline([len(fr) for fr in mwstd[::-1]],
                                [fr.rf for fr in mwstd[::-1]],
                                bc_type="natural",
                                extrapolate=False)
    interpolator.mwstd = mwstd
    return interpolator


def gel(samples=None,
        gel_length=600,
        margin=50,
        interpolator=interpolator(mwstd=_mwstd)):
    """docstring."""
    max_intensity = 256
    lane_width = 50
    lanesep = 10
    start = 10
    samples = samples or [interpolator.mwstd]
    width = int(60 + (lane_width + lanesep) * len(samples))
    lanes = _np.zeros((len(samples), gel_length), dtype=int)
    image = _Image.new("RGB", (width, gel_length), "#ddd")
    draw = _ImageDraw.Draw(image)
    draw.rectangle((0, 0, (width, gel_length)), fill=(0, 0, 0))
    scale = (gel_length - margin) / interpolator(min(interpolator.x))

    for labelsource in samples[0]:
        peak_centre = (interpolator(len(labelsource))) * scale - 5 + start
        label = f"{len(labelsource):<5} -"
        draw.text((2, peak_centre), label, fill=(255, 255, 255))

    for lane_number, lane in enumerate(samples):
        for band in lane:
            log = _math.log(len(band), 10)
            height = (band.m() / (240 * log)) * 1e10
            peak_centre = interpolator(len(band)) * scale + start
            max_spread = 10
            if len(band) <= 50:
                peak_centre += 50
                max_spread *= 4
                max_intensity /= 10
            band_spread = max_spread / log
            for i in range(max_spread, 0, -1):
                y1 = peak_centre - i
                y2 = peak_centre + i
                intensity = (
                    height * _math.exp(
                        -float(((y1 - peak_centre) ** 2)) / (2 * (band_spread ** 2))
                    ) * max_intensity)
                for y in range(int(y1), int(y2)):
                    try:
                        lanes[lane_number][y] += intensity
                    except IndexError:
                        pass

    for i, lane in enumerate(lanes):
        max_intensity = _np.amax(lanes[i])
        if max_intensity > 256:
            lanes[i] = _np.multiply(lanes[i], 256)
            lanes[i] = _np.divide(lanes[i], max_intensity)

    for i, lane in enumerate(lanes):
        x1 = 50 + i * (lane_width + lanesep)
        x2 = x1 + lane_width
        for y, intensity in enumerate(lane):
            y1 = y
            y2 = y + 1
            draw.rectangle((x1, y1, x2, y2), fill=(intensity,
                                                   intensity,
                                                   intensity))

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
