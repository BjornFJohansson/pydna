#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2013 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

"""
This module contain experimental functionality for smulation of agarose gel electrophoresis of DNA.

Note
----

This code is in a very early stage of development and documentation.
Read the source code to get an idea for how it works.

"""
try:
    import pygame
except ImportError:
    pass
import math
import webbrowser


def gamma(x, t, k=1):
    b = 1.0 / t
    x = x * k
    c = 0.0
    for i in range(0, k):
        c += (math.exp(-b * x) * (b * x) ** i) / math.factorial(i)
    return c


def testColour(c):
    """ Convert integer to colour, which saturates if > 255 """

    if c <= 255:
        return (c, c, c)
    else:
        return (255, 255, 255)


# def display(width, height, image):
#    """ Create a Pygame image of dimensions width x height and show image """
#
#    screen = pygame.display.set_mode((width, height))
#    screen.blit(image, (60, 30))
#    pygame.display.flip()
# running = True
# while running:
#    for event in pygame.event.get():
#        if event.type == pygame.QUIT:
#            running = False


class Gel:
    """A Gel contains lanes into which DNA can be loaded and run.
    The Gel is exposured to see the location and intensity of DNA."""

    LANE_WIDTH = 30.0
    LANE_HEIGHT = 6
    LANE_MARGIN = 6
    BORDER = 40

    def __init__(self, size, agarose=1):
        self.y = size[0]
        self.x = (
            2 * Gel.BORDER
            + size[1] * (Gel.LANE_WIDTH + Gel.LANE_MARGIN)
            - Gel.LANE_MARGIN
        )

        self.agarose = (
            agarose  # Should test whether value is reasonable (say, 0.5% - 3%)
        )
        self.optimum_DNA_length = (
            2000 / agarose ** 3
        )  # Best separate based on agarose hole size

        self.samples = []
        self.lanes = []

        for n in range(size[1]):
            self.samples.append([])
            self.lanes.append([])

        for lane in self.lanes:
            for n in range(Gel.BORDER, self.y + 3):
                lane.append(0.0)

    def loadSample(self, lane, sample):
        """Add list containing tuple of DNA length and concentrations to lane in gel. """

        for dna in sample:
            strand = {
                "length": float(dna[0]),
                "conc": dna[1] * 1.0 / Gel.LANE_HEIGHT,
                "position": Gel.BORDER + 1,
            }
            self.samples[lane - 1].append(strand)

    def run(self, time=30.0, voltage=20.0):
        """ Move loaded DNA down the gel at a rate dependent on voltage, DNA length and agarose concentration. """

        max_dist = 0.25 * time * voltage
        for sample in self.samples:
            for dna in sample:
                g = gamma(dna["length"] / 20, int(self.optimum_DNA_length / 20))
                dna["position"] += max_dist * g

    def expose(self, exposure=0.1, aperture=2):
        """Returns an image of the gel with the DNA highlighted"""

        c1, c2 = exposure * 100, exposure * 200
        fill_colour = (c1, c1, c2)
        image = pygame.Surface((self.x + 2, self.y + 2))
        pygame.draw.rect(image, fill_colour, (1, 1, self.x, self.y), 0)

        edge_colour = (140, 140, 150)
        edges = pygame.Surface((self.x + 2, self.y + 2))
        edges.set_colorkey((0, 0, 0))
        edges.set_alpha(120)
        pygame.draw.rect(image, edge_colour, (0, 0, self.x + 2, self.y + 2), 1)

        self._findDNAConcentrations(c1)

        x = Gel.BORDER + 1

        for lane in self.lanes:
            for position in range(len(lane)):
                if lane[position] > 0:

                    brightness = lane[position] * exposure / aperture + c1
                    colour1 = testColour(brightness)
                    colour2 = testColour(brightness * 0.6)
                    colour3 = testColour(brightness * 0.3)

                    # draw bands
                    pygame.draw.line(
                        image,
                        colour3,
                        (x - 1, position + Gel.BORDER),
                        (x + Gel.LANE_WIDTH + 1, position + Gel.BORDER),
                    )
                    pygame.draw.line(
                        image,
                        colour2,
                        (x, position + Gel.BORDER),
                        (x + Gel.LANE_WIDTH, position + Gel.BORDER),
                    )
                    pygame.draw.line(
                        image,
                        colour1,
                        (x + 1, position + Gel.BORDER),
                        (x + Gel.LANE_WIDTH - 1, position + Gel.BORDER),
                    )

            # draw well
            pygame.draw.rect(
                edges, edge_colour, (x, Gel.BORDER, Gel.LANE_WIDTH, Gel.LANE_HEIGHT), 1
            )
            x += Gel.LANE_WIDTH + Gel.LANE_MARGIN

        image.blit(edges, (0, 0))
        return image

    def _findDNAConcentrations(self, background):
        """Determines where in the concentration of DNA in every part of the gel"""

        length = len(self.lanes[0])

        for x in range(len(self.samples)):
            for dna in self.samples[x]:
                for y in range(Gel.LANE_HEIGHT - 2):
                    pos = int(dna["position"]) + y
                    if pos < length - 4:
                        # Very crude way to create blurred line
                        self.lanes[x][pos - 2] += 0.06 * dna["conc"] * dna["length"]
                        self.lanes[x][pos - 1] += 0.12 * dna["conc"] * dna["length"]
                        self.lanes[x][pos] += 0.2 * dna["conc"] * dna["length"]
                        self.lanes[x][pos + 1] += 0.12 * dna["conc"] * dna["length"]
                        self.lanes[x][pos + 2] += 0.06 * dna["conc"] * dna["length"]


def AgaroseGel(samples, *args, **kwargs):

    myGel = Gel(size=(360, len(samples)), agarose=1.0)

    for i, sample in enumerate(samples):
        try:
            sample = [(len(x), 100) for x in sample]
        except:
            pass
        myGel.loadSample(1 + i, sample)

    myGel.run(time=120)
    my_image = myGel.expose(exposure=0.01)
    pygame.image.save(my_image, "test_run.jpg")
    webbrowser.open("test_run.jpg")

    # display(360, 360, my_image)
    # from PIL import Image
    # img = Image.open('test_run.png')
    # img.show()
    # display(360, 360, my_image)
    # import sys;sys.exit()
    # Or save image as a PNG
    # import os
    # os.startfile('test_run.png')


if __name__ == "__main__":

    from pydna import *

    lambda_ = read("../tests/lambda.gb")
    from Bio.Restriction import PstI

    samples = sorted(lambda_.cut(PstI), key=len)
    AgaroseGel(
        [
            samples,
        ]
    )

    LADDER_1KB = [
        (300, 250),
        (500, 150),
        (700, 100),
        (1000, 75),
        (1500, 50),
        (2000, 40),
        (2500, 35),
        (3000, 30),
        (4000, 25),
        (5000, 30),
        (6000, 15),
        (8000, 12),
        (10000, 15),
        (14000, 4),
        (24000, 2),
    ]

    # Samples are list of tuples in the form (length, concentration)
    sample1 = [(400, 200), (500, 200), (900, 200)]
    sample2 = [(400, 200), (900, 200), (1500, 200)]
    sample3 = [(4000, 20), (4100, 20), (4900, 10)]

    # Set up gel, giving its size and % agarose
