#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.


# from abc import ABC, abstractmethod
# import re
# from pydna.utils import rc


"""
Nicking & USER cloning


CGAuGTCGACTTAGATCTCACAGGCTTTTTTCAAGaCGGCCTTGAATTCAGTCATTTGGATCCGGCCGATC
GCTACAGCTGAATCTAGAGTGTCCGAAAAAAGTTCTGCCGGAACTTAAGTCAGTAAACCTAGGCCGGCuAG

1. Digest both strands
2. Collect all linear ssDNA
3. Anneal all combinations
4. Keep ones present in original molecule
5. Rank by stability
6.

"""
