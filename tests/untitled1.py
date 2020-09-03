#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on Wed Apr  8 09:54:33 2020 @author: bjorn

from textwrap import dedent

fig = r"""
        |98°C|98°C      |    |tmf:71.6
        |____|____      |    |tmr:75.3
        |30s |10s \ 72°C|72°C|15s/kb
        |    |     \____|____|GC 81%
        |    |      0: 0|5min|65bp
        """[
    1:
]

fig = dedent(fig)

fig = r"""
        |98°C|98°C               |    |tmf:32.6
        |____|_____          72°C|72°C|tmr:39.6
        |30s |10s  \ 35.6°C _____|____|15s/kb
        |    |      \______/ 0: 0|5min|GC 14%
        |    |       10s         |    |55bp
        """[
    1:
]

fig = dedent(fig)
