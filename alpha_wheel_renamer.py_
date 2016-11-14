#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wheel renamer
"""

import os, glob, pathlib

oldname = pathlib.Path( min( glob.glob("dist/*.whl") ) ).name

import re
posint = '(0|[1-9]\d*)'
tpl_string_re = ('^'
                 '({name}-)?'
                 '(0|[1-9]\d*\.)?'            # N
                 '(0|[1-9]\d*\.)?'        # (.N)*
                 '(0|[1-9]\d*)?'
                 '(a|b|rc)?(0|[1-9]\d*)?' # [{a|b|rc}N]
                 '(\.post(0|[1-9]\d*))?'  # [.postN]
                 '(\.dev{postdev})?'   # [.devN]
                 '(.*)'
                 '(\.whl)'
                 '$')
string_re = tpl_string_re.format(posint=posint, postdev=posint, name="pydna")
pep440re = re.compile(string_re)
m = pep440re.match(oldname)
print(m.groups())
print("".join(s or "" for s in m.groups()))

