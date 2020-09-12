#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The pretty_str class is similar to str but has a _repr_pretty_ method for for nicer string output in the IPython shell"""


class pretty_str(str):
    """ Thanks to Min RK, UC Berkeley for this"""

    def _repr_pretty_(self, p, cycle):
        p.text(self)


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
