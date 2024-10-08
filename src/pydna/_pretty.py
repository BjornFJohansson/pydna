#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Classes for nicer Jupyter output.

The pretty_str class is similar to str but has a _repr_pretty_ method
for for nicer string output in the IPython shell and Jupyter notebook.
"""

from prettytable import PrettyTable as _Pt
from prettytable import MARKDOWN as _md

from copy import copy as _copy
from typing import List as _List


class pretty_str(str):
    """Thanks to Min RK, UC Berkeley for this."""

    def _repr_pretty_(self, p, cycle):
        p.text(self)


class PrettyTable(_Pt):
    """docstring."""

    def lol(self) -> _List[list]:
        """docstring."""
        return [self._field_names] + self._rows

    def __repr__(self) -> str:
        """docstring."""
        return self.get_string()

    def _repr_markdown_(self) -> pretty_str:
        c = _copy(self)
        c.set_style(_md)
        return pretty_str(c.get_string())


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
