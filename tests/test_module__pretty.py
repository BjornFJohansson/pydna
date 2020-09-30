#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_repr_pretty_():
    from pydna._pretty import pretty_str
    from unittest.mock import MagicMock

    pp = MagicMock()
    s = pretty_str("123")
    t = s._repr_pretty_(pp, None)
    pp.text.assert_called_with("123")


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
