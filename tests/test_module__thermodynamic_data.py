#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_dHBr():
    from pydna import _thermodynamic_data

    assert hasattr(_thermodynamic_data.dHBr, "setdefault")
    assert len(_thermodynamic_data.dHBr) == 16


def test_dSBr():
    from pydna import _thermodynamic_data

    assert hasattr(_thermodynamic_data.dSBr, "setdefault")
    assert len(_thermodynamic_data.dSBr) == 16


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
