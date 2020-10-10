#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from unittest import mock


def test_editor_wo_features(monkeypatch):
    import subprocess
    from pydna import editor

    Popen = mock.MagicMock(name="subprocess.Popen")
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    from pydna.dseqrecord import Dseqrecord

    argument = Dseqrecord("ggatcc")
    editor.ape(argument)
    assert Popen.called


def test_editor_with_feature_label(monkeypatch):
    import subprocess
    from pydna import editor

    Popen = mock.MagicMock(name="subprocess.Popen")
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    from pydna.dseqrecord import Dseqrecord

    argument = Dseqrecord("ggatcc")
    argument.add_feature(2, 4, label="lbl")
    editor.ape(argument)
    assert Popen.called


def test_editor_with_feature_note(monkeypatch):
    import subprocess
    from pydna import editor

    Popen = mock.MagicMock(name="subprocess.Popen")
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    from pydna.dseqrecord import Dseqrecord

    argument = Dseqrecord("ggatcc")
    argument.add_feature(2, 4, note="nte")
    del argument.features[0].qualifiers["label"]
    editor.ape(argument)
    assert Popen.called


def test_editor_with_feature_wo_label_and_note(monkeypatch):
    import subprocess
    from pydna import editor

    Popen = mock.MagicMock(name="subprocess.Popen")
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    from pydna.dseqrecord import Dseqrecord

    argument = Dseqrecord("ggatcc")
    argument.add_feature(2, 4)
    del argument.features[0].qualifiers["label"]
    editor.ape(argument)
    assert Popen.called


if __name__ == "__main__":
    pytest.main([__file__, "-x", "-vv", "-s"])
