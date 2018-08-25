#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import subprocess
subprocess      # shut up spyder!
from pydna.dseqrecord import Dseqrecord
from unittest import mock

def test_editor_wo_features(monkeypatch):
    from pydna import editor
    Popen = mock.MagicMock(name="subprocess.Popen")    
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    argument = Dseqrecord("ggatcc")
    editor.ape(argument)
    assert Popen.called

def test_editor_with_feature_label(monkeypatch):
    from pydna import editor
    Popen = mock.MagicMock(name="subprocess.Popen")    
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    argument = Dseqrecord("ggatcc")
    argument.add_feature(2,4, label = "lbl")
    editor.ape(argument)
    assert Popen.called

    
def test_editor_with_feature_note(monkeypatch):
    from pydna import editor
    Popen = mock.MagicMock(name="subprocess.Popen")    
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    argument = Dseqrecord("ggatcc")
    argument.add_feature(2,4, note = "nte")
    del argument.features[0].qualifiers["label"]
    editor.ape(argument)
    assert Popen.called



def test_editor_with_feature_wo_label_and_note(monkeypatch):
    from pydna import editor
    Popen = mock.MagicMock(name="subprocess.Popen")    
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    argument = Dseqrecord("ggatcc")
    argument.add_feature(2,4)
    del argument.features[0].qualifiers["label"]
    editor.ape(argument)
    assert Popen.called


if __name__ == '__main__':
    pytest.main([__file__, "-vv", "-s","--cov=pydna","--cov-report=html"])


