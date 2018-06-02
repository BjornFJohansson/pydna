#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import subprocess
subprocess      # shut up spyder!

from unittest import mock

def test_editor_wo_features(monkeypatch):
    Popen = mock.MagicMock(name="subprocess.Popen")    
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    
    from pydna.dseqrecord import Dseqrecord
    
    argument = Dseqrecord("ggatcc")


    from pydna import editor
    
    editor.ape(argument)
    Popen.assert_called()

def test_editor_with_feature_label(monkeypatch):
    Popen = mock.MagicMock(name="subprocess.Popen")    
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    
    from pydna.dseqrecord import Dseqrecord
    
    argument = Dseqrecord("ggatcc")
    argument.add_feature(2,4, label = "lbl")

    from pydna import editor
    
    editor.ape(argument)
    Popen.assert_called()
    
def test_editor_with_feature_note(monkeypatch):
    Popen = mock.MagicMock(name="subprocess.Popen")    
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    
    from pydna.dseqrecord import Dseqrecord
    
    argument = Dseqrecord("ggatcc")
    argument.add_feature(2,4, note = "nte")
    del argument.features[0].qualifiers["label"]
    from pydna import editor
    
    editor.ape(argument)
    Popen.assert_called()


def test_editor_with_feature_wo_label_and_note(monkeypatch):
    Popen = mock.MagicMock(name="subprocess.Popen")    
    monkeypatch.setattr("subprocess.Popen", Popen)
    monkeypatch.setenv("pydna_ape", "path/to/ape")
    
    from pydna.dseqrecord import Dseqrecord
    
    argument = Dseqrecord("ggatcc")
    argument.add_feature(2,4)
    del argument.features[0].qualifiers["label"]
    from pydna import editor
    
    editor.ape(argument)
    Popen.assert_called()
    
    
if __name__ == '__main__':
    pytest.main([__file__, "-vv", "-s","--cov=pydna","--cov-report=html"])


