#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


def test_drive(monkeypatch):

    from pydrive2.auth import GoogleAuth as _GoogleAuth
    from pydrive2.auth import (ServiceAccountCredentials
                               as _ServiceAccountCredentials)
    from pydrive2.drive import GoogleDrive as _GoogleDrive

    from unittest.mock import MagicMock, PropertyMock, seal

    # class mock_GoogleDrive(object):
    #     def __init__(self, auth=None):
    #         self.auth = auth

    mock_GA = MagicMock(spec=_GoogleAuth)
    mock_GD = MagicMock(spec=_GoogleDrive)
    mock_SAC = MagicMock(spec=_ServiceAccountCredentials)

    mock_GA.credentials = ""

    monkeypatch.setattr("pydna.myprimers_gdoc._GoogleAuth",
                        mock_GA)
    monkeypatch.setattr("pydna.myprimers_gdoc._GoogleDrive",
                        mock_GD)
    monkeypatch.setattr("pydna.myprimers_gdoc._ServiceAccountCredentials",
                        mock_SAC)

    from pydna import myprimers_gdoc

    primerlist = myprimers_gdoc.get_primer_list_from_gdoc()

    mock_GA.assert_called()
    assert mock_GA.credentials == ""

    mock_GD.assert_called()
    mock_SAC.assert_not_called()

    import os
    from pathlib import Path

    pth = Path(os.getenv("pydna_config_dir"))/"service_account.json"

    call = (str(pth),
            ['https://www.googleapis.com/auth/drive'])

    mock_SAC.from_json_keyfile_name.assert_called_with(*call)




if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
