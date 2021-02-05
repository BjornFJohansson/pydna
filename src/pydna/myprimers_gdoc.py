#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""docstring.

You need pydrive2 ( https://github.com/iterative/PyDrive2 ).
"""
import os as _os
from pydrive2.auth import GoogleAuth as _GoogleAuth
from pydrive2.auth import (ServiceAccountCredentials
                           as _ServiceAccountCredentials)
from pydrive2.drive import GoogleDrive as _GoogleDrive
from pydna.parsers import parse_primers as _parse_primers


def get_primer_list_from_gdoc(title=_os.environ["pydna_primersgdoc"],
                              mime="application/vnd.google-apps.document",
                              dir_=_os.environ["pydna_config_dir"],
                              jsonfn="service_account.json"):
    """Assumes that a google service account is used.

    See instructions at the docs for gdrive
    https://gspread.readthedocs.io/en/latest/oauth2.html

    In step 7, put the service_account.json file into the pydna config
    directory. You can find this directory by:

        pydna.open_config_folder()

    """
    JSON_FILE = _os.path.join(dir_, jsonfn)
    scope = ['https://www.googleapis.com/auth/drive']
    gauth = _GoogleAuth()
    gauth.credentials = (_ServiceAccountCredentials
                         .from_json_keyfile_name(JSON_FILE,
                                                 scope))
    gauth.auth_method = 'service'

    fl = _GoogleDrive(gauth).ListFile(
        {'q': f"title = '{title}' and mimeType='{mime}'"}).GetList()

    content = fl.pop(0).GetContentString(mimetype="text/plain")

    lines = []
    for line in content.splitlines():
        if not line.startswith("#"):
            lines.append(line)

    return _parse_primers("\n".join(lines))[::-1]


if __name__ == "__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"] = str(cache)
