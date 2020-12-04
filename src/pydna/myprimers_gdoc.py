#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2020 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
""" docstring """

import os as _os
import copy as _copy
from pydna._pretty import pretty_str as _pretty_str
from pydrive2.auth import GoogleAuth as _GoogleAuth
from pydrive2.auth import (ServiceAccountCredentials
                           as _ServiceAccountCredentials)
from pydrive2.drive import GoogleDrive as _GoogleDrive
from pydna.parsers import parse_primers as _parse_primers


def drive(dir_=_os.environ["pydna_config_dir"],
          jsonfn="service_account.json"):
    """docstring."""
    JSON_FILE = _os.path.join(dir_, jsonfn)
    gauth = _GoogleAuth()
    scope = ['https://www.googleapis.com/auth/drive']
    gauth.credentials = (_ServiceAccountCredentials
                         .from_json_keyfile_name(JSON_FILE,
                                                 scope))
    gauth.auth_method = 'service'
    return _GoogleDrive(gauth)


def get_primer_list_from_gdoc(title=_os.environ["pydna_primersgdoc"],
                              mime="application/vnd.google-apps.document",
                              drive=drive):
    """Assumes that a google service account is used.

    See instructions at the docs for gdrive
    https://gspread.readthedocs.io/en/latest/oauth2.html

    In step 7, put the service_account.json file into the pydna config
    directory. You can find this directory by:

        pydna.open_config_folder()

    """
    fl = drive().ListFile(
        {'q': f"title = '{title}' and mimeType='{mime}'"}).GetList()

    content = fl.pop(0).GetContentString(mimetype="text/plain")
    return _parse_primers(content)[::-1]


def prepend_primer_list(primerlist,
                        title=_os.environ["pydna_primersgdoc"],
                        mime="application/vnd.google-apps.document",
                        drive=drive):
    """docstring."""

    primers = get_primer_list_from_gdoc(title, mime, drive())

    number_of_existing_primers = len(primers)

    newprimers = []

    for i,p in zip(range(number_of_existing_primers+len(primerlist)-1,
                         number_of_existing_primers-1,
                         -1),
                 primerlist):
        suff = p.id.split(str(i))[-1]
        suff.lstrip("_")
        newprimer = _copy.copy(p)
        newprimer.id = f"{i}_{suff}"
        newprimers.append(newprimer)

    return _pretty_str("\n".join([p.format("fasta") for p in newprimers]))


if __name__ == "__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"] = "nocache"
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cache"] = str(cache)
