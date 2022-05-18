#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import gspread
import os
from pydrive2.auth import GoogleAuth
from pydrive2.auth import ServiceAccountCredentials
from pydrive2.drive import GoogleDrive

from pydna.myprimers import PrimerList
from pydna.myprimers import check_primer_numbers
from pydna.parsers import parse_primers

p = PrimerList()

check_primer_numbers()

def get_google_doc_text(title, mime="application/vnd.google-apps.document"):
    """Assumes that a google service account is used.

    See instructions at the docs for gdrive
    https://gspread.readthedocs.io/en/latest/oauth2.html
    """
    JSON_FILE = os.path.join( gspread.auth.get_config_dir(),
                             'service_account.json' )
    gauth = GoogleAuth()
    scope = ['https://www.googleapis.com/auth/drive']
    gauth.credentials = ServiceAccountCredentials.from_json_keyfile_name(JSON_FILE,
                                                                         scope)
    gauth.auth_method = 'service'
    drive = GoogleDrive(gauth)

    fl = drive.ListFile({'q': f"title = '{title}' and mimeType='{mime}'"}).GetList()

    return fl.pop(0).GetContentString(mimetype="text/plain")


primer_string = get_google_doc_text("PRIMERS")

p_gdoc = parse_primers(primer_string)
p_gdoc.reverse()

for p1,p2 in zip(p, p_gdoc):
    if str(p1.seq).lower() != str(p2.seq).lower():
        print(p1)
        print(p2)
        print("----")
