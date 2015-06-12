#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pydna.download import Genbank as _Genbank
import sys
import os

import percache
import appdirs

data_dir = appdirs.user_data_dir("pydna_helper")
if not os.path.isdir(data_dir):
    os.mkdir(data_dir)
cache = percache.Cache( os.path.join( data_dir, "genbank" ) , livesync=True)

#print os.path.join( data_dir, "genbank" ) #/home/bjorn/.local/share/pydna_helper

def _get_proxy_from_global_settings():
    """Get proxy settings."""
    if sys.platform.startswith('linux'):
        try:
            from gi.repository import Gio
        except ImportError:
            return ''
        mode = Gio.Settings.new('org.gnome.system.proxy').get_string('mode')
        if mode == 'none' or mode == 'auto':
            return None
        http_settings = Gio.Settings.new('org.gnome.system.proxy.http')
        host = http_settings.get_string('host')
        port = http_settings.get_int('port')
        if http_settings.get_boolean('use-authentication'):
            username = http_settings.get_string('authentication_user')
            password = http_settings.get_string('authentication_password')
        else:
            username = password = None
            return 'http://{}:{}'.format(host, port)

#print _get_proxy_from_global_settings()
_gbloader = _Genbank("bjornjobb@gmail.com", proxy = _get_proxy_from_global_settings())

@cache
def gb(item):
    return _gbloader.nucleotide(item)
