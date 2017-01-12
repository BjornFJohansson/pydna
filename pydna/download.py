#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''Provides a class for downloading sequences from genbank.
'''
import logging as _logging
_module_logger = _logging.getLogger("pydna."+__name__)

from pydna.utils  import memorize as _memorize
from pydna._pretty import pretty_str as _pretty_str
import os        as _os
import requests  as _requests
import textwrap  as _textwrap           

@_memorize("download_text")
def download_text(url):
    _module_logger.info("#### DOWNLOAD TEXT ####")
    _module_logger.info("url = %s", url)
    req=_requests.get(url)
    _module_logger.info("url = %s", str(req))
    result = _textwrap.dedent( req.text ).strip()
    result = result.replace('\r\n', '\n').replace('\r', '\n')
    _module_logger.info("result[:160] = %s", result[:160])
    return _pretty_str(result)
    
if __name__=="__main__":
    cache = _os.getenv("pydna_cache")
    _os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod(verbose=True)
    _os.environ["pydna_cache"]=cache
    
    

