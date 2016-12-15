#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''Provides a class for downloading sequences from genbank.
'''
import pickle
import shelve
import os
import requests
import textwrap

def download_text(url):
    cached  = False
    refresh = False
    cache = shelve.open(os.path.join(os.environ["pydna_data_dir"], "web"), protocol=pickle.HIGHEST_PROTOCOL, writeback=False)
    key = str(url)

    if os.environ["pydna_cache"] in ("compare", "cached"):
        try:
            cached = cache[key]
        except KeyError:
            if os.environ["pydna_cache"] == "compare":
                raise Exception("no result for this key!")
            else:
                refresh = True

    if refresh or os.environ["pydna_cache"] in ("compare", "refresh", "nocache"):
        result = requests.get(url).text
    if os.environ["pydna_cache"] == "compare":
        if result!=cached:
            module_logger.warning('download error')

    if refresh or os.environ["pydna_cache"] == "refresh":
        cache = shelve.open(os.path.join(os.environ["pydna_data_dir"],"genbank"), protocol=pickle.HIGHEST_PROTOCOL, writeback=False)
        cache[key] = result

    elif cached and os.environ["pydna_cache"] not in ("nocache", "refresh"):
        result = cached

    cache.close()
    
    result = textwrap.dedent(result).strip()
    result = result.replace( '\r\n', '\n')
    result = result.replace( '\r',   '\n')
    return result

if __name__=="__main__":
    import os
    cache = os.getenv("pydna_cache")
    os.environ["pydna_cache"]="nocache"
    import doctest
    doctest.testmod()
    os.environ["pydna_cache"]=cache
