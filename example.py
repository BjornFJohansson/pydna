# -*- coding: utf-8 -*-
import requests
from pydna.design import primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.tm import tm_default
import time


def primer_tm_neb(primer, conc=0.5, prodcode="q5-0"):
    """Calculates a single primers melting temp from NEB.

    Parameters
    ----------
    primer1 : str
    conc : float
    prodcode : str
        find product codes on nebswebsite: https://tmapi.neb.com/docs/productcodes

    Returns
    -------
    tm : int
        primer melting temperature

    """

    url = "https://tmapi.neb.com/tm"

    params = {"seq1": primer, "conc": conc, "prodcode": prodcode}

    # Mimic a browser request with this headers
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        "Accept-Encoding": "gzip, deflate, br",
        "Accept-Language": "en-US,en;q=0.5",
        "Connection": "keep-alive",
        "Referer": "https://example.com",
    }
    res = requests.get(url, params=params, headers=headers)
    r = res.json()
    print("making a request")
    if r["success"]:
        return r["data"]["tm1"]
    else:
        print("request failed")
        print(r["error"][0])


start_time = time.time()
result = primer_design(
    Dseqrecord("atgtcgtATGaaaccgttatcgatcatatgtGcgaaatgtcgcgcgtcatctacgtatcatcgatctactTAAacgtgta"),
    limit=8,
    target_tm=60,
    estimate_function=tm_default,
    tm_func=primer_tm_neb,
)
print("result", result.seq)
print("--- %s seconds ---" % (time.time() - start_time))


start_time = time.time()
result = primer_design(
    Dseqrecord("atgtcgtATGaaaccgttatcgatcatatgtGcgaaatgtcgcgcgtcatctacgtatcatcgatctactTAAacgtgta"),
    limit=8,
    target_tm=60,
    tm_func=primer_tm_neb,
)
print("result", result.seq)
print("--- %s seconds ---" % (time.time() - start_time))
