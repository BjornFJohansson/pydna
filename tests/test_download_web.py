#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest

from pydna.download import download_text
    
def test_web():
    cachevar = os.environ["pydna_cache"]
    os.environ["pydna_cache"] = "cached"
    tx = download_text("http://www.textfiles.com/holiday/brownose.xms")    
    text = 'Rudolph the Brown Nose Reindeer\n\nAll the reindeer were milling around in the barn a few days before XMAS.\nSanta came out to check up on them and Rudolph began sucking up by\ntelling him how nice his new suit looked.  Santa says, "Well, thank you\nfor noticing, Rudolph" while eyeing the other reindeer. Rudolph then\nstarts complimenting Santa on how trim and slim he looks and how he must\nsurely have lost some weight and begun working out...again Santa thanks\nRudolph.  The other reindeer are muttering, "Sheesh, does this guy never\ngive up?" After Santa left, the other reindeer began badgering Rudolph\nfor his blatant kissing up.  Rudolph gives them a haughty look and says,\n"Well, just remember when we go flying what YOU GUYS will be looking\nat!"'
    assert tx==text    
    os.environ["pydna_cache"] = cachevar
    
    

if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])
    

    
