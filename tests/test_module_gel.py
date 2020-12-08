import pytest


def test_gel():

    from pydna.gel import gel

    from pydna.ladders import PennStateLadder

    from pydna.dseqrecord import Dseqrecord

    mygel = gel([PennStateLadder, [Dseqrecord("A"*2000)]])

    from PIL import Image

    im = Image.open("pydna_gel.png")

    import numpy as np

    frame1 = np.asarray(mygel)
    frame2 = np.asarray(im)

    assert frame1.all() == frame2.all()


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
