import pytest


def test_fakeseq():

    from pydna.fakeseq import FakeSeq

    fs = FakeSeq(1000)

    assert fs.m() == 3.08979e-08

    assert len(fs) == 1000

    assert repr(fs) == 'FakeSeq(1000)'

    assert fs > FakeSeq(999)


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
