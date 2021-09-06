import pytest


def test_Seq():

    from pydna.seq import Seq

    from pydna.codon import rare_codons

    assert Seq("atgtaa").gc() == 0.167

    assert Seq("atgtaa").cai() == 1.0

    lol = {}

    lol["sce"] = [['cds',
                   'len',
                   'cai',
                   'gc',
                   'sta',
                   'stp',
                   'n-end',
                   'CGA',
                   'CGG',
                   'CGC',
                   'CCG',
                   'CTC',
                   'GCG',
                   'rare'],
                  ['ATG...TAA', 8.0, 0.219, 0.708, 1.0, 0.47,
                   '2 min', 1, 1, 1, 1, 1, 1, 0.75]]

    lol["eco"] = [['cds',
                   'len',
                   'cai',
                   'gc',
                   'sta',
                   'stp',
                   'n-end',
                   'CGA',
                   'CGG',
                   'CGC',
                   'CCG',
                   'CTC',
                   'GCG',
                   'rare'],
                  ['ATG...TAA', 10.0, 0.387, 0.5, 1.0, 0.47,
                   '2 min', 1, 1, 0, 0, 0, 0, 0.2]]

    for organism, codons in rare_codons.items():
        s = Seq("atg"+"".join(codons)+"taa")
        slices = s.rarecodons(organism=organism)
        for slc in slices:
            assert s[slc].upper() in codons
        assert s.express().lol() == lol[organism]

    assert Seq("atgtaa").startcodon() == 1.0
    assert Seq("atgtaa").stopcodon() == 0.47

    s = Seq("atgaaattttaa")

    assert s.orfs(2) == [Seq('atgaaattttaa')]

    s2 = Seq("ctgaaattttaa")

    assert s2.orfs(2) == []
    
    

if __name__ == "__main__":
    args = [
    __file__,
    "--cov=pydna",
    "--cov-append",
    "--cov-report=html:../htmlcov",
    "--cov-report=xml",
    "--capture=no",
    "--durations=10",
    "--import-mode=importlib",
    "--nbval",
    "--current-env",
    "--doctest-modules",
    "--capture=no",
    "-vvv"]    
    pytest.main(args)


#    >>> import warnings
#    >>> from pydna import _PydnaWarning
#    >>> warnings.simplefilter('ignore', _PydnaWarning)

# obj.name = "1111111111111111"


#    def __init__(self, seq, id="<unknown id>", name="<unknown name>",
#                 description="<unknown description>", dbxrefs=None,
#                 features=None, annotations=None,
#                 letter_annotations=None):
#        """Create a SeqRecord.
