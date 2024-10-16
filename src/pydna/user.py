# -*- coding: utf-8 -*-
import re
from pydna.dseq import Dseq
from pydna.tm import tm_product
from typing import List


class USER:
    """
    A class representing a USER cloning reaction.
    In vitro, the USER enzyme recognizes dU bases in DNA sequences and creates
    an abasic site (single nucleotide gap). Then, an AP liase cleaves the DNA strand
    at the abasic site. Due to the high instability of the short oligonucleotide
    upstream of the cleavage site, the oligonucleotide detaches and leaves a 3' overhang.
    For optimal USER cloning, sequence upstream of the dU should be 6 to 10 nucleotides long.

    Parameters
    ----------
    size : int, optional
        The size of the pattern (default is 7).
    max_size : int, optional
        The maximum size of the pattern (default is 11).

    Examples
    --------

    >>> from pydna.user import USER
    >>> from pydna.dseq import Dseq
    >>> us=USER()
    >>> dna = Dseq("CGTCGCuCACACGT")
    >>>
    >>>

    """

    size = 6
    pattern = f"([ACGT]{{{size-1}}}U)"
    site = "N" * (size - 1) + "U"
    fst5 = size + 1  # First 5' cut
    fst3 = None
    ovhg = fst5 - 1  # This is a placeholder

    def __init__(self, size: int = 7, max_size: int = 11):
        """
        Initialize a USER enzyme with a pattern size.
        """
        self.size = size
        self.fst5 = size + 1
        # TODO: Properly implement max_size, requires change in search function (finditer vs match)
        self.max_size = max_size
        self.pattern = f"([ACGT]{{{size-1}}}U)"
        self.site = f"N{{{size-1}}}U"
        self.compsite = re.compile(f"(?=(?P<USER>{self.pattern}))", re.UNICODE)
        print("USER enzyme initialized with pattern size:", size)

    def search(self, dna, linear=True):
        """
        Search function for a USER enzyme that returns cut sites in the sense strand only.

        Parameters
        ----------
        dna : Dseq
            Dseq object representing the DNA sequence to search for USER site.
        linear : bool, optional
            If True, the search is performed on the input sequence.
            If False, the search is performed on the sequence + sequence[1:]. (default is True)

        Returns
        -------
        list
            A list of the positions of the USER target sites.

        """
        # TODO: Deal with circular DNA, should not be supported?
        # TODO: Should this support subsequent USER sites? i.e.
        #       0         10        20        30        40        50
        #       012345678901234567890123456789012345678901234567890
        #
        #       CGTCGCuTTCAGCACGTuGCTAGCGAGCGTAGTCTGACGTGCATC
        #       ------u1 -> u1 site recognized, cleaved and oligo detached
        #              ----------u2 -> u2 site recognized and cleaved. Once u1 is cleaved and the oligo detached, oligo between u1 and u2 is detached
        #       -----------------u2* -> If u1 is not cleaved, u2* is recognized, cleaved, but oligo may not detach
        #
        #       USER + AP liase should cleave any dU site
        #       However, spontaneous detachment of the oligonucleotide
        #       upstream depends on its length. For optimal USER cloning,
        #       it should be 6 to 10 nucleotides long.
        #       Should we support this weird use case? -> Finditer
        #       or should we just return the first cut site? -> Match
        #       >>> Dseq("CGTCGCuCACACGT")
        #
        results = []
        for mobj in self.compsite.finditer(dna.watson.upper()):
            cut = mobj.start() + self.fst5
            results.append(cut)

        return results

    def products(self, dseq: Dseq) -> List[Dseq]:
        """
        Generate USER products from the given Dseq.

        Parameters
        ----------
        dseq : Dseq
            The Dseq object representing the DNA sequence to generate USER products from.

        Returns
        -------
        List[Dseq]
            A list of Dseq objects representing the USER products.

        """
        results = []
        # USER sites in the forward strand
        for wcut in dseq.get_cutsites(self):
            wcut = wcut[0][0]  # Get position
            watson_user = dseq[wcut:]  # Get the sequence from the cutsite

            # USER sites in the reverse strand
            for ccut in Dseq(dseq.crick).get_cutsites(self):
                ccut = ccut[0][0]  # Get position
                crick_user = Dseq(dseq.crick)[ccut:]  # Get the sequence from the cutsite
                result = Dseq(str(watson_user), crick=str(crick_user), ovhg=wcut)

                # Calculate the Tm of the double stranded portion
                ds_portion = dseq[wcut:-ccut]
                tm = tm_product(ds_portion)
                results.append((result, tm))

        sorted_results = [x[0] for x in sorted(results, key=lambda x: x[1], reverse=True)]
        return sorted_results

    def __repr__(self):
        return f"ssUSER({self.site})"

    def __str__(self):
        return f"ssUSER({self.site})"


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
