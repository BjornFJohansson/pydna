#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test parse
"""

import pytest
from pydna.dseqrecord import Dseqrecord
from pydna.parsers import parse, parse_primers
from pydna.amplify import pcr, Anneal


def test_set_primer_footprint():
    f, r = parse_primers(
        """>ForwardPrimer
                            gctactacacacgtactgactg

                            >ReversePrimer
                            tgtggttactgactctatcttg"""
    )

    t = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    ampl = pcr((f, r), t)

    assert len(ampl.forward_primer.footprint) == 22
    assert len(ampl.reverse_primer.footprint) == 22

    ampl.set_forward_primer_footprint(15)
    ampl.set_reverse_primer_footprint(15)

    assert len(ampl.forward_primer.footprint) == 15
    assert len(ampl.reverse_primer.footprint) == 15


def test_string_arguments():
    f0, r0 = parse_primers(
        """>ForwardPrimer
                            gctactacacacgtactgactg

                            >ReversePrimer
                            tgtggttactgactctatcttg"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    f = str(f0.seq)
    r = str(r0.seq)
    t = str(t0.seq)

    assert str(pcr((f, r), t).seq) == "gctactacacacgtactgactgcctccaagatagagtcagtaaccaca"


def test_Seq_arguments():
    from Bio.Seq import Seq

    f0, r0 = parse_primers(
        """>ForwardPrimer
                            gctactacacacgtactgactg

                            >ReversePrimer
                            tgtggttactgactctatcttg"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    f = Seq(str(f0.seq))
    r = Seq(str(r0.seq))
    t = Seq(str(t0.seq))

    assert str(pcr(f, r, t).seq) == "gctactacacacgtactgactgcctccaagatagagtcagtaaccaca"


def test_Dseq_arguments():
    from pydna.dseq import Dseq

    f0, r0 = parse_primers(
        """>ForwardPrimer
                            gctactacacacgtactgactg

                            >ReversePrimer
                            tgtggttactgactctatcttg"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    f = Dseq(str(f0.seq))
    r = Dseq(str(r0.seq))
    t = Dseq(str(t0.seq))

    assert str(pcr(f, r, t).seq) == "gctactacacacgtactgactgcctccaagatagagtcagtaaccaca"


def test_wrong_argument_type():
    with pytest.raises(TypeError):
        pcr(1, 2, 3)


def test_no_primers_anneal():
    f0, r0 = parse_primers(
        """>ForwardPrimer
                             gctacta

                             >ReversePrimer
                             tgtggtt"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    f = f0
    r = r0
    t = t0

    with pytest.raises(ValueError):
        pcr(f, r, t)


def test_no_fwdprimer_anneal():
    f0, r0 = parse_primers(
        """>ForwardPrimer
                             gctact

                             >ReversePrimer
                             tgtggttactgactctatcttg"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    f = f0
    r = r0
    t = t0

    with pytest.raises(ValueError):
        pcr(f, r, t)


def test_no_revprimer_anneal():
    f0, r0 = parse_primers(
        """>ForwardPrimer
                             gctactacacacgtactgactg

                             >ReversePrimer
                             tgtggtt"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    f = f0
    r = r0
    t = t0

    with pytest.raises(ValueError):
        pcr(f, r, t)


def test_Primer_arguments():
    f0, r0 = parse_primers(
        """>ForwardPrimer
           gctactacacacgtactgactg

           >ReversePrimer
           tgtggttactgactctatcttg"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    f = f0
    r = r0
    t = t0

    assert str(pcr(f, r, t).seq) == "gctactacacacgtactgactgcctccaagatagagtcagtaaccaca"


def test_feature_label():
    f0, r0 = parse_primers(
        """>ForwardPrimer
           gctactacacacgtactgactg

           >ReversePrimer
           tgtggttactgactctatcttg"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")
    t0.add_feature()

    f = f0
    r = r0
    t = t0

    assert str(pcr(f, r, t).seq) == "gctactacacacgtactgactgcctccaagatagagtcagtaaccaca"


def test_feature_note():
    f0, r0 = parse_primers(
        """>ForwardPrimer
           gctactacacacgtactgactg

           >ReversePrimer
           tgtggttactgactctatcttg"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")
    t0.add_feature()
    del t0.features[0].qualifiers["label"]
    t0.features[0].qualifiers["note"] = ["note"]

    f = f0
    r = r0
    t = t0

    obj = pcr(f, r, t)
    assert str(obj.seq) == "gctactacacacgtactgactgcctccaagatagagtcagtaaccaca"
    assert obj.name == "note"


def test_Amplicon_argument():
    f0, r0 = parse_primers(
        """>ForwardPrimer
           gctactacacacgtactgactg

           >ReversePrimer
           tgtggttactgactctatcttg"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    f = f0
    r = r0
    t = t0

    ampl = pcr(f, r, t)

    assert str(ampl.seq) == "gctactacacacgtactgactgcctccaagatagagtcagtaaccaca"

    amplicon_from_amplicon = pcr(ampl)

    assert str(amplicon_from_amplicon.seq) == "gctactacacacgtactgactgcctccaagatagagtcagtaaccaca"


def test_pcr_not_specific():
    f0, r0 = parse_primers(
        """>ForwardPrimer
           gctactacacacgtactgactg

           >ReversePrimer
           tgtggttactgactctatcttg"""
    )

    t0 = Dseqrecord("gctactacacacgtactgactgtgctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    f = f0
    r = r0
    t = t0

    with pytest.raises(ValueError):
        pcr(f, r, t)


def test_too_short_primers():
    f, r = parse_primers(
        """>ForwardPrimer
           gctactacacacgtactgactg

           >ReversePrimer
           tgtggttactgactctatcttg"""
    )

    t = Dseqrecord("gctactacacacgtactgactgcctccaagatagagtcagtaaccaca")

    ann = Anneal((f, r), t, limit=22)

    assert ann.report() == (
        "Template name 48 bp linear limit=22:\n"
        "ForwardPrimer anneals forward (--->) at 22\n"
        "ReversePrimer anneals reverse (<---) at 26"
    )

    assert repr(ann) == "Reaction(products = 1)"

    p = ann.products[0]

    assert str(p.seq) == str(t.seq)

    ann = Anneal((f, r), t, limit=23)

    assert ann.products == []

    assert ann.report() == (
        "Template name 48 bp linear limit=23:\n" "No forward primers anneal...\n" "No reverse primers anneal..."
    )
    assert repr(ann) == "Reaction(products = 0)"


def test_circ_pcr():
    # A pathological example

    """
       <-----------------------------------------------------58->

       <----------------------------33->
         ccaccaccaccaccaccaccaccaccaccacggaggggccgtggtggacgagggcc
                                        |||||||||||||||||||||||||
    -->TCAGGTGAGGCGGAACCAACCCTCCTGGCCATGggaggggccgtggtggacgagggccccacaggcg-->
               ||||||||||||||||||||||||||
               ggcggaaccaaccctcctggccATGgactacaaagacgatgacgacaagcttgcggc

       <----------------------------------------------------------------------42->
       <-----------------------------------------------------25-><------------17->
         ccaccaccaccaccaccaccaccaccaccacggaggggccgtggtggacgagggcc
                                        |||||||||||||||||||||||||
                                        ggaggggccgtggtggacgagggccccacaggcgTCAGGTGAGGCGGAACCAACCCTCCTGGCCATG
                                                                                  ||||||||||||||||||||||||||
                                                                                  ggcggaaccaaccctcctggccATGgactacaaagacgatgacgacaagcttgcggc
    """

    s = parse(
        """
    >MCT4_flaghis_rv
    gccgcaagcttgtcgtcatcgtctttgtagtcCATggccaggagggttggttccgcc
    >MCT4_flaghis_fw
    ccaccaccaccaccaccaccaccaccaccacggaggggccgtggtggacgagggcc""",
        ds=False,
    )

    t = parse(
        """
    >hej circular
    TCAGGTGAGGCGGAACCAACCCTCCTGGCCATGggaggggccgtggtggacgagggccccacaggcg
    """
    )

    p = pcr(s, t)

    assert (
        str(p.seq).lower()
        == "ccaccaccaccaccaccaccaccaccaccacggaggggccgtggtggacgagggccccacaggcgTCAGGTGAGGCGGAACCAACCCTCCTGGCCATGgactacaaagacgatgacgacaagcttgcggc".lower()
    )


def test_pcr():
    """test pcr"""

    raw = []

    raw.append(
        (
            "ldseguid-q4nWrUFRjs0uWGk-byoK2NOpDQs",
            parse(
                """
    >524_pFA6aF (29-mer)
    cacatacgatttaggtgacactatagaac

    >523_AgTEF1tpR (21-mer)
    ggttgtttatgttcggatgtg
    """
            ),
            parse("pAG25.gb"),
        )
    )

    raw.append(
        (
            "ldseguid-1Ih52rBdWeFKdaF2psCkhbn0yoM",
            parse(
                """
    >lowgc_f
    TTTCACTAGTTACTTGTAGTCGacgtgccatctgtgcagacaaacgcatcaggatat

    >lowgc_r
    AAGTTGGAAATCTAGCTTTTCTTgacgtcagcggccgcattgcaca
    """
            ),
            parse("pCAPs.gb"),
        )
    )

    raw.append(
        (
            "ldseguid-q4nWrUFRjs0uWGk-byoK2NOpDQs",
            parse(
                """
    >524_pFA6aF (29-mer)
    cacatacgatttaggtgacactatagaac

    >523_AgTEF1tpR (21-mer)
    ggttgtttatgttcggatgtg
    """
            ),
            parse("pAG25.gb"),
        )
    )

    raw.append(
        (
            "ldseguid-2KtDVV5G_coqA0MkVhPB6FpMnUo",
            parse(
                """
    >ForwardPrimer1
    gctactacacacgtactgactg

    >ReversePrimer
    tgtggttactgactctatcttg

    LOCUS       sequence_50_bp            46 bp    DNA     circular UNK 08-FEB-2013
    DEFINITION  sequence_50_bp circular
    ACCESSION   sequence_50_bp
    VERSION     sequence_50_bp
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
    ORIGIN
            1 ccaagataga gtcagtaacc acagctacta cacacgtact gactgt
    //

    """
            ),
        )
    )

    raw.append(
        (
            "ldseguid-2KtDVV5G_coqA0MkVhPB6FpMnUo",
            parse(
                """
    >ForwardPrimer2
    gctactacacacgtactgactg

    >ReversePrimer
    tgtggttactgactctatcttg

    LOCUS       template                  46 bp    DNA     circular UNK 15-OCT-2012
    DEFINITION  template circular
    ACCESSION   template
    VERSION     template
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
    ORIGIN
            1 ccaagataga gtcagtaacc acagctacta cacacgtact gactgt
    //
    """
            ),
        )
    )
    raw.append(
        (
            "ldseguid-2KtDVV5G_coqA0MkVhPB6FpMnUo",
            parse(
                """
    >ForwardPrimer3
    gctactacacacgtactgactg

    >ReversePrimer
    tgtggttactgactctatcttg

    LOCUS       template                  46 bp    DNA     circular UNK 08-FEB-2013
    DEFINITION  template circular
    ACCESSION   template
    VERSION     template
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
    ORIGIN
            1 tccaagatag agtcagtaac cacagctact acacacgtac tgactg
    //
    """
            ),
        )
    )

    raw.append(
        (
            "ldseguid-2KtDVV5G_coqA0MkVhPB6FpMnUo",
            parse(
                """
    >ForwardPrimer4
    gctactacacacgtactgactg

    >ReversePrimer
    tgtggttactgactctatcttg

    LOCUS       template                  46 bp    DNA     circular UNK 15-OCT-2012
    DEFINITION  template circular
    ACCESSION   template
    VERSION     template
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
    ORIGIN
            1 gtccaagata gagtcagtaa ccacagctac tacacacgta ctgact
    //
    """
            ),
        )
    )

    raw.append(
        (
            "ldseguid--7vB-gVadO8mAUViEiYZWdCzPhE",
            parse(
                """
    >fw1
    cacatacgatttaggtgacactatagaac
    >rv
    ggttgtttatgttcggatgtg

    LOCUS       tm                        50 bp    DNA     circular UNK 15-OCT-2012
    DEFINITION  tm circular
    ACCESSION   tm
    VERSION     tm
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
    ORIGIN
            1 cacatccgaa cataaacaac ccacatacga tttaggtgac actatagaac
    //
    """
            ),
        )
    )

    raw.append(
        (
            "ldseguid--7vB-gVadO8mAUViEiYZWdCzPhE",
            parse(
                """
    >fw2
    cacatacgatttaggtgacactatagaac
    >rv
    ggttgtttatgttcggatgtg

    LOCUS       tm                        50 bp    DNA     circular UNK 15-OCT-2012
    DEFINITION  tm circular
    ACCESSION   tm
    VERSION     tm
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
    ORIGIN
            1 acatccgaac ataaacaacc cacatacgat ttaggtgaca ctatagaacc
    //

    """
            ),
        )
    )

    raw.append(
        (
            "ldseguid--7vB-gVadO8mAUViEiYZWdCzPhE",
            parse(
                """
    >fw3
    cacatacgatttaggtgacactatagaac
    >rv
    ggttgtttatgttcggatgtg
    LOCUS       tm                        50 bp    DNA     circular UNK 15-OCT-2012
    DEFINITION  tm circular
    ACCESSION   tm
    VERSION     tm
    KEYWORDS    .
    SOURCE      .
      ORGANISM  .
                .
    FEATURES             Location/Qualifiers
    ORIGIN
            1 ccacatccga acataaacaa cccacatacg atttaggtga cactatagaa
    //

    """
            ),
        )
    )

    raw.append(
        (
            "ldseguid-IPAo_O5u6Q3-5ToDBIIOLLC_I88",
            parse(
                """
    >f_Eric_Ma
    ARATGAGTCTTCTRACCGAGGTCG
    >r_Eric_Ma
    TGAAAAGACATCYTCAAGYYTCTG
    >templ
    AGCAAAAGCAGGTAGATATTGAAAAATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCGCAGAGACTTGAAGATGTCTCTGCAGGGAAGAACACTGATCTCGAGGCTCTCATGGAATGGCTAAAGACAAGACCAATCCTGTCACCTCTGACTAAGGGGATTTTAGGGTTTGTGTTCACGCTCACCGTGCCCAGTGAGCGAGGACTGCAGCGTAGACGCTTTGTCCAGAATGCCTTAAATGGGAATGGAGACCCAAACAACATGGACAGGGCAGTCAAACTATACAGGAAGCTGAAAAGAGAGATAACATTCCATGGGGCTAAAGAGGTTGCACTCAGCTATTCAACCGGTGCACTTGCCAGTTGCATGGGTCTCATATACAACAGGATGGGAACGGTAACCACAGAAGTAGCTTTTGGCCTGGTGTGTGCCACTTGTGAGCAGATTGCTGACTCACAGCATCGATCTCACAGACAGATGGTGACTACCACCAACCCACTAATCAGGCATGAAAACAGAATGGTGCTGGCCAGCACTACAGCTAAGGCTATGGAGCAGATGGCTGGATCGAGTGAACAGGCAGCGGAAGCCATGGAGGTTGCTAGTCAGGCTAGGCAGATGGTGCAGGCAATGAGGACAATTGGGACTCACCCTAGCTCCAGTGCCGGTCTGAAAGATGATCTTCTTGAAAATTTGCAGGCCTACCAGAAGCGGATGGGAGTGCAAATGCAGCGATTCAAGTGATCCTCTCGTTATTGCCGCAAGTATCATTGGGATCTTGCACTTGATATTGTGGATTCTTGATCGTCCTTTCTTCAAATGTATTTATCGTCGCCTTAAATACGGTTTGAAAAGAGGGCCTTCTACGGAAGGAGTGCCTGAGTCTATGAGGGAAGAGTATCGGCAGGAACAGCAGAGTGCTGTGGATGTTGACGATGGTCATTTTGTCAACATAGAGCTGGAGTAAAAAACTACCTTGTTTCTACT
    """
            ),
        )
    )

    for key, tst in enumerate(raw):
        assert tst[0] in pcr(tst[1:]).seguid()


def test_shifts():
    from pydna.parsers import parse
    from pydna.parsers import parse_primers
    from pydna.amplify import pcr

    # from pydna.amplify import nopcr

    s = parse(
        """
    >MCT4_flaghis_rv
    gccgcaagcttgtcgtcatcgtctttgtagtcCATggccaggagggttggttccgcc
    >MCT4_flaghis_fw
    ccaccaccaccaccaccaccaccaccaccacggaggggccgtggtggacgagggcc""",
        ds=False,
    )

    t = parse(
        """
    >hej circular
    TCAGGTGAGGCGGAACCAACCCTCCTGGCCATGggaggggccgtggtggacgagggccccacaggcg
    """
    )

    p = pcr(s, t)

    assert (
        str(p.seq).lower()
        == "ccaccaccaccaccaccaccaccaccaccacggaggggccgtggtggacgagggccccacaggcgtcaggtgaggcggaaccaaccctcctggccatggactacaaagacgatgacgacaagcttgcggc"
    )

    f, r, t = parse(
        """
    #A

    >ForwardPrimer
    actacacacgtactgactg

    >ReversePrimer
    ggttactgactctatcttg

    >MyTemplate
    gctactacacacgtactgactGcctcCaagatAgagtcagtaaccaca"""
    )
    a = pcr(f, r, t)
    assert str(a.seq).lower() == "actacacacgtactgactGcctcCaagatAgagtcagtaacc".lower()

    f, r, t = parse(
        """
    #B

    >ForwardPrimer
    actacacacgtactgactg

    >ReversePrimer
    ggttactgactctatcttg

    >MyTemplate circular
    gctactacacacgtactgactgcctccaagatagagtcagtaaccaca"""
    )

    #    actacacacgtactgactg>
    #    |||||||||||||||||||
    # gctactacacacgtactgactgcctccaagatagagtcagtaaccaca
    #
    # cgatgatgtgtgcatgactgacggaggttctatctcagtcattggtgt
    #                           |||||||||||||||||||
    #                          <gttctatctcagtcattgg

    b = pcr(f, r, t)

    assert a.seq == b.seq

    a.template = None
    b.template = None

    assert a == b

    f, r, t = parse(
        """
    #C

    >ForwardPrimer
    actacacacgtactgactg

    >ReversePrimer
    ggttactgactctatcttg

    >MyTemplate circular
    cgtactgactgcctccaagatagagtcagtaaccacagctactacaca"""
    )
    c = pcr(f, r, t)

    assert b.seq == c.seq

    f, r, t = parse(
        """
    #D

    >ForwardPrimer
    actacacacgtactgactg

    >ReversePrimer
    ggttactgactctatcttg

    >MyTemplate48 circular
    tccaagatagagtcagtaaccacagctactacacacgtactgactgcc"""
    )

    """
    012345678901234567890123456789012345678901234567
    ------------27--------------
                               actacacacgtactgactg
    tccaagatagagtcagtaaccacagctactacacacgtactgactgcc
      caagatagagtcagtaacc
    --
    012345678901234567890123456789012345678901234567

    ------------23----------
    actacacacgtactgactg
    actacacacgtactgactgcctccaagatagagtcagtaaccacagct
                           caagatagagtcagtaacc

    ------------------------------------------
    actacacacgtactgactgcctccaagatagagtcagtaacc

    """

    d = pcr(f, r, t)

    assert c.seq == d.seq

    f, r, t = parse(
        """
    #E

    >ForwardPrimer
    actacacacgtactgactg

    >ReversePrimer
    ggttactgactctatcttg

    >MyTemplate circular
    gagtcagtaaccacagctactacacacgtactgactGcctccaagata"""
    )
    e = pcr(f, r, t)

    assert d.seq == e.seq

    f, r, t = parse(
        """
    #F

    >ForwardPrimer
    actacacacgtactgactg

    >ReversePrimer
    ggttactgactctatcttg

    >MyTemplate circular
    actacacacgtactgactGcctccaagatagagtcagtaaccacagct"""
    )
    f = pcr(f, r, t)


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
