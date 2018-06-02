#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2013 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

from pydna.dseqrecord import Dseqrecord as _Dseqrecord

class GenbankFile(_Dseqrecord):

    def __init__(self, record, *args, path=None, **kwargs):
        super().__init__(record, *args, **kwargs)
        self.path=path

    def __repr__(self):
        '''returns a short string representation of the object'''
        return "File({})({}{})".format(self.id, {True:"-", False:"o"}[self.linear],len(self))
        
    def _repr_pretty_(self, p, cycle):
        '''returns a short string representation of the object'''
        p.text("File({})({}{})".format(self.id, {True:"-", False:"o"}[self.linear],len(self)))
            
    def _repr_html_(self):
        return "<a href='{path}' target='_blank'>{path}</a><br>".format(path=self.path)

    def reverse_complement(self):
        answer = type(self)(super().reverse_complement(), path = self.path) 
        return answer
    
    rc = reverse_complement

if __name__=="__main__":
#    import os as _os
#    cached = _os.getenv("pydna_cached_funcs", "")
#    _os.environ["pydna_cached_funcs"]=""
#    import doctest
#    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
#    _os.environ["pydna_cached_funcs"]=cached   
    from pydna.readers import read

    a=read('''
LOCUS       PRS416                  4898 bp    DNA     circular SYN 24-MAY-1995
DEFINITION  Yeast centromere vector pRS416 with URA3 marker, complete sequence.
ACCESSION   U03450
VERSION     U03450.1
KEYWORDS    .
SOURCE      Cloning vector pRS416
  ORGANISM  Cloning vector pRS416
            other sequences; artificial sequences; vectors.
REFERENCE   1  (bases 1 to 4898)
  AUTHORS   Sikorski,R.S. and Hieter,P.
  TITLE     A system of shuttle vectors and yeast host strains designed for
            efficient manipulation of DNA in Saccharomyces cerevisiae
  JOURNAL   Genetics 122 (1), 19-27 (1989)
   PUBMED   2659436
REFERENCE   2  (bases 1 to 4898)
  AUTHORS   Christianson,T.W., Sikorski,R.S., Dante,M., Shero,J.H. and
            Hieter,P.
  TITLE     Multifunctional yeast high-copy-number shuttle vectors
  JOURNAL   Gene 110 (1), 119-122 (1992)
   PUBMED   1544568
REFERENCE   3  (bases 1 to 4898)
  AUTHORS   Stillman,D.J.
  TITLE     Direct Submission
  JOURNAL   Submitted (11-NOV-1993) David J. Stillman, Dept. of Cellular, Viral
            and Molecular Biology, University of Utah Medical Center, Salt Lake
            City, UT 84132 USA
FEATURES             Location/Qualifiers
     source          1..4898
                     /organism="Cloning vector pRS416"
                     /mol_type="genomic DNA"
                     /db_xref="taxon:31841"
ORIGIN      
        1 tcgcgcgttt cggtgatgac ggtgaaaacc tctgacacat gcagctcccg gagacggtca
       61 cagcttgtct gtaagcggat gccgggagca gacaagcccg tcagggcgcg tcagcgggtg
      121 ttggcgggtg tcggggctgg cttaactatg cggcatcaga gcagattgta ctgagagtgc
      181 accataccac agcttttcaa ttcaattcat catttttttt ttattctttt ttttgatttc
      241 ggtttctttg aaattttttt gattcggtaa tctccgaaca gaaggaagaa cgaaggaagg
      301 agcacagact tagattggta tatatacgca tatgtagtgt tgaagaaaca tgaaattgcc
      361 cagtattctt aacccaactg cacagaacaa aaacctgcag gaaacgaaga taaatcatgt
      421 cgaaagctac atataaggaa cgtgctgcta ctcatcctag tcctgttgct gccaagctat
      481 ttaatatcat gcacgaaaag caaacaaact tgtgtgcttc attggatgtt cgtaccacca
      541 aggaattact ggagttagtt gaagcattag gtcccaaaat ttgtttacta aaaacacatg
      601 tggatatctt gactgatttt tccatggagg gcacagttaa gccgctaaag gcattatccg
      661 ccaagtacaa ttttttactc ttcgaagaca gaaaatttgc tgacattggt aatacagtca
      721 aattgcagta ctctgcgggt gtatacagaa tagcagaatg ggcagacatt acgaatgcac
      781 acggtgtggt gggcccaggt attgttagcg gtttgaagca ggcggcagaa gaagtaacaa
      841 aggaacctag aggccttttg atgttagcag aattgtcatg caagggctcc ctatctactg
      901 gagaatatac taagggtact gttgacattg cgaagagcga caaagatttt gttatcggct
      961 ttattgctca aagagacatg ggtggaagag atgaaggtta cgattggttg attatgacac
     1021 ccggtgtggg tttagatgac aagggagacg cattgggtca acagtataga accgtggatg
     1081 atgtggtctc tacaggatct gacattatta ttgttggaag aggactattt gcaaagggaa
     1141 gggatgctaa ggtagagggt gaacgttaca gaaaagcagg ctgggaagca tatttgagaa
     1201 gatgcggcca gcaaaactaa aaaactgtat tataagtaaa tgcatgtata ctaaactcac
     1261 aaattagagc ttcaatttaa ttatatcagt tattacccta tgcggtgtga aataccgcac
     1321 agatgcgtaa ggagaaaata ccgcatcagg aaattgtaaa cgttaatatt ttgttaaaat
     1381 tcgcgttaaa tttttgttaa atcagctcat tttttaacca ataggccgaa atcggcaaaa
     1441 tcccttataa atcaaaagaa tagaccgaga tagggttgag tgttgttcca gtttggaaca
     1501 agagtccact attaaagaac gtggactcca acgtcaaagg gcgaaaaacc gtctatcagg
     1561 gcgatggccc actacgtgaa ccatcaccct aatcaagttt tttggggtcg aggtgccgta
     1621 aagcactaaa tcggaaccct aaagggagcc cccgatttag agcttgacgg ggaaagccgg
     1681 cgaacgtggc gagaaaggaa gggaagaaag cgaaaggagc gggcgctagg gcgctggcaa
     1741 gtgtagcggt cacgctgcgc gtaaccacca cacccgccgc gcttaatgcg ccgctacagg
     1801 gcgcgtcgcg ccattcgcca ttcaggctgc gcaactgttg ggaagggcga tcggtgcggg
     1861 cctcttcgct attacgccag ctggcgaaag ggggatgtgc tgcaaggcga ttaagttggg
     1921 taacgccagg gttttcccag tcacgacgtt gtaaaacgac ggccagtgag cgcgcgtaat
     1981 acgactcact atagggcgaa ttgggtaccg ggccccccct cgaggtcgac ggtatcgata
     2041 agcttgatat cgaattcctg cagcccgggg gatccactag ttctagagcg gccgccaccg
     2101 cggtggagct ccagcttttg ttccctttag tgagggttaa ttgcgcgctt ggcgtaatca
     2161 tggtcatagc tgtttcctgt gtgaaattgt tatccgctca caattccaca caacatagga
     2221 gccggaagca taaagtgtaa agcctggggt gcctaatgag tgaggtaact cacattaatt
     2281 gcgttgcgct cactgcccgc tttccagtcg ggaaacctgt cgtgccagct gcattaatga
     2341 atcggccaac gcgcggggag aggcggtttg cgtattgggc gctcttccgc ttcctcgctc
     2401 actgactcgc tgcgctcggt cgttcggctg cggcgagcgg tatcagctca ctcaaaggcg
     2461 gtaatacggt tatccacaga atcaggggat aacgcaggaa agaacatgtg agcaaaaggc
     2521 cagcaaaagg ccaggaaccg taaaaaggcc gcgttgctgg cgtttttcca taggctccgc
     2581 ccccctgacg agcatcacaa aaatcgacgc tcaagtcaga ggtggcgaaa cccgacagga
     2641 ctataaagat accaggcgtt tccccctgga agctccctcg tgcgctctcc tgttccgacc
     2701 ctgccgctta ccggatacct gtccgccttt ctcccttcgg gaagcgtggc gctttctcat
     2761 agctcacgct gtaggtatct cagttcggtg taggtcgttc gctccaagct gggctgtgtg
     2821 cacgaacccc ccgttcagcc cgaccgctgc gccttatccg gtaactatcg tcttgagtcc
     2881 aacccggtaa gacacgactt atcgccactg gcagcagcca ctggtaacag gattagcaga
     2941 gcgaggtatg taggcggtgc tacagagttc ttgaagtggt ggcctaacta cggctacact
     3001 agaaggacag tatttggtat ctgcgctctg ctgaagccag ttaccttcgg aaaaagagtt
     3061 ggtagctctt gatccggcaa acaaaccacc gctggtagcg gtggtttttt tgtttgcaag
     3121 cagcagatta cgcgcagaaa aaaaggatct caagaagatc ctttgatctt ttctacgggg
     3181 tctgacgctc agtggaacga aaactcacgt taagggattt tggtcatgag attatcaaaa
     3241 aggatcttca cctagatcct tttaaattaa aaatgaagtt ttaaatcaat ctaaagtata
     3301 tatgagtaaa cttggtctga cagttaccaa tgcttaatca gtgaggcacc tatctcagcg
     3361 atctgtctat ttcgttcatc catagttgcc tgactccccg tcgtgtagat aactacgata
     3421 cgggagggct taccatctgg ccccagtgct gcaatgatac cgcgagaccc acgctcaccg
     3481 gctccagatt tatcagcaat aaaccagcca gccggaaggg ccgagcgcag aagtggtcct
     3541 gcaactttat ccgcctccat ccagtctatt aattgttgcc gggaagctag agtaagtagt
     3601 tcgccagtta atagtttgcg caacgttgtt gccattgcta caggcatcgt ggtgtcacgc
     3661 tcgtcgtttg gtatggcttc attcagctcc ggttcccaac gatcaaggcg agttacatga
     3721 tcccccatgt tgtgcaaaaa agcggttagc tccttcggtc ctccgatcgt tgtcagaagt
     3781 aagttggccg cagtgttatc actcatggtt atggcagcac tgcataattc tcttactgtc
     3841 atgccatccg taagatgctt ttctgtgact ggtgagtact caaccaagtc attctgagaa
     3901 tagtgtatgc ggcgaccgag ttgctcttgc ccggcgtcaa tacgggataa taccgcgcca
     3961 catagcagaa ctttaaaagt gctcatcatt ggaaaacgtt cttcggggcg aaaactctca
     4021 aggatcttac cgctgttgag atccagttcg atgtaaccca ctcgtgcacc caactgatct
     4081 tcagcatctt ttactttcac cagcgtttct gggtgagcaa aaacaggaag gcaaaatgcc
     4141 gcaaaaaagg gaataagggc gacacggaaa tgttgaatac tcatactctt cctttttcaa
     4201 tattattgaa gcatttatca gggttattgt ctcatgagcg gatacatatt tgaatgtatt
     4261 tagaaaaata aacaaatagg ggttccgcgc acatttcccc gaaaagtgcc acctgggtcc
     4321 ttttcatcac gtgctataaa aataattata atttaaattt tttaatataa atatataaat
     4381 taaaaataga aagtaaaaaa agaaattaaa gaaaaaatag tttttgtttt ccgaagatgt
     4441 aaaagactct agggggatcg ccaacaaata ctacctttta tcttgctctt cctgctctca
     4501 ggtattaatg ccgaattgtt tcatcttgtc tgtgtagaag accacacacg aaaatcctgt
     4561 gattttacat tttacttatc gttaatcgaa tgtatatcta tttaatctgc ttttcttgtc
     4621 taataaatat atatgtaaag tacgcttttt gttgaaattt tttaaacctt tgtttatttt
     4681 tttttcttca ttccgtaact cttctacctt ctttatttac tttctaaaat ccaaatacaa
     4741 aacataaaaa taaataaaca cagagtaaat tcccaaatta ttccatcatt aaaagatacg
     4801 aggcgcgtgt aagttacagg caagcgatcc gtcctaagaa accattatta tcatgacatt
     4861 aacctataaa aataggcgta tcacgaggcc ctttcgtc
//

    ''')
    
    b=GenbankFile(a)
    
    print( b )