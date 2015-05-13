#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test empty
'''

import unittest

from pydna import *

class test_empty(unittest.TestCase):

    def test_MXblaster1(self):
        ''' test MXblaster1'''

        primer = parse("primers.fas", ds=False)
        primer     = primer[::-1]
        primer     = primer[37:]

        for i,p in enumerate(primer):
            assert int(p.id.split("_")[0]) == i

        ''' These are PCRs to get the genes and the terminator-promoters '''
        AgTEFp          = pcr(primer[524],primer[523],read("pAG25.gb"))
        hph             = pcr(primer[502],primer[501],read("pAG32.gb"))
        KlLEU2tt        = pcr(primer[520],primer[519],read("KlLEU2tt.gb"))

        ''' The Gal1 promoter-ISceI fragment is made in two steps '''
        gal1_ISceI_1    = pcr(primer[234],primer[316],read("pGSHU 7180bp.gb"))
        gal1_ISceI_2    = pcr(primer[562],primer[234],gal1_ISceI_1)
        AgTEFt          = pcr(primer[522],primer[521],read("pAG25.gb"))

        ''' load pCAPs and pCAPs-pSU0 sequences as Dseqrecord objects '''
        pCAPs      = read("pCAPs-AjiI.gb")
        pCAPs_pSU0 = read("pCAPs-pSU0.gb")

        # cut the pCAPs vectors for cloning
        from Bio.Restriction import EcoRV, ZraI
        pCAPs_ZraI              = pCAPs.cut(ZraI).pop()
        pCAPs_PCR_prod          = pcr(primer[492],primer[493], pCAPs)
        pCAPs_EcoRV             = pCAPs.cut(EcoRV).pop()
        pCAPs_pSU0_E_Z, stuffer = pCAPs_pSU0.cut(EcoRV, ZraI)


        # make the pCAPs clones, six altoghether
        pCAPs_ZraI_AgTEFp       = (pCAPs_ZraI + AgTEFp ).looped()

        pCAPs_PCR_prod_hph      = (pCAPs_PCR_prod  + hph ).looped()
        pCAPs_EcoRV_KlLEU2tt    = (pCAPs_EcoRV     + KlLEU2tt ).looped()

        pCAPs_ZraI_KlLEU2tt         = (pCAPs_ZraI      + KlLEU2tt ).looped()
        pCAPs_PCR_prod_gal1_ISceI_2 = (pCAPs_PCR_prod  + gal1_ISceI_2 ).looped()
        pCAPs_EcoRV_AgTEFt          = (pCAPs_EcoRV     + AgTEFt ).looped()

        # PCR each clone for the assembly in pCAPs

        A_AgTEFp_b     = pcr([primer[167],primer[493]], pCAPs_ZraI_AgTEFp)
        B_hph_c        = pcr([primer[467],primer[468]], pCAPs_PCR_prod_hph)
        C_KlLEU2tt_d   = pcr([primer[492],primer[166]], pCAPs_EcoRV_KlLEU2tt)

        #Homologous recombination of the two tp-gene-tp building blocks

        a=Assembly( ( A_AgTEFp_b,
                      B_hph_c,
                      C_KlLEU2tt_d,
                      pCAPs_pSU0_E_Z) , limit=28)

        YPK0_AgTEFp_hph_KlLEU2tt = a.circular_products[0]

        AgTEFp_hph_KlLEU2tt_2 = pcr(primer[166],primer[167], YPK0_AgTEFp_hph_KlLEU2tt)

        A_KlLEU2tt_b   = pcr([primer[167],primer[567]], pCAPs_ZraI_KlLEU2tt)
        B_gal1_ISceI_c = pcr([primer[467],primer[468]], pCAPs_PCR_prod_gal1_ISceI_2)
        C_AgTEFt_d     = pcr([primer[568],primer[166]], pCAPs_EcoRV_AgTEFt)

        a=Assembly(( A_KlLEU2tt_b,
                     B_gal1_ISceI_c,
                     C_AgTEFt_d,
                     pCAPs_pSU0_E_Z), limit=25)
        YPK0_KlLEU2tt_gal1_ISceI_AgTEFt = a.circular_products[0]

        KlLEU2tt_gal1_ISceI_AgTEFt_2 =  pcr(primer[166],primer[167], YPK0_KlLEU2tt_gal1_ISceI_AgTEFt)

        a=Assembly(( AgTEFp_hph_KlLEU2tt_2,
                     KlLEU2tt_gal1_ISceI_AgTEFt_2,
                     pCAPs_pSU0_E_Z), limit=61)

        pCAPs_MX4blaster1 = a.circular_products[0]

        pCAPs_MX4blaster1=pCAPs_MX4blaster1.synced("tcgcgcgtttcggtgatgacggtgaaaacc")

        self.assertTrue(  pCAPs_MX4blaster1.seguid()=="X9WqaNk2lw6FbZlJr995MaDfn-M" )

        from Bio.Restriction import AjiI, AgeI

        AX023560 = read("AX023560.gb")

        GAL10prom_slice= slice(AX023560.features[1].location.start,
                               AX023560.features[1].location.end)

        GAL10prom = AX023560[GAL10prom_slice]

        self.assertTrue( GAL10prom.seq == AX023560.features[1].extract(AX023560).seq )

        GIN11M86 = read("GIN11M86.gb")

        GAL_GIN = pcr(primer[592],primer[593], GAL10prom + GIN11M86)

        self.assertTrue( GAL_GIN.seguid() == "7Lkfw8dsz9_kkBU3XXnz4KAON3A" )

        self.assertTrue( pCAPs.seguid() =="-XHU8OxITyHGTl9XtMrJ4NvEv3o" )

        pCAPs_GAL_GIN = ( pCAPs.cut(AjiI).pop() + GAL_GIN ).looped()

        self.assertTrue( pCAPs_GAL_GIN.seguid() == "T1eWCPIXPlq2HriSfpFSNnGwmd4" )

        GAL_GIN2 = pcr(primer[592], primer[467], pCAPs_GAL_GIN)

        self.assertTrue( GAL_GIN2.seguid() =="zdIU4vjdfOxLkTTnKzIxhphnewg" )

        self.assertTrue( pCAPs_MX4blaster1.seguid() =="X9WqaNk2lw6FbZlJr995MaDfn-M" ) # 9772bp__a

        pCAPs_MX4blaster1_AgeI = pCAPs_MX4blaster1.cut(AgeI).pop()

        pCAPs_MX4blaster1_AgeI.seq = pCAPs_MX4blaster1_AgeI.seq.fill_in()

        a=Assembly([GAL_GIN2, pCAPs_MX4blaster1_AgeI], limit=30)

        pCAPs_MX4blaster2 = a.circular_products[0]

        pCAPs_MX4blaster2 = pCAPs_MX4blaster2.synced("tcgcgcgtttcggtgatgacggtgaaaacc")
        
        self.assertTrue( len(pCAPs_MX4blaster2) == 10566 )
        pCAPs_MX4blaster2_old = read("./pMX4blaster2_old.gb")
        self.assertTrue( len(pCAPs_MX4blaster2_old) == 10566 )
        self.assertEqual( pCAPs_MX4blaster2_old.seguid(), "7B4KKAeM2x8npjkp5U942rtMbB8" )
        self.assertTrue( eq(pCAPs_MX4blaster2, pCAPs_MX4blaster2_old))
        self.assertEqual( pCAPs_MX4blaster2.seguid(), "7B4KKAeM2x8npjkp5U942rtMbB8" )

if __name__ == '__main__':
    unittest.main()
