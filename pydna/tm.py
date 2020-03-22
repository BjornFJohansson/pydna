#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.

'''This module provide functions for melting temperature calculations.'''


import math as _math
from Bio.SeqUtils import MeltingTemp as _mt


def tm_default( seq,
                check=True,
                strict=True,
                c_seq=None,
                shift=0,
                nn_table=_mt.DNA_NN4,
                tmm_table=None,
                imm_table=None,
                de_table=None,
                dnac1=250,
                dnac2=250,
                selfcomp=False,
                Na=40,
                K=0,
                Tris=75.0,
                Mg=1.5,
                dNTPs=0.8,
                saltcorr=7,
                func = _mt.Tm_NN):
    return func(seq,
                check=check,
                strict=strict,
                c_seq=c_seq,
                shift=shift,
                nn_table=nn_table,
                tmm_table=tmm_table,
                imm_table=imm_table,
                de_table=de_table,
                dnac1=dnac1,
                dnac2=dnac2,
                selfcomp=selfcomp,
                Na=Na,
                K=K,
                Tris=Tris,
                Mg=Mg,
                dNTPs=dNTPs,
                saltcorr=saltcorr)


def tm_dbd( seq,
            check=True,
            strict=True,
            c_seq=None,
            shift=0,
            nn_table=_mt.DNA_NN3,
            tmm_table=None,
            imm_table=None,
            de_table=None,
            dnac1=250,
            dnac2=250,
            selfcomp=False,
            Na=50,
            K=0,
            Tris=0,
            Mg=1.5,
            dNTPs=0.8,
            saltcorr=1,
            func = _mt.Tm_NN):
    return func(seq,
                check=check,
                strict=strict,
                c_seq=c_seq,
                shift=shift,
                nn_table=nn_table,
                tmm_table=tmm_table,
                imm_table=imm_table,
                de_table=de_table,
                dnac1=dnac1,
                dnac2=dnac2,
                selfcomp=selfcomp,
                Na=Na,
                K=K,
                Tris=Tris,
                Mg=Mg,
                dNTPs=dNTPs,
                saltcorr=saltcorr)


def tm_product( seq,
                check=True,
                strict=True,
                valueset=7,
                userset=None,
                Na=50,
                K=0,
                Tris=0,
                Mg=0,
                dNTPs=0,
                saltcorr=0,
                mismatch=True,
                func=_mt.Tm_GC):
    return func(seq,
                check=check,
                strict=strict,
                valueset=valueset,
                userset=userset,
                Na=Na,
                K=K,
                Tris=Tris,
                Mg=Mg,
                dNTPs=dNTPs,
                saltcorr=saltcorr,
                mismatch=mismatch)


def ta_default(fp, rp, seq, tm=tm_default, tm_product=tm_product):
    # Ta calculation according to
    # Rychlik, Spencer, and Rhoads, 1990, Optimization of the anneal
    # ing temperature for DNA amplification in vitro
    # http://www.ncbi.nlm.nih.gov/pubmed/2243783
    # The formula described uses the length and GC content of the product and
    # salt concentration (monovalent cations)
    return 0.3*min((tm(fp),tm(rp)))+0.7*tm_product(seq)-14.9


def ta_dbd(fp, rp, seq, tm=tm_dbd, tm_product=None):
    # Ta calculation according to
    # Rychlik, Spencer, and Rhoads, 1990, Optimization of the anneal
    # ing temperature for DNA amplification in vitro
    # http://www.ncbi.nlm.nih.gov/pubmed/2243783
    # The formula described uses the length and GC content of the product and
    # salt concentration (monovalent cations)
    return min((tm(fp),tm(rp)))+3


def Q5(primer:str,*args,**kwargs):
    '''For Q5 Ta they take the lower of the two Tms and add 1C 
    (up to 72C). For Phusion they take the lower of the two 
    and add 3C (up to 72C). 
    '''
    raise NotImplementedError


def tmbresluc(primer:str, *args, primerc=500.0, saltc=50, **kwargs):
    '''Returns the tm for a primer using a formula adapted to polymerases
    with a DNA binding domain, such as the Phusion polymerase.

    Parameters
    ----------

    primer : string
        primer sequence 5'-3'

    primerc : float
       primer concentration in nM), set to 500.0 nm by default.

    saltc : float, optional
       Monovalent cation concentration in mM, set to 50.0 mM by default.

    thermodynamics : bool, optional
        prints details of the thermodynamic data to stdout. For
        debugging only.

    Returns
    -------
    tm : float
        the tm of the primer
        
    '''

    from . import _thermodynamic_data

    saltc = float(saltc)/1000
    pri  = primerc/1E9
    dS = -12.4
    dH = -3400

    STR = primer.lower();

    for i in range(len(STR)-1):
        n1=ord(STR[i])
        n2=ord(STR[i+1])
        dH += _thermodynamic_data.dHBr[n1 - 97][n2 - 97]
        dS += _thermodynamic_data.dSBr[n1 - 97][n2 - 97]

    tm = (dH / (1.9872 * _math.log(pri / 1600) + dS) + (16.6 * _math.log(saltc)) / _math.log(10)) - 273.15

    return tm 
 

if __name__=="__main__":
    import os as _os
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"]=""
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"]=cached
