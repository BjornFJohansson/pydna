#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright 2013-2018 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
import math as _math

def tmstaluc98(primer,*args, dnac=50, saltc=50, **kwargs):
    '''Returns the melting temperature (Tm) of the primer using
    the nearest neighbour algorithm. Formula and thermodynamic data
    is taken from SantaLucia 1998 [#]_. This implementation gives the same
    answer as the one provided by Biopython (See Examples).

    Thermodynamic data used:

    =====  ====  ====
    pair   dH    dS
    =====  ====  ====
    AA/TT  7.9   22.2
    AT/TA  7.2   20.4
    TA/AT  7.2   21.3
    CA/GT  8.5   22.7
    GT/CA  8.4   22.4
    CT/GA  7.8   21.0
    GA/CT  8.2   22.2
    CG/GC  10.6  27.2
    GC/CG  9.8   24.4
    GG/CC  8.0   19.9
    =====  ====  ====

    Parameters
    ----------
    primer : string
        Primer sequence 5'-3'

    Returns
    -------
    tm : float
        tm of the primer

    References
    ----------
    .. [#] SantaLucia J Jr. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proc Natl Acad Sci U S A 1998;95:1460–5.

    Examples
    --------

    >>> from pydna.tm import tmstaluc98
    >>> from Bio.SeqUtils.MeltingTemp import Tm_NN
    >>> tmstaluc98("ACGTCATCGACACTATCATCGAC")
    54.55597724052518
    >>> Tm_NN("ACGTCATCGACACTATCATCGAC")
    54.55597724052518


    '''

    nntermsl={  "AA": (7.9  , 22.2),
                "TT": (7.9  , 22.2),
                "AT": (7.2  , 20.4),
                "TA": (7.2  , 21.3),
                "CA": (8.5  , 22.7),
                "TG": (8.5  , 22.7),
                "GT": (8.4  , 22.4),
                "AC": (8.4  , 22.4),
                "CT": (7.8  , 21.0),
                "AG": (7.8  , 21.0),
                "GA": (8.2  , 22.2),
                "TC": (8.2  , 22.2),
                "CG": (10.6 , 27.2),
                "GC": (9.8  , 24.4),
                "GG": (8    , 19.9),
                "CC": (8    , 19.9),
                "A" : (0    , 0   ),
                "C" : (0    , 0   ),
                "G" : (0    , 0   ),
                "T" : (0    , 0   )  }

    helixinit = {   "G": (-0.1 ,2.8),
                    "C": (-0.1 ,2.8),
                    "A": (-2.3, -4.1),
                    "T": (-2.3, -4.1) }
    primer = primer.upper()
    dH, dS = helixinit[primer[0]]
    H ,  S = helixinit[primer[-1]]
    dH = dH+H
    dS = dS+S
    for p in range(len(primer)):
        dn = primer[p:p+2]
        H,S = nntermsl[dn]
        dH+=H
        dS+=S
    R = 1.987 # universal gas constant in Cal/degrees C*Mol
    k = (dnac/4.0)*1e-9
    dS = dS-0.368*(len(primer)-1)*_math.log(float(saltc)/1e3)
    tm = ((1000* (-dH))/(-dS+(R * (_math.log(k)))))-273.15
    return tm

def tmbreslauer86(primer, *args, dnac=500.0, saltc=50, thermodynamics=False, **kwargs):
    '''Returns the melting temperature (Tm) of the primer using
    the nearest neighbour algorithm. Formula and thermodynamic data
    is taken from Breslauer 1986.

    These data are no longer widely used.


    Breslauer 1986, table 2 [#]_

    =====  ===== ====   ===
    pair   dH    dS     dG
    =====  ===== ====   ===
    AA/TT  9.1   24.0   1.9
    AT/TA  8.6   23.9   1.5
    TA/AT  6.0   16.9   0.9
    CA/GT  5.8   12.9   1.9
    GT/CA  6.5   17.3   1.3
    CT/GA  7.8   20.8   1.6
    GA/CT  5.6   13.5   1.6
    CG/GC  11.9  27.8   3.6
    GC/CG  11.1  26.7   3.1
    GG/CC  11.0  26.6   3.1
    =====  ===== ====   ===

    Parameters
    ----------
    primer : string
        Primer sequence 5'-3'

    Returns
    -------
    tm : float


    References
    ----------
    .. [#] K.J. Breslauer et al., “Predicting DNA Duplex Stability from the Base Sequence,” Proceedings of the National Academy of Sciences 83, no. 11 (1986): 3746.


    Examples
    ---------

    >>> from pydna.tm import tmbreslauer86
    >>> tmbreslauer86("ACGTCATCGACACTATCATCGAC")
    64.28863985851899


    '''

    nntermbr={  "AA": (9.1   ,24.0   ,1.9),
                "TT": (9.1   ,24.0   ,1.9),
                "AT": (8.6   ,23.9   ,1.5),
                "TA": (6.0   ,16.9   ,0.9),
                "CA": (5.8   ,12.9   ,1.9),
                "TG": (5.8   ,12.9   ,1.9),
                "GT": (6.5   ,17.3   ,1.3),
                "AC": (6.5   ,17.3   ,1.3),
                "CT": (7.8   ,20.8   ,1.6),
                "AG": (7.8   ,20.8   ,1.6),
                "GA": (5.6   ,13.5   ,1.6),
                "TC": (5.6   ,13.5   ,1.6),
                "CG": (11.9  ,27.8   ,3.6),
                "GC": (11.1  ,26.7   ,3.1),
                "GG": (11.0  ,26.6   ,3.1),
                "CC": (11.0  ,26.6   ,3.1),
                "A" : (0     , 0     ,0),
                "C" : (0     , 0     ,0),
                "G" : (0     , 0     ,0),
                "T" : (0     , 0     ,0),     }
    dH=3.4
    dS=12.4
    dG=0
    primer = primer.upper()
    for p in range(len(primer)):
        dn = primer[p:p+2]
        H,S,G = nntermbr[dn]
        dG+=G
        dH+=H
        dS+=S

    R = 1.9872          # universal gas constant in Cal/degrees C*Mol
    k = dnac*1E-9/2.0
    dH = dH - 5
    dS = dS-0.368*(len(primer)-1)*_math.log(float(saltc)/1E3) # SantaLucia salt correction formula
    tm = 1000 * -dH /(-dS + R * _math.log(k) )  - 273.15 # degrees Celsius

    if thermodynamics:
        return tm,dH,dS
    else:
        return tm


def tmbresluc(primer, *args, primerc=500.0, saltc=50, thermodynamics=False, **kwargs):
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

    tm,dH,dS : tuple
        tm and dH and dS used for the calculation

    '''

    from . import _thermodynamic_data

    saltc = float(saltc)/1000
    pri  = primerc/10E7
    dS = -12.4
    dH = -3400

    STR = primer.lower();

    for i in range(len(STR)-1):
        n1=ord(STR[i])
        n2=ord(STR[i+1])
        dH += _thermodynamic_data.dHBr[n1 - 97][n2 - 97]
        dS += _thermodynamic_data.dSBr[n1 - 97][n2 - 97]

    tm = (dH / (1.9872 * _math.log(pri / 1600) + dS) + (16.6 * _math.log(saltc)) / _math.log(10)) - 273.15

    if thermodynamics:
        return tm,dH,dS
    else:
        return tm

def basictm(primer, *args, **kwargs):
    '''Returns the melting temperature (Tm) of the primer using
    the basic formula. This function returns the same value as
    the Biopython Bio.SeqUtils.MeltingTemp.Tm_Wallace

    | Tm = (wA+xT)*2 + (yG+zC)*4 assumed 50mM monovalent cations
    |
    | w = number of A in primer
    | x = number of T in primer
    | y = number of G in primer
    | z = number of C in primer

    Parameters
    ----------
    primer : string
        Primer sequence 5'-3'

    Returns
    -------
    tm : int

    Examples
    --------
    >>> from pydna.tm import basictm
    >>> basictm("ggatcc")
    20
    >>>

    '''
    primer = str(primer).lower()
    return (primer.count("a") + primer.count("t"))*2 + (primer.count("g") + primer.count("c"))*4

# http://www.promega.com/techserv/tools/biomath/calc11.htm#melt_results        
        

def Q5(*args,**kwargs):
    '''For Q5 Ta they take the lower of the two Tms and add 1C 
    (up to 72C). For Phusion they take the lower of the two 
    and add 3C (up to 72C). 
    '''
    raise NotImplementedError


if __name__=="__main__":
    import os as _os
    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"]=""
    import doctest
    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"]=cached
