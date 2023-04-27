#                  _
#  _ __  _   _  __| |_ __   __ _
# | '_ \| | | |/ _` | '_ \ / _` |
# | |_) | |_| | (_| | | | | (_| |
# | .__/ \__, |\__,_|_| |_|\__,_|
# |_|    |___/
#
#                            ..,.
#                  .*/((((((((###((((((((/*
#              *((((#. ,                 *#(((/*
#           *(((( ,*    ,                     #(((*
#        *(((/,***(                         (%%%*#((/.
#      ,((# ,**(                               &%%(.((/
#     /(( ,**/              #@*                  &%  *((*
#    ((#.**(           ( .@ &  (,./*                   #(/
#   ((#,***          &***@&,/%*,,,,,/                   #(/
#  /(#.**/           /****,,*(#,,,,,%                (   ((*
# *(( ,*(          &         %***,*                      (((
# *(( ***        #..        %****@      ,(@/.             ((*
# *(#.**,                  /****/*(@%***#%                ((*
# *(( **/              %*(/////*******@                   ((*
# ,((/,*/            *//*#/////////*&.,,                 #((
#  *(( **(           %(//////****,,,,,.(*               .((*
#   /(# **(          (****#&@@@@&#*,,,,,,           @@/ ((*
#    *(( ,**          ,&//////////***%&           #@@%*((*
#     *((#.**/                                  (@@@.#((.
#       *((#.***(                             @@@& #((*
#         *(((,,                          @@@@@(((((*
#            *(((#      %@@@*     /@@@@@@@@# #(((*
#               */((((#, (##&@@@@@%##/ /#((((*,
#                    ,*/(((((((((((((((/*.

# Cloning with Python

myinsertdna = "gtctgtgttgtt"

myplasmid = "ccatttgtatgttcagctaaCCCGGGcttctacccatcccccgaag"

def gc(s):
    gc_fraction = 0  # your code here!
    return gc_fraction

myplasmidGC = gc(myplasmid)

SmaI = "CCCGGG"

myplasmid.find(SmaI)

myrecplasmid = "?"














import pydna

pydna.logo()

from pydna.dseq import Dseq

Dseq("aaa", "ttt")

Dseq("aaa", "ttt", ovhg=0)

Dseq("aaa", "ttt", ovhg=1)

Dseq("aaa", "ttt", ovhg=-1)

Dseq("aaa", "ttt", linear=False, ovhg=0)

Dseq("aaa", "ttt", circular=True, ovhg=0)

Dseq("ggaaa", "ccttt")



from pydna.dseqrecord import Dseqrecord

Dseqrecord("GGATCC")

from Bio.Restriction import BamHI

a, b = Dseqrecord("GGATCC").cut(BamHI)

a, b

a.seq

b.seq

(a+b).seq

(b+a).seq

from Bio.Restriction import EcoRI

Dseqrecord("GAATTC")



from pydna.genbank import genbank

gbfile = genbank("CS570233.1")

