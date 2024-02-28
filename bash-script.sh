#!/usr/bin/env bash

rrep --recursive --fixed-strings "ldseguid=5v2YtLDJRhA0_O_Ut6HUXr_4EIc" "ldseguid=2oLvuMsAc-fgAISC7b3XJBuejGI" . --include=*.py
rrep --recursive --fixed-strings "cdseguid=_nYOPeDWBC7OaB2Arnt5x5fAFZM" "cdseguid=fLRW6jTMqqwZCbIFOGyl1TdDbrs" . --include=*.py
rrep --recursive --fixed-strings "cdseguid=mIbMOKcoiRlnUs7ZaYwBxQnFego" "cdseguid=cJQ-aN4c5sLYq7FcNLKr5ZVSdSM" . --include=*.py
rrep --recursive --fixed-strings "" "" . --include=*.py
rrep --recursive --fixed-strings "" "" . --include=*.py
rrep --recursive --fixed-strings "" "" . --include=*.py
rrep --recursive --fixed-strings "" "" . --include=*.py
rrep --recursive --fixed-strings "" "" . --include=*.py
