#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 21 15:41:42 2021

@author: bjorn


https://en.wikipedia.org/wiki/Cre-Lox_recombination

13bp	      8bp	   13bp
ATAACTTCGTATA-NNNTANNN-TATACGAAGTTAT


Name	    13 bp  	        8 bp  	    13 bp 
            Recognition     Spacer      Recognition 
            Region          Region      Region

Wild-Type	ATAACTTCGTATA	ATGTATGC	TATACGAAGTTAT
lox 511	    ATAACTTCGTATA	ATGTATaC	TATACGAAGTTAT
lox 5171	ATAACTTCGTATA	ATGTgTaC	TATACGAAGTTAT
lox 2272	ATAACTTCGTATA	AaGTATcC	TATACGAAGTTAT
M2	        ATAACTTCGTATA	AgaaAcca	TATACGAAGTTAT
M3	        ATAACTTCGTATA	taaTACCA	TATACGAAGTTAT
M7	        ATAACTTCGTATA	AgaTAGAA	TATACGAAGTTAT
M11	        ATAACTTCGTATA	cgaTAcca	TATACGAAGTTAT
lox 71	    TACCGTTCGTATA	NNNTANNN	TATACGAAGTTAT
lox 66	    ATAACTTCGTATA	NNNTANNN	TATACGAACGGTA

"""