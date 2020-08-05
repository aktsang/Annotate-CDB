#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 19:52:42 2020

@author: tsanga


Frameshifts: call the first frameshift mutation, not just the first residue in the list.
Frameshifts: negate snp calls for frameshifted constructs which have them.
Frameshifts: many called frameshifts seem to have equal ins/del events, which should cancel out. 

Polymorphisms: don't call mutations that occur in the stop codon (i.e. for 376 this is bases 1282-1284'
"""

