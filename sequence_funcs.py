#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 17:58:52 2020

@author: tsanga

Called by main.py, updateMutants.py, and annotations.py
"""


def checkSequence(seqtext):
    seqlen = len(seqtext)
    if seqlen > 1:
        seqfound = True
    else:
        seqfound = False
    #print(seqlen)
    return seqfound  # Boolean out


# check if assembled sequences are the same. Sequences already converted 
# to upper case by main.py
def compareSeq(map2refseq, denovoseq):
    map2refseq = map2refseq.upper()
    denovoseq = denovoseq.upper()
    if map2refseq == denovoseq:
        return True
    else: # if there is a common mutation but the sequences disagree
        return False # more code to be written
    

def codon_index(dnaseq):
    codonlist=[]
    if len(dnaseq)%3 == 0:
        length = len(dnaseq)
    elif len(dnaseq)%3 == 1:
        length = len(dnaseq)-3+len(dnaseq)%3
    else:
        length = len(dnaseq)-2
        
    for i in range(0, length, 3):
        codonlist.append(dnaseq[i:i+3])
    
    return codonlist
    

# If commonvar exists but the assembled sequences are not identical, 
def seqMismatch(map2refseq, denovoseq):
    print('find the discrepancy')
    

def findFrameshift(mutations): #pass a mutation list into this function
    # Frameshifts should appear in the denovo_var field
    if mutations is not None:
        import re
        mutRegex = re.compile(r'[A-Z]\d+[A-Z]|[A-Z]\d+\*')
        aanumRegex = re.compile(r'[A-Z](\d+)') # search pattern for residue number
        # Parse for residue numbers, e.g. 399, 400, 402, 403, 404... 
        mutNames = mutRegex.findall(mutations)
        residues = aanumRegex.findall(mutations)
        intervalList = []
        frameshift = False
        # make a list of intervals between residue numbers, e.g. [1, 2, 1, 1....]
        for i in range(0, len(residues)-1):
            intervalList.append(int(residues[i+1])-int(residues[i]))
        
        if len(intervalList) > 3: # search only mutation lists longer than 3 substitutions for frameshifts. Arbitrary threshold.
            for j in range(0, len(intervalList)):
                if intervalList[j] < 3: # if the amino acids are less than 3 positions apart
                    frameshift = True
                else:
                    frameshift = False
                    break
        
        # return frameshift as Boolean, and the first found mutation of the frameshift. 
        if frameshift == True:
            return frameshift, mutNames[0]
    
    
### DNA translation dictionary and function
# input a dna sequence
# output an amino acid sequence

# borrowed from:
# https://towardsdatascience.com/starting-off-in-bioinformatics-turning-dna-sequences-into-protein-sequences-c771dc20b89f


def translate(dna):
    from reference_sequences import protein
    dna = dna.upper()
    prot_seq = ''
    
    for i in range(0, len(dna)-3+len(dna)%3, 3):
        codon = dna[i:i+3]
        prot_seq += protein[dna[i:i+3]]
        if protein[codon] == "*":
            break
            
    return prot_seq