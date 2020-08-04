#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 17:58:52 2020

@author: tsanga

Called by main.py, updateMutants.py, and annotations.py
"""
from cdb_globals import *

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



# if comparing map2ref and denovoseq, pass denovoseq as sequence 1
def simple_alignment(seq1, seq2):
    from Bio import pairwise2
    gap_open_penalty = -2
    gap_extend_penalty = -1
    alignment = pairwise2.align.globalxs(seq1, seq2, gap_open_penalty, gap_extend_penalty)
    while len(alignment) > 1:
        if gap_open_penalty < -10: # something is clearly very wrong with the data at this point
                break
        gap_open_penalty -= 1
        gap_extend_penalty -=0.5
        alignment = pairwise2.align.globalxs(seq1, seq2, gap_open_penalty, gap_extend_penalty)
        
        refseq = alignment[0].seqA # calling each named tuple item returned by pairwise2
        testseq = alignment[0].seqB
        
        alignment_comment = ''
        
        for j in range(len(refseq)):
    
            if refseq[j] != testseq[j]:
                if testseq[j] == '-':
                    alignment_comment += str(refseq[j]).lower() + str(j+1) + 'del, '
                elif refseq[j] == '-': # and testseq[j+1] == refseq[j]:
                        alignment_comment += ' ins' + str(j+1) + ', ' #str(refseq[j+1]).lower() + str(j+2) + ', '
                else:
                    alignment_comment += 'snp-' + str(refseq[j]).lower() + str(j+1) + str(testseq[j]).lower() + ', '
            
    print(alignment_comment)
    return alignment_comment



### find single nucleotide polymorphisms
# exclude the intended SDM codons
def polymorphisms(scaffold, denovoseq, goodMutations):
    import re

    aanumRegex = re.compile(r'[A-Z](\d+)') # search pattern for residue number
    
    if goodMutations is not None:
        
        from reference_sequences import scaffold_dna
        ref_dna = scaffold_dna[str(scaffold)] # get the reference dna sequence
        ref_dna_codons = codon_index(ref_dna) # break the sequence into codons
        
        numGoodMuts = len(goodMutations)
        sdm_bases = []
        for k in range (0, numGoodMuts):
            sdm_sites = aanumRegex.findall(goodMutations[k]) # find the sdm residue numbers to evaluate later
                
        # since this is looking for miscellaneous dna polymorphisms, 
        # make a list of sdm base numbers to exclude from analysis
            if sdm_sites is not None:
                for i in range(len(sdm_sites)):
                    codonbasenum = 3*int(sdm_sites[i])
                    sdm_bases.append(codonbasenum-2)
                    sdm_bases.append(codonbasenum-1)
                    sdm_bases.append(codonbasenum)
                
            # print(sdm_bases)
        
            # trim scaffold sequence at 3' end if necessary (scaffold probably has the stop codon, while denovoseq omits it)
            if len(ref_dna) > len(denovoseq):
                while len(ref_dna) > len(denovoseq):
                    ref_dna = ref_dna[:len(ref_dna)-1]
                    
            
            # create the alignment
            from Bio import pairwise2
            gap_open_penalty = -2
            gap_extend_penalty = -1
            alignment = pairwise2.align.globalxs(ref_dna, denovoseq, gap_open_penalty, gap_extend_penalty)
            while len(alignment) > 1:
                if gap_open_penalty < -10: # something is clearly very wrong with the data at this point
                        break
                gap_open_penalty -= 1
                gap_extend_penalty -=0.5
                alignment = pairwise2.align.globalxs(ref_dna, denovoseq, gap_open_penalty, gap_extend_penalty)
            
            refseq = alignment[0].seqA # calling each named tuple item returned by pairwise2
            testseq = alignment[0].seqB
            
            polymorphism_comment = ''
            
            # comment the mutations
            # due to ubiquity, some mutations are silenced. See reference_sequences
            
            from reference_sequences import suppressed_snps
            # from cdb_globals import dna_mutlist
            for j in range(len(refseq)):
                if j+1 not in sdm_bases:
                    if refseq[j] != testseq[j]:
                        if testseq[j] == '-':
                            mut_comment = str(refseq[j]).lower() + str(j+1) + 'del, '
                            if mut_comment not in suppressed_snps:
                                polymorphism_comment += mut_comment
                            dna_mutlist.append(mut_comment)
                        elif refseq[j] == '-': # and testseq[j+1] == refseq[j]:
                            mut_comment = ' ins' + str(j+1) + ', '
                            if mut_comment not in suppressed_snps:
                                polymorphism_comment +=  mut_comment #str(refseq[j+1]).lower() + str(j+2) + ', '
                            dna_mutlist.append(mut_comment)
                        else:
                            mut_comment = 'snp-' + str(refseq[j]).lower() + str(j+1) + str(testseq[j]).lower() + ', '
                            if mut_comment not in suppressed_snps:
                                polymorphism_comment += mut_comment
                            dna_mutlist.append(mut_comment)
                    
            # print(polymorphism_comment)
                    
            # return alignment, polymorphism_comment
            return polymorphism_comment