#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 14:24:21 2020

@author: tsanga

Called by main.py
"""
import re

# This function creates a comment with called mutations and extra undesired mutations ("xmutations")
def annotate(origComment, goodMutations, extraMutations, map2refvar, denovovar, denovoseq, scaffold):
    
    from reference_sequences import scaffold_protein
    from reference_sequences import scaffold_dna
    from sequence_funcs import codon_index
    
    mutRegex = re.compile(r'[A-Z]\d+[A-Z]|[A-Z]\d+\*') # search pattern for any mutation (e.g. N391X, T293Y, L405*)
    mutResidueRegex = re.compile(r'[A-Z]\d+([A-Z]|\*)') # search pattern to return the mutation residue, e.g. X, Y, or *
    aanumRegex = re.compile(r'[A-Z](\d+)') # search pattern for residue number
    
    newComment = '_DS: '
    
    scaffoldSeq = scaffold_dna[str(scaffold)]
    reflength = len(scaffoldSeq)
    
    
    if goodMutations is not None:
        numGoodMuts = len(goodMutations)
        for i in range(0, numGoodMuts):
            mutation = goodMutations[i] # gets mutation name (e.g. T392Y) as string
            mutResidue = mutResidueRegex.findall(mutation) # gets mutated residue (e.g. Y or *) as list
            residueNum = aanumRegex.findall(mutation) # finds the residue number (e.g. 392) as list
            residueIndex = int(residueNum[0])# convert to integer
            
            ### make an indexable list of codons from denovo sequence to get the mutation codons
            denovoseq_codon_index = codon_index(denovoseq) # produces a list of codons for denovoseq
            mutCodon = denovoseq_codon_index[residueIndex-1] # -1 because python starts counting at zero
            
            ### full or partial length
            
            newComment += mutation + ', ' + mutCodon + ', '
        
        # moved out of the above for loop  8/17/20
        if len(denovoseq) > 0.9 * reflength: # here, having more than 90% of the gene sequence length is considered full length
            seqlength = 'full-length'
            newComment += seqlength  + ', '
        
        else:
            seqlength = 'partial-length'
            
            # Perform a local alignment
            from Bio import pairwise2
            gap_open_penalty = -2
            gap_extend_penalty = -1
            partialAlignment = pairwise2.align.localxs(scaffoldSeq, denovoseq, gap_open_penalty, gap_extend_penalty)
            while len(partialAlignment) > 1: # run pairwise2 with increasing stringency while more than one alignment is returned
                if gap_open_penalty < -20: # something is clearly very wrong with the data at this point
                    break
                gap_open_penalty -= 1
                gap_extend_penalty -=0.5
                partialAlignment = pairwise2.align.localxs(scaffoldSeq, denovoseq, gap_open_penalty, gap_extend_penalty)
            
            if len(partialAlignment) > 0:
                refAligned = partialAlignment[0].seqA
                denovoAligned = partialAlignment[0].seqB
                start_base = partialAlignment[0].start
                end_base = partialAlignment[0].end
            else:
                start_base = 'N/A' # REVISE THIS CODE
                end_base = 'N/A'
                
            newComment += seqlength + ', ' + 'start-base ' + str(start_base) + ', ' + 'end-base ' + str(end_base) + ', '
            
        
    if extraMutations is not None:
        numExtraMuts = len(extraMutations)
        for j in range(0, numExtraMuts):
            xmutation = extraMutations[j] # gets mutation name (e.g. T392Y) as string
            xmutResidue = mutResidueRegex.findall(xmutation) # gets mutated residue (e.g. Y or *) as list
            xresidueNum = aanumRegex.findall(xmutation) # finds the residue number (e.g. 392) as list
            xresidueIndex = int(xresidueNum[0]) # convert to integer
            
            newComment += 'sub-' + xmutation + ', '
        
    finalComment = origComment + newComment
    return finalComment


