#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 14:24:21 2020

@author: tsanga

Called by main.py
"""
import re

def annotate(origComment, goodMutations, extraMutations, map2refvar, denovovar, denovoseq, scaffold):
    
    from reference_sequences import scaffold_protein
    from reference_sequences import scaffold_dna
    from sequence_funcs import codon_index
    
    mutRegex = re.compile(r'[A-Z]\d+[A-Z]|[A-Z]\d+\*') # search pattern for any mutation (e.g. N391X, T293Y, L405*)
    mutResidueRegex = re.compile(r'[A-Z]\d+([A-Z]|\*)') # search pattern to return the mutation residue, e.g. X, Y, or *
    aanumRegex = re.compile(r'[A-Z](\d+)') # search pattern for residue number
    
    newComment = ': '
    
    # as written, this will error if denovo seqeuence is not a full length gene sequence
    
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
            scaffoldSeq = scaffold_dna[str(scaffold)]
            reflength = len(scaffoldSeq)
            
            if len(denovoseq) > 0.9 * reflength:
                seqlength = 'full-length'
            else:
                seqlength = 'partial-length'
                
            if seqlength == 'partial-length':
                #look for the first 3 codons of partial length seq to find start-base
                # look for the last 3 codons of partial length seq to find end-base
                start_search = denovoseq[0:8]
                end_search = denovoseq[len(denovoseq)-9 : len(denovoseq)]
                
                start_base = scaffoldSeq.find(start_search)
                end_base = scaffoldSeq.find(end_search) + 9
            
            ### Build the comment
            if seqlength == 'full-length':
                newComment += mutation + ', ' + mutCodon + ', ' + seqlength  + ', '
            else:
                newComment += mutation + ', ' + mutCodon + ', ' + seqlength + ', ' \
                    + 'start-base ' + str(start_base) + ', ' + 'end-base ' + str(end_base)
        
    if extraMutations is not None:
        numExtraMuts = len(extraMutations)
        for j in range(0, numExtraMuts):
            xmutation = extraMutations[j] # gets mutation name (e.g. T392Y) as string
            xmutResidue = mutResidueRegex.findall(xmutation) # gets mutated residue (e.g. Y or *) as list
            xresidueNum = aanumRegex.findall(xmutation) # finds the residue number (e.g. 392) as list
            xresidueIndex = int(xresidueNum[0]) # convert to integer
            
            # xmutCodon = denovoseq_codon_index[xresidueIndex-1] # -1 because python starts counting at zero
            
            #if len(denovoseq) == len(scaffoldSeq): # a simple codon change, seq length is the same
                # 
            # determine amino acid substitutions
            # determine frameshifts
            # determine snps
            # determine indels
            
            newComment += 'sub-' + xmutation + ', '
            
            # find index of extra mutation
            
            # get original residue by xresidueNuym
            
    # Check for a frameshift in denovovar
    from sequence_funcs import findFrameshift
    frameshiftList = findFrameshift(denovovar)
    
    if frameshiftList is not None:
        if frameshiftList[0] == True:
            newComment += 'frameshift-' + frameshiftList[1] + ', '
        
    finalComment = origComment + newComment
    return finalComment

    
    #Search for an existing comment (there shouldn't be one if no Sanger sequence exists)
    prevComment = origComment.find(':')
    if prevComment == -1: # condition when no comment is found
        print('')
    
    # get the mutation(s) found
    # get the codon(s) found
    # get full length or partial length
    # list additional mutations
    # start-base number if applicable
    # end-base number if partial length