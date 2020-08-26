#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 17:58:52 2020

@author: tsanga

Called by main.py, updateMutants.py, and annotations.py
"""
from cdb_globals import *
import re

def checkSequence(seqtext):
    if seqtext is not None:
        seqlen = len(seqtext)
    else:
        seqlen = 1
        
    if seqlen > 1:
        seqfound = True
    else:
        seqfound = False
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
    

def compareSanger(constructname, commonvar, map2refseq, denovoseq, scaffold, searchtext):
    mutRegex = re.compile(r'[A-Z]\d+[A-Z]') # original mutation regex
    # mutRegex = re.compile(r'[A-Z]\d+[A-Z]')
    stopregex = re.compile(r'[A-Z]\d+STOP')
    asteriskregex = re.compile(r'[A-Z]\d+\*')
    origRegex = re.compile(r'([A-Z]\d+)[A-Z]') # original residue regex
    indexRegex = re.compile(r'[A-Z](\d+)') # residue number regex
    
    constructname_mutations = []
    
    # pass constructname through regex to find mutations
    if 'STOP' in constructname:
        constructname_stops = stopregex.findall(constructname) #mutRegex will find 'T559STOP' as 'T559S', so need to clear out stops first.
        for h in range(len(constructname_stops)): # delete stop mutations
            temp_constructname = constructname.replace(constructname_stops[h], '')
        
        # find all non-stop mutations
        constructname_mutations = mutRegex.findall(temp_constructname) # find all mutations except stops
        
        for g in range(len(constructname_stops)):
            constructname_mutations.append(constructname_stops[h])
    else:
        constructname_mutations = mutRegex.findall(constructname)
        
    commonvarMut = []
    # if commonvar contains asterisk, convert to 'STOP'
    if commonvar is not None:
        findasterisk = commonvar.find('*')
        if findasterisk > -1:
            commonvar = commonvar.replace('*', 'STOP')
            commonvarMut_stop = stopregex.findall(commonvar)
            for f in range(len(commonvarMut_stop)):
                temp_commonvar = commonvar.replace(commonvarMut_stop[f], '')
            
            commonvarMut = mutRegex.findall(temp_commonvar)
            
            for e in range(len(commonvarMut_stop)):
                commonvarMut.append(commonvarMut_stop[e])
        else:
            commonvarMut = mutRegex.findall(commonvar)
    else: commonvarMut = []
    
    # concatenate constructname mutations into same format as commonvar
    sangervar = ''
    if len(constructname_mutations) > 0:
        for i in range(len(constructname_mutations)):
            # read constructname for lowercase stops, convert to upper
            findlowerstop = constructname_mutations[i].find('stop')
            # print(findlowerstop)
            if findlowerstop > -1:
                constructname_mutations[i] = constructname_mutations[i].upper()
                
            if i < len(constructname_mutations)-1:
                sangervar += constructname_mutations[i] + ','
            else:
                sangervar += constructname_mutations[i]
                
        # # get map2refFound and denovoFound (the found mutations from sequence)
        # map2refTranslation = translate(map2refseq) # translate the dna sequences provided
        # denovoTranslation = translate(denovoseq)
        
        # # get the position(s) of interest from constructname
        # positions = indexRegex.findall(constructname)
        # # get the original residue (e.g. A123) from constructname
        # origResidue = origRegex.findall(constructname)
        
        # map2refFound = [] # stores individual mutation residue letters (e.g. A, Y, T)
        # denovoFound = []
        # map2refFoundMut = [] # stores the concatenated new mutation  (e.g. L330A)
        # denovoFoundMut = []
    map2refString = ''
    denovoString = ''
        
    foundMutations = findMutations(scaffold, map2refseq, denovoseq, searchtext)
        
    map2refMut = foundMutations[0]
    denovoMut = foundMutations[1]
        
    # create spreadsheet strings
    for i in range(len(map2refMut)):
        if i < len(map2refMut)-1:
            map2refString += map2refMut[i] + ','
        else:
            map2refString += map2refMut[i]
        
    for j in range(len(denovoMut)):
        if j < len(denovoMut)-1:
            denovoString += denovoMut[j] + ','
        else:
            denovoString += denovoMut[j]
        
        # if len(positions) > 0:
        # # for each position id, find the mutation from map2refTranslation and denovoTranslation
        #     for i in range(0, len(positions)):
        #         aaindex = int(positions[i])
        #         map2refFound.append(map2refTranslation[aaindex-1]) # stores the amino acid residue letter
        #         denovoFound.append(denovoTranslation[aaindex-1])
            
        #     for j in range(0, len(map2refFound)):
        #         map2refNewMut = origResidue[j] + str(map2refFound[j])
        #         map2refFoundMut.append(map2refNewMut)
                
        #         # create string for spreadsheet
        #         if j < len(map2refFound) -1:
        #             map2refFoundString += origResidue[j] + str(map2refFound[j]) + ','
        #         else:
        #             map2refFoundString += origResidue[j] + str(map2refFound[j])
                
        #     for k in range(0, len(denovoFound)):
        #         denovoNewMut = origResidue[j] + str(map2refFound[j])
        #         denovoFoundMut.append(map2refNewMut)
        #         if k < len(denovoFound) - 1:
        #             denovoFoundString += origResidue[k] + str(denovoFound[k]) + ','
        #         else:
        #             denovoFoundString += origResidue[k] + str(denovoFound[k])
                    
    if commonvar is None and map2refMut == denovoMut:
        commonvarMut = denovoMut
    
    print(commonvarMut)
    print(denovoMut)
        
    
    if commonvar is not None:
        
        if len(constructname_mutations) > 1:
            for d in range(len(constructname_mutations)):
                if constructname_mutations[d] in map2refMut or constructname_mutations[d] in denovoMut:
                    print("condition 1a - combo")
                    return sangervar, True, map2refString, denovoString
                else: 
                    return sangervar, False, map2refString, denovoString
        elif sangervar == map2refString or sangervar == denovoString or sangervar == commonvar:
            print("condition 1b - combo")
            return sangervar, True, map2refString, denovoString
        
        # these evaluations are good for single mutations
        if sangervar in commonvarMut:
            print('condition 1')
            return sangervar, True, map2refString, denovoString
        else:
            print('condition 2')
            return sangervar, False, map2refString, denovoString
        
    # where commonvar is None
    elif sangervar in map2refMut or sangervar in denovoMut:
        print('condition 3')
        return sangervar, True, map2refString, denovoString
    
    elif len(constructname_mutations) > 1:
        
        for d in range(len(constructname_mutations)):
            if constructname_mutations[d] in map2refMut or constructname_mutations[d] in denovoMut:
                print('condition 4')
                return sangervar, True, map2refString, denovoString
            else: 
                return sangervar, False, map2refString, denovoString
        
        
        
    else:
        print('condition 5')
        return sangervar, False, map2refString, denovoString
        
    
            
    print('COMPARE SANGER')
    print(constructname_mutations)
    print(map2refMut)
    print(denovoMut)
    print(commonvar)
    print(sangervar)
        

    
    


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



def findMutations(scaffold, map2refseq, denovoseq, searchtext):
    #
    # get scaffold sequence from scaffold number input
    from reference_sequences import scaffold_protein
    scaffold_prot_seq = scaffold_protein[scaffold]
    
    # translate the sequences
    map2refTranslation = translate(map2refseq)
    denovoTranslation = translate(denovoseq)
    
    # align the map2ref sequences
    from Bio import pairwise2
    gap_open_penalty = -1
    gap_extend_penalty = -0.5
    # full length alignment
    if len(map2refTranslation) > 0.95 * len(scaffold_prot_seq):
        map2refAlignment = pairwise2.align.globalxs(scaffold_prot_seq, map2refTranslation, gap_open_penalty, gap_extend_penalty)
        while len(map2refAlignment) > 1 or gap_open_penalty > -10: 
            gap_open_penalty -= 1
            gap_extend_penalty -= 1
            map2refAlignment = pairwise2.align.globalxs(scaffold_prot_seq, map2refTranslation, gap_open_penalty, gap_extend_penalty)
            print(searchtext + ' Global: Trying gap penalty ' + str(gap_open_penalty))
            
            if gap_open_penalty <= -10: # break and go to local alignment
                gap_penalty_break = 1
                break
        
        if gap_penalty_break == 1:
            map2refAlignment = pairwise2.align.localxs(scaffold_prot_seq, map2refTranslation, gap_open_penalty, gap_extend_penalty)
            while len(map2refAlignment) > 1:
                gap_open_penalty -= 1
                gap_extend_penalty -= 1
                map2refAlignment = pairwise2.align.localxs(scaffold_prot_seq, map2refTranslation, gap_open_penalty, gap_extend_penalty)
                print(searchtext + ' Local: Trying gap penalty ' + str(gap_open_penalty))
                if gap_open_penalty < -20:
                    return None
                    break
            
    
    else: # partial length alignment
        map2refAlignment = pairwise2.align.localxs(scaffold_prot_seq, map2refTranslation, gap_open_penalty, gap_extend_penalty)
        while len(map2refAlignment) > 1:
            gap_open_penalty -= 1
            gap_extend_penalty -= 1
            map2refAlignment = pairwise2.align.localxs(scaffold_prot_seq, map2refTranslation, gap_open_penalty, gap_extend_penalty)
            print(searchtext + ' Local: Trying gap penalty ' + str(gap_open_penalty))
            if gap_open_penalty < -20:
                    return None
                    break
            
        
    # align the denovo sequences
    from Bio import pairwise2
    gap_open_penalty = -1
    gap_extend_penalty = -0.5
    # full length alignment
    if len(denovoTranslation) > 0.95 * len(scaffold_prot_seq):
        denovoAlignment = pairwise2.align.globalxs(scaffold_prot_seq, denovoTranslation, gap_open_penalty, gap_extend_penalty)
        while len(denovoAlignment) > 1 or gap_open_penalty > -10:
            gap_open_penalty -= 1
            gap_extend_penalty -= 1
            denovoAlignment = pairwise2.align.globalxs(scaffold_prot_seq, denovoTranslation, gap_open_penalty, gap_extend_penalty)
            print(searchtext + ' Global: Trying gap penalty ' + str(gap_open_penalty))
            
            if gap_open_penalty <= -10: # break and go to local alignment
                gap_penalty_break = 1
                break
            
        if gap_penalty_break == 1:
            denovoAlignment = pairwise2.align.localxs(scaffold_prot_seq, denovoTranslation, gap_open_penalty, gap_extend_penalty)
            while len(denovoAlignment) > 1:
                gap_open_penalty -= 1
                gap_extend_penalty -= 1
                denovoAlignment = pairwise2.align.localxs(scaffold_prot_seq, denovoTranslation, gap_open_penalty, gap_extend_penalty)
                print(searchtext + ' Local: Trying gap penalty ' + str(gap_open_penalty))
                if gap_open_penalty < -20:
                    return None
                    break
    
    else: # partial length alignment
        denovoAlignment = pairwise2.align.localxs(scaffold_prot_seq, denovoTranslation, gap_open_penalty, gap_extend_penalty)
        while len(denovoAlignment) > 1:
            gap_open_penalty -= 1
            gap_extend_penalty -= 1
            denovoAlignment = pairwise2.align.localxs(scaffold_prot_seq, denovoTranslation, gap_open_penalty, gap_extend_penalty)
            print(searchtext + ' Local: Trying gap penalty ' + str(gap_open_penalty))
            if gap_open_penalty < -20:
                return None
                break
    
    map2refseqA = map2refAlignment[0].seqA
    map2refseqB = map2refAlignment[0].seqB
    
    denovoseqA = denovoAlignment[0].seqA
    denovoseqB = denovoAlignment[0].seqB
    
    map2refMut = []
    denovoMut = []
    
    for i in range(len(map2refTranslation)):
        if map2refseqA[i] != map2refseqB[i]:
            if map2refseqB[i] == '*':
                mutation = map2refseqA[i] + str(i+1) + 'STOP'
            else:
                mutation = map2refseqA[i] + str(i+1) + map2refseqB[i]
            map2refMut.append(mutation)
            
    for i in range(len(denovoTranslation)):
        if denovoseqA[i] != denovoseqB[i]:
            if denovoseqB[i] == '*':
                mutation = denovoseqA[i] + str(i+1) + 'STOP'
            else:
                mutation = denovoseqA[i] + str(i+1) + denovoseqB[i]
            denovoMut.append(mutation)
            
    return map2refMut, denovoMut


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
        
        frameshiftBreak = 0 # gets incremented when a reside meets frameshift criteria. 
        
        frameshiftThreshold = 5 # when frameshiftBreak gets to this arbitrary threshold, 
        # then this "trips" the frameshift condition True no matter what. This prevents other spaced 
        # out mutations from nullifying the frameshift condition.
        
        frameshiftTrueList = [] # tracks whether frameshifts have been found yet. 
        
        if len(intervalList) > 3: # search only mutation lists longer than 3 substitutions for frameshifts. Arbitrary threshold.
            for j in range(0, len(intervalList)):
                if intervalList[j] < 3: # if the amino acids are less than 3 positions apart
                    if len(frameshiftTrueList)  == 0: 
                            firstFrameshiftSite = mutNames[j] # store the first site of the frameshift.
                    frameshift = True
                    frameshiftBreak += 1
                    frameshiftTrueList.append(1)
                
                elif frameshiftBreak >= 5:
                    frameshift = True
                
                else:
                    frameshift = False
                    frameshiftBreak -= 1
                    
        
        # return frameshift as Boolean, and the first found mutation of the frameshift. 
        if frameshiftBreak >= 5 or frameshift == True:
            return frameshift, firstFrameshiftSite
    
    
### DNA translation dictionary and function
# input a dna sequence
# output an amino acid sequence

# borrowed from:
# https://towardsdatascience.com/starting-off-in-bioinformatics-turning-dna-sequences-into-protein-sequences-c771dc20b89f


def iterativeAlignment(seq, refseq):
    # align the denovo sequences
    from Bio import pairwise2
    gap_open_penalty = -1
    gap_extend_penalty = -0.5
    # full length alignment
    if len(seq) > 0.95 * len(refseq):
        denovoAlignment = pairwise2.align.globalxs(refseq, seq, gap_open_penalty, gap_extend_penalty)
        while len(denovoAlignment) > 1 or gap_open_penalty > -10:
            gap_open_penalty -= 1
            gap_extend_penalty -= 1
            denovoAlignment = pairwise2.align.globalxs(refseq, seq, gap_open_penalty, gap_extend_penalty)
            
            if gap_open_penalty <= -10: # break and go to local alignment
                gap_penalty_break = 1
                break
            
        if gap_penalty_break == 1:
            denovoAlignment = pairwise2.align.localxs(refseq, seq, gap_open_penalty, gap_extend_penalty)
            while len(denovoAlignment) > 1:
                gap_open_penalty -= 1
                gap_extend_penalty -= 1
                denovoAlignment = pairwise2.align.localxs(refseq, seq, gap_open_penalty, gap_extend_penalty)
    
    else: # partial length alignment
        denovoAlignment = pairwise2.align.localxs(refseq, seq, gap_open_penalty, gap_extend_penalty)
        while len(denovoAlignment) > 1:
            gap_open_penalty -= 1
            gap_extend_penalty -= 1
            denovoAlignment = pairwise2.align.localxs(refseq, seq, gap_open_penalty, gap_extend_penalty)
    
    denovoTransAlign = denovoAlignment[0].seqA
    refTransAlign = denovoAlignment[0].seqB
    return denovoTransAlign, refTransAlign
    

def translate(dna):
    from reference_sequences import protein
    dna = dna.upper()
    prot_seq = ''
    
    for i in range(0, len(dna)-3+len(dna)%3, 3):
        codon = dna[i:i+3]
        if codon.find('N') > -1: # can occur where nucelotide sequence is missing bases, 'n' gets filled in. 
            prot_seq += 'x'
        else:
            prot_seq += protein[dna[i:i+3]]
            
            if protein[codon] == "*" or i > len(dna)-3:
                break
            
        
            
    return prot_seq



# if comparing map2ref and denovoseq, pass denovoseq as sequence 1
def simple_alignment(denovoseq, map2refseq):
    from Bio import pairwise2
    gap_open_penalty = -2
    gap_extend_penalty = -1
    alignment = pairwise2.align.globalxs(denovoseq, map2refseq, gap_open_penalty, gap_extend_penalty)
    while len(alignment) > 1:
        if gap_open_penalty < -20: # something is clearly very wrong with the data at this point
                break
        gap_open_penalty -= 1
        gap_extend_penalty -=0.5
        alignment = pairwise2.align.globalxs(denovoseq, map2refseq, gap_open_penalty, gap_extend_penalty)
        
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
        # ref_dna_codons = codon_index(ref_dna) # break the sequence into codons
        
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
            # if len(ref_dna) > len(denovoseq):
            #     while len(ref_dna) > len(denovoseq):
            #         ref_dna = ref_dna[:len(ref_dna)-1] 
                    
            
            # create the alignment
            from Bio import pairwise2
            gap_open_penalty = -2
            gap_extend_penalty = -1
            alignment = pairwise2.align.globalxs(ref_dna, denovoseq, gap_open_penalty, gap_extend_penalty)
            while len(alignment) > 1: # run pairwise2 with increasing stringency while more than one alignment is returned
                if gap_open_penalty < -20: # something is clearly very wrong with the data at this point
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
            for j in range(len(refseq)-6): # denovo assembly gives variable length toward the end, ignore these bases. 
                if j+1 not in sdm_bases:
                    if refseq[j] != testseq[j]:
                        if testseq[j] == '-':
                            mut_comment = str(refseq[j]).lower() + str(j+1) + 'del, '
                            if mut_comment not in suppressed_snps[scaffold]:
                                polymorphism_comment += mut_comment
                            dna_mutlist.append(mut_comment) # tracking all mutations
                        elif refseq[j] == '-': # and testseq[j+1] == refseq[j]:
                            mut_comment = ' ins' + str(j+1) + ', '
                            if mut_comment not in suppressed_snps[scaffold]:
                                polymorphism_comment +=  mut_comment #str(refseq[j+1]).lower() + str(j+2) + ', '
                            dna_mutlist.append(mut_comment)
                        else:
                            mut_comment = 'snp-' + str(refseq[j]).lower() + str(j+1) + str(testseq[j]).lower() + ', '
                            if mut_comment not in suppressed_snps[scaffold]:
                                polymorphism_comment += mut_comment
                            dna_mutlist.append(mut_comment)
                    
            # check polymorphism comment for a frameshift like condition and limit it to first site
            
                    
            # return alignment, polymorphism_comment
            return polymorphism_comment