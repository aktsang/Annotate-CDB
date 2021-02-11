#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v. 20200717_0039
Created on Mon Jul 13 19:15:40 2020

@author: tsanga

Called by main.py
"""


import re


def updateMutations(constructname, map2refvar, map2refseq, denovovar, denovoseq, commonvar, scaffold):
    
    # # mutRegex = re.compile(r'[A-Z]\d+[A-Z]|[A-Z]\d+\*') # search pattern for any mutation (e.g. N391X, T293Y, L405*)
    mutstarRegex = re.compile(r'[A-Z]\d+[A-Z]|[A-Z]\d+\*') # search pattern for any mutation (e.g. N391X, T293Y, L405*)
    aaRegex = re.compile(r'[A-Z]\d+') # search pattern for residue number
    resRegex = re.compile(r'([A-Z])\d+') # search pattern for first residue letter only, when used with .findall
    aanumRegex = re.compile(r'[A-Z](\d+)')
    stopstarRegex = re.compile(r'\*') # search for asterisk stop notation
    newMutRegex = re.compile(r'[A-Z]\d+([A-Z])|[A-Z]\d+(\*)')

    
# Look for a Sanger sequence and mutation in CDB. If it exists, do nothing to the clone name
    print('UPDATE MUTATIONS')    
    
    newname2 = '' # the new construct name with sequenced mutations
    mutList = [] # a list of the updated (desired) mutations
    miscMutList = [] # a list of unexpected mutations
    
    
    if scaffold ==  '514':
        from sequence_funcs import igabasnfr_translation
        cvar_copy = commonvar #create another copy of commonvar to parse in parallel with wild type searching later
        commonvar = igabasnfr_translation(commonvar)
        
        # make a list of mutations in cvar_copy
        cvar_muts = mutstarRegex.findall(cvar_copy)
            

    if commonvar is not None:

        cname_match = mutstarRegex.findall(constructname) # get list of mutations in CDB
        commonvar_match = mutstarRegex.findall(commonvar) # get list of mutations called by DeepSeq
        
        # convert '*' characters to 'STOP'
        for q in range(len(commonvar_match)):
            find_asterisk = stopstarRegex.search(commonvar_match[q])
            if find_asterisk is not None:
                if find_asterisk.group() == '*':
                    commonvar_match[q] = commonvar_match[q].replace('*', 'STOP')
                    
        print(cname_match)
        print(commonvar_match)
        
        num_ds_Matches = len(commonvar_match) # get number of mutations found
        num_cdb_Matches = len(cname_match) # get number of mutations expected
        
        # string splitting on constructname
        if len(cname_match) > 0: # found mutation sites in constructname
            firstMut = constructname.find(cname_match[0])
            lastMut = constructname.find(cname_match[num_cdb_Matches-1])
            lastMutLength = len(cname_match[num_cdb_Matches-1])
            lastMutChar = lastMut + lastMutLength
            name1 = constructname[:firstMut]
            name2 = constructname[firstMut:lastMutChar]
            name3 = constructname[lastMutChar:]
            
            newname2 = name2
        
            ### Compare each residue from deep seq to each residue in cdb
            # Deep seq may find multiple mutations. Is any of them the one we expect? 
            
            cname_aa_match = aaRegex.findall(name2) # find residue numbers from name2 fragment of constructname
            commonvar_aa_match = aaRegex.findall(commonvar) # find residue numbers from commonvar
            
            matched_cname_mut = [] # to store matched constructname mutations
            unmatched_cname_mut = [] # store unmatched constructname mutations, everything can be assumed wild type
            

          
            # put all mutation sites from the database constructname into unmatched_cname_mut[]
            # these will be moved if they get matched with mutations found by sequencing.
            for g in range(len(cname_match)): #cname_aa_match
                unmatched_cname_mut.append(cname_match[g]) #cname_aa_match[g]
                print(unmatched_cname_mut)

            for r in range(0, num_ds_Matches): # for each deep seq identified mutation...
                for s in range(0, len(cname_match)): # cname_aa_match ...for each construct db intended mutation...
                    # print('r: ' + str(r) + ', ' + commonvar_match[r])
                    # print('s: ' + str(s) + ', ' + cname_match[s])
                    if commonvar_aa_match[r] == cname_aa_match[s]: # if the residue number (e.g. L405) is the same...
                        # newname2 += commonvar_match[r] + ' ' # add the deep seq identified mutation to the name.
                        newname2 = newname2.replace(cname_match[s], commonvar_match[r])
                        # print(newname2)
                        mutList.append(commonvar_match[r]) # store all the deep seq identified mutations in mutList
                        
                        matched_cname_mut.append(cname_match[s]) #cname_aa_match
                        
                        if cname_match[s] in unmatched_cname_mut: #cname_aa_match
                            unmatched_cname_mut.remove(cname_match[s]) #cname_aa_match
                            
                            
            if len(unmatched_cname_mut) > 0:
                
                # translate the sequences to look for wild type
                # check for wild type identity
                from reference_sequences import scaffold_protein
                scaffoldTranslated = scaffold_protein[scaffold]
                from sequence_funcs import translate
                denovoTranslated = translate(denovoseq)
                from sequence_funcs import iterativeAlignment
                denovoAlignment = iterativeAlignment(denovoTranslated, scaffoldTranslated)
                
                denovoAlignedSeq = denovoAlignment[0]
                scaffoldAlignedSeq = denovoAlignment[1]
                
                
                for f in range(len(unmatched_cname_mut)): # for each unmatched mutation from constructname (e.g. I29)
                    origResidueId = resRegex.findall(unmatched_cname_mut[f]) # e.g. ['I']
                    origResidueIndex = aanumRegex.findall(unmatched_cname_mut[f]) # e.g. ['29']
                    # print(origResidueId)
                    # print(origResidueIndex)
                    
                    if scaffold == '514':
                        # residues in unmatched_cname_mut are currently in jonny's numbering.
                        # convert to linear for indexing with the experimental sequence. 
                        from reference_sequences import igabasnfr_reverse_mutated
                        jonnyMutId = origResidueId[0] + origResidueIndex[0] # concatenate site ('I29')
                        linearMut = igabasnfr_reverse_mutated.get(jonnyMutId) # run reverse lookup
                        
                        # if a match is found, parse and use the linear index number for wild type searching.
                        if linearMut is not None:
                            linResidueIndex = aanumRegex.findall(linearMut)
                        else:
                            linResidueIndex = None
                    
                    # checking for wild type identity
                    if scaffold == '514':
                        if linResidueIndex is not None:
                            myindex = int(linResidueIndex[0])-1 # use the linear index number 
                        else:
                            myindex = int(origResidueIndex[0])-1
                    else:
                        myindex = int(origResidueIndex[0])-1 # for non-iGABA clones, this is the default. 
                        
                    if denovoAlignedSeq[myindex] == scaffoldAlignedSeq[myindex]:
                        if scaffold == '514':
                            # this is required to feed jonny's numbering system into annotations.py
                            wildtypemut = origResidueId[0] + origResidueIndex[0] + origResidueId[0]
                        else:
                            wildtypemut = origResidueId[0] + str(myindex+1) + origResidueId[0]
                        sdm_site = unmatched_cname_mut[f]
                        newname2 = newname2.replace(sdm_site, wildtypemut)
                        mutList.append(wildtypemut)
                                # print(newname2)
                                # print(denovoAlignedSeq[myindex])
                                # print(scaffoldAlignedSeq[myindex])
                    
                    
                        # if r < len(commonvar_aa_match)-1:
                        #     # newname2 += commonvar_match[r] + ' ' # add the deep seq identified mutation to the name.
                        #     newname2 = newname2.replace(cname_match[s], commonvar_match[r])
                        #     print(newname2)
                        #     mutList.append(commonvar_match[r]) # store all the deep seq identified mutations in mutList
                        # else: 
                        #     # newname2 += commonvar_match[r]
                        #     newname2 = newname2.replace(cname_match[s], commonvar_match[r])
                        #     print(newname2)
                        #     mutList.append(commonvar_match[r])
                                        
            # 2nd loop to check against mutList, everything else goes to miscMutList
            for t in range(0, num_ds_Matches):
                if commonvar_match[t] not in mutList:
                    # mutList.append(commonvar_match[t]) # just throw it in mutList now so annotations() finds the codon
                    miscMutList.append(commonvar_match[t]) # store all the deep seq identified mutations in miscMutList
                    
            
            # put misc mutants into the new construct name
            # if len(miscMutList) > 0:
            #     newname2 += ' '
            #     for u in range(len(miscMutList)):
            #         if u < len(miscMutList)-1:
            #             newname2 += miscMutList[u] + ' '
            #         else:
            #             newname2 += miscMutList[u]
                        
            if len(newname2) > 1:
                newconstructname = name1 + newname2 + name3 # if newname2 was changed then create a new construct name. 
            else:
                newconstructname = constructname
    
            return newconstructname, mutList, miscMutList
        
        elif constructname.find('eppcr') > -1 and commonvar is not None:
            
            splitpos = constructname.find('polyA') + 5
            name1 = constructname[:splitpos]
            name2 = ''
            name3= constructname[splitpos:]
            
            # some CDB entries have a space after polyA which will be captured in name3
            if name3[0] == ' ': 
                name3 = name3[1:]
        
            commonvar_match = mutstarRegex.findall(commonvar)
        
            for v in range(len(commonvar_match)):
                if v < len(commonvar_match) - 1:
                    name2 += commonvar_match[v] + ' '
                else:
                    name2 += commonvar_match[v]
                    
                mutList.append(commonvar_match[v])
                
            # replace asterisks with 'STOP'
            if name2.find('*') > -1:
                name2 = name2.replace('*', 'STOP')
                    
            newconstructname = name1 + ' ' + name2 + name3
            
            return newconstructname, mutList, miscMutList
            
        
        else: 
            # if no SDM sites were indicated in the constructname, just return 
            # everything in the miscMutList
            for t in range(0, num_ds_Matches):
                miscMutList.append(commonvar_match[t])
                
            
            return constructname, mutList, miscMutList
        
        
            
        
        
        # if len(newname2) > 1:
        #     return newconstructname, mutList, miscMutList
        # else:
        #     return None
       
        

# No code yet to handle map2ref disagreeing with denovo