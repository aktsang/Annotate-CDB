#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v. 20200717_0039
Created on Mon Jul 13 19:15:40 2020

@author: tsanga

Called by main.py
"""


import re


def updateMutations(constructname, map2refvar, denovovar, commonvar):
    
    mutRegex = re.compile(r'[A-Z]\d+[A-Z]|[A-Z]\d+\*') # search pattern for any mutation (e.g. N391X, T293Y, L405*)
    aaRegex = re.compile(r'[A-Z]\d+') # search pattern for residue number
    
# TO DO Look for a Sanger sequence and mutation in CDB. If it exists, do nothing to the clone name
    
    if commonvar is not None:
        cname_match = mutRegex.findall(constructname) # get list of mutations in CDB
        commonvar_match = mutRegex.findall(commonvar) # get list of mutations called by DeepSeq
        
        num_ds_Matches = len(commonvar_match) # get number of mutations found
        num_cdb_Matches = len(cname_match) # get number of mutations expected
        
        # string splitting on constructname
        firstMut = constructname.find(cname_match[0])
        lastMut = constructname.find(cname_match[num_cdb_Matches-1])
        lastMutLength = len(cname_match[num_cdb_Matches-1])
        lastMutChar = lastMut + lastMutLength        

        name1 = constructname[:firstMut]
        name2 = constructname[firstMut:lastMutChar]
        name3 = constructname[lastMutChar:]
        
        ### Compare each residue from deep seq to each residue in cdb
        # Deep seq may find multiple mutations. Is any of them the one we expect? 
        
        cname_aa_match = aaRegex.findall(name2) # find residue numbers from name2 fragment of constructname
        commonvar_aa_match = aaRegex.findall(commonvar) # find residue numbers from commonvar
        
        newname2 = '' # the new construct name with sequenced mutations
        mutList = [] # a list of the updated (desired) mutations
        miscMutList = [] # a list of unexpected mutations
        
        
        for r in range(0, num_ds_Matches): # for each deep seq identified mutation...
            for s in range(0, len(cname_aa_match)): # for each construct db intended mutation...
                if commonvar_aa_match[r] == cname_aa_match[s]: # if the residue number (e.g. L405) is the same...
                    newname2 += commonvar_match[r] # add the deep seq identified mutation to the name.
                    mutList.append(commonvar_match[r]) # store all the deep seq identified mutations in mutList
                else:
                    miscMutList.append(commonvar_match[r]) # store all the deep seq identified mutations in mutList
                    
        if len(newname2) > 1:
            newconstructname = name1 + newname2 + name3 # if newname2 was changed then create a new construct name. 

        
        if len(newname2) > 1:
            return newconstructname, mutList, miscMutList
        else:
            return None
       
        

# No code yet to handle map2ref disagreeing with denovo