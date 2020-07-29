#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 19:49:02 2020

@author: tsanga

Called by main.py
"""

# This function takes in one construct name containing an underscore or dash between mutations
# and converts it to a space. The new name is returned by the function. 

# Some Construct DB combo mutation entries list mutations separated by underscore or dash, 
# instead of a space. Examples:
    
#     N391X_T392X
#     N391X-T392X
#     D145X L146X


import re

def fixComboNames(constructname):
    uscoreRegex = re.compile(r'([A-Z]\d+X_)+') # search pattern for mutation with underscore (e.g. N391X_T392X)
    dashRegex = re.compile(r'([A-Z]\d+X-)+') #search pattern for mutations separated by dashes (e.g. N391X-T392X)
    mutRegex = re.compile(r'[A-Z]\d+[A-Z]') # search pattern for any mutation (e.g. N391X, T293Y)
    
    uscore_mo = uscoreRegex.search(constructname) #underscore regex match object
    dash_mo = dashRegex.search(constructname) # dash regex match object
    
    # initialize flag variables
    uflag = False
    dflag = False
    
    if uscore_mo is not None:
        uflag = True
        
    if dash_mo is not None:
        dflag = True
        
    if uflag == True or dflag == True:
        matches = mutRegex.findall(constructname) # returns a list of the unknown mutations
        numMatches = len(matches) # find the number of mutations
        
        firstMut = constructname.find(matches[0])
        lastMut = constructname.find(matches[numMatches-1])
        lastMutLen = len(matches[numMatches-1])
        lastMutChar = lastMut + lastMutLen

        name1 = constructname[:firstMut]
        name2 = constructname[firstMut:lastMutChar]
        name3 = constructname[lastMutChar:]
        
        if uflag == True:
            newname2 = name2.replace('_', ' ')
        elif dflag == True:
            newname2 = name2.replace('-', ' ')
        
        newconstructname = name1 + newname2 + name3
        
        return newconstructname

        
        # print('All of the mutations found: ' + str(matches))
        # print('first mutation position: ' + str(firstMut))
        # print('last mutation position: ' + str(lastMut))
        # print('character length of the last mutation: ' + str(lastMutLen))
        # print('name part 1: ' + name1)
        # print('name part 2: ' + name2)
        # print('name part 3: ' + name3)
        
        # print('The new name2 is: ' + newname2)
        # print('The new constructname is: ' + newconstructname)


            


