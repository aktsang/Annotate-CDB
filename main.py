    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 08:25:17 2020

@author: tsanga
"""


# PURPOSE:
#     Clone information was downloaded to a spreadsheet from the GENIE Construct Database.
#     Update the construct database spreadsheet in a new excel file with results from 
#     deep sequencing shown in the master summary file. 


import re
import os
import openpyxl

from openpyxl import Workbook
from openpyxl import load_workbook
from openpyxl.styles import PatternFill
from openpyxl.styles.colors import Color
from Bio import pairwise2

from cdb_globals import *


### file inputs

# location of sequencing results file from Quantitative Genomics
ms_file_loc = '/Users/tsanga/Documents/code/deepseq_update_cdb/testfiles/master_summary_QG_20190116_copy.xlsx'

#location of CDB file
cdb_database_loc = '/Users/tsanga/Documents/code/deepseq_update_cdb/testfiles/ConstructDB_combined_testing.xlsx'

# information on the master_summary from Quantitative Genomics
samplename_col = 'A'
plate_col = 'B'
well_col = 'C'
map2ref_var_col = 'V'
denovo_var_col = 'W'
common_var_col = 'X'
map2ref_seq_col = 'AC'
denovo_seq_col = 'AE'

# information on the CDB Lookup Table
constructname_col = 'A'
sequence_col = 'B'
construct_id_col = 'D'
categorynumber_col = 'E'
construct_loc_col = 'F'
comments_col = 'I'


# load lookup workbooks
master_summary_file = load_workbook(ms_file_loc)
ms= master_summary_file.active
msrows = ms.max_row + 1

cdb_file = load_workbook(cdb_database_loc)
cdb = cdb_file.active
cdbrows = cdb.max_row + 1


### Results workbook initialization
results_wb = Workbook()

# create worksheets for each type of result
results_ws0 = results_wb.active # the default active sheet
results_ws1 = results_wb.create_sheet("Combos Renamed")
results_ws2 = results_wb.create_sheet("Frameshifts")
results_ws3 = results_wb.create_sheet("Unidentical Sequences")
results_ws4 = results_wb.create_sheet("Mutation Count")

# intial values of the excel results summary
results_ws1_row = 1
results_ws2_row = 1
results_ws3_row = 1
results_ws4_row = 4

results_ws0.title = "Summary"
results_ws0['A1'].value = 'Worksheet name'
results_ws0['B1'].value = 'Description'
results_ws0['A2'].value = 'Combos Renamed'
results_ws0['B2'].value = 'Combo mutation names that contained underscores or dashes, which were corrected to spaces.'
results_ws0['A3'].value = 'Frameshifts'
results_ws0['B3'].value = 'Frameshifts found in the master summary data.'
results_ws0['A4'].value = 'Mapped and Denovo Mismatch'
results_ws0['B4'].value = 'Called mutations with mismatches between mapped and de novo assembly gene sequences.'
results_ws0['A5'].value = 'Mutation Count'
results_ws0['B5'].value = 'All non-SDM mutations found in the dataset.'

results_ws1['A1'].value = 'Old construct name'
results_ws1['B1'].value = 'New construct name'

results_ws2['A1'].value = 'Construct'
results_ws2['B1'].value = 'Frameshift identified'
results_ws2['C1'].value = 'Comment'

results_ws3['A1'].value = 'Construct'
results_ws3['B1'].value = 'common_var'
results_ws3['C1'].value = 'map2ref'
results_ws3['D1'].value = 'denovo'

results_ws4['A1'].value = 'Total number of entries: '
results_ws4['A2'].value = 'Number of constructs updated: '
results_ws4['A4'].value = 'Mutation'
results_ws4['B4'].value = 'Number of constructs containing'
# the total number of constructs processed is added to the workbook at the end of the program. 

# newNames_wb = Workbook()
# newNames_ws = newNames_wb.active

### fix underscores and dashes in mutation list

from fix_combo_names import fixComboNames
for i in range(2, cdbrows):
    constructname = cdb[samplename_col + str(i)].value
    newconstructname = fixComboNames(constructname)
    
    if newconstructname is not None: # if a change was made by fixComboNames
        cdb[samplename_col + str(i)].value = newconstructname # updates the CDB name
        
        # The following code summarizes the changes in a spreadsheet. 
        results_ws1_row = results_ws1.max_row + 1
        results_ws1['A' + str(results_ws1_row)].value = constructname
        results_ws1['B' + str(results_ws1_row)].value = newconstructname

# results_ws1.save('/Users/tsanga/Documents/code/deepseq_update_cdb/testfiles/output/results_summary.xlsx')    
# results_ws1.close()
# newNames_wb.save('/Users/tsanga/Documents/code/deepseq_update_cdb/testfiles/output/Fixed_Combo_Names_summary.xlsx')    
# newNames_wb.close()

print('Fix combo names complete.')


### Update the ConstructDB list using the master summary sequencing data. 

platewellRegex = re.compile(r'A\d\d\d\d_[A-H]\d\d')
    
# loop through master summary 
for j in range(2, msrows):

    platename = ms[plate_col + str(j)].value # get plate name from row j
    wellname = ms[well_col + str(j)].value # get well name from row j
    
    searchtext = platename + '_' + wellname # the plate_well we're looking for, e.g. A0123_A01

    # Search the CDB list for the plate and well listing.
    for k in range(2, cdbrows):
        pwmatch = platewellRegex.search(str(cdb[construct_loc_col + str(k)].value)) 
        if pwmatch is not None:
            cdb_loc_compare = pwmatch.group() # extract the text of the matching pattern
            if cdb_loc_compare == searchtext: # compare the matched text to the desired plate and well
                #print('Found match for ' + searchtext + ' at index ' + str(k) + ' (' + cdb_loc_compare + ')')
                cdbindex = k # save the row number of this entry
                break # exit cdb search loop when the entry is found. 

    
    ### Check for existing Sanger sequence
    from sequence_funcs import checkSequence
    seqFound = checkSequence(cdb[sequence_col + str(cdbindex)].value)
    
    if seqFound == False: # no Sanger sequence found, so update the record

        ### Name and Sequence update
        constructname = cdb[samplename_col + str(cdbindex)].value
        map2refvar = ms[map2ref_var_col + str(j)].value
        denovovar = ms[denovo_var_col + str(j)].value
        commonvar = ms[common_var_col + str(j)].value
        
        if commonvar is not None: # For now, only work on records with a called mutation
            
            ### update name
            from updateMutants import updateMutations
            newconstructname = updateMutations(constructname, map2refvar, denovovar, commonvar)
            
            if newconstructname is not None:
                cdb[samplename_col + str(cdbindex)].value = newconstructname[0]
                
                # TESTING ONLY - color the changed construct names yellow
                cdb[samplename_col + str(cdbindex)].fill = PatternFill(fill_type="solid", fgColor="ffff00")
                cdb[comments_col + str(cdbindex)].fill = PatternFill(fill_type="solid", fgColor="ffff00")
                print('Name updated: ' + newconstructname[0])
                constructs_updated += 1
                
            ### udpate sequence
                # read map2ref and denovo sequences
                map2refseq = ms[map2ref_seq_col + str(j)].value
                map2refseq = map2refseq.upper() # convert to all upper case
                denovoseq = ms[denovo_seq_col + str(j)].value
                denovoseq = denovoseq.upper() # convert to all upper case
                # pass to function
                from sequence_funcs import compareSeq
                seqident = compareSeq(map2refseq, denovoseq)
                if seqident == True:
                    cdb[sequence_col + str(cdbindex)].value = denovoseq # if assembled sequences are identical, plug in denovo sequence
                    # make the cell light blue
                    cdb[sequence_col + str(cdbindex)].fill = PatternFill(fill_type="solid", fgColor="00e5ff")
                else:
                    # common var called but sequences are different.
                    # make the cell orange
                    cdb[sequence_col + str(cdbindex)].fill = PatternFill(fill_type="solid", fgColor="ffbf00")
                    # add to spreadsheet
                    results_ws3_row = results_ws3.max_row + 1
                    results_ws3['A' + str(results_ws3_row)].value = constructname
                    results_ws3['B' + str(results_ws3_row)].value = commonvar
                    results_ws3['C' + str(results_ws3_row)].value = map2refseq
                    results_ws3['D' + str(results_ws3_row)].value = denovoseq
                    
                ### Annotation functions
                origComment = cdb[comments_col + str(cdbindex)].value
                goodMutations = newconstructname[1]
                extraMutations = newconstructname[2]
                scaffold = str(cdb[categorynumber_col + str(cdbindex)].value)
                
                from annotations import annotate
                finalComment = annotate(origComment, goodMutations, extraMutations, map2refvar, denovovar, denovoseq, scaffold)
                
                from sequence_funcs import findFrameshift
                frameshiftList = findFrameshift(denovovar)
                if frameshiftList is not None:
                    if frameshiftList[0] == True:
                        finalComment += ' ' + 'frameshift-' + frameshiftList[1] + ', '
                        
                from sequence_funcs import polymorphisms
                other_polymorphisms = polymorphisms(scaffold, denovoseq, goodMutations)
                if other_polymorphisms is not None:
                    finalComment += other_polymorphisms
                        
                cdb[comments_col + str(cdbindex)].value = finalComment

        
        else: # commonvar is None
        
            ### Frameshifts without a commonvar    
            origComment = cdb[comments_col + str(cdbindex)].value
        
            from sequence_funcs import findFrameshift
            frameshiftList = findFrameshift(denovovar)
            if frameshiftList is not None:
                if frameshiftList[0] == True:
                    finalComment = origComment + ': frameshift-' + frameshiftList[1] + ', ' # simply states the residue at which the frameshift starts.
            else:
                finalComment = origComment
                    
            # If no commonvar exists, no need to look for problems besides frameshifts.
            
            # from sequence_funcs import polymorphisms
            # other_polymorphisms = polymorphisms(scaffold, denovoseq, goodMutations)
            # if other_polymorphisms is not None:
            #     finalComment += other_polymorphisms
            
            cdb[comments_col + str(cdbindex)].value = finalComment
                    
        
        ### save frameshift results to results_ws2
        if frameshiftList is not None:
            results_ws2_row = results_ws2.max_row + 1
            results_ws2['A' + str(results_ws2_row)].value = constructname
            results_ws2['B' + str(results_ws2_row)].value = denovovar
            results_ws2['C' + str(results_ws2_row)].value = finalComment
            
        
        
        ### Other single nucleotide polymorphisms and insertion/deletions:
        # from sequence_funcs import polymorphisms
        # other_polymorphisms = polymorphisms(scaffold, denovoseq, goodMutations)
        # if other_polymorphisms is not None:
        #     finalComment += other_polymorphisms
    
    
    elif seqFound == True: # found a Sanger sequence, don't modify to CDB record
        # TESTING ONLY - color the unchanged construct records green
        cdb[samplename_col + str(cdbindex)].fill = PatternFill(fill_type="solid", fgColor="48ff00")
        cdb[sequence_col + str(cdbindex)].fill = PatternFill(fill_type="solid", fgColor="48ff00")
            
        
# tally dna mutations in dna_mutlist

# first traverse the dna_mutlist and add unique values to a new list
unique_muts = []
for l in range(0, len(dna_mutlist)):
    if dna_mutlist[l] not in unique_muts:
        unique_muts.append(dna_mutlist[l])

# add up the mutations
for m in range(0, len(unique_muts)):
    occurrences = dna_mutlist.count(unique_muts[m])
    results_ws4_row += 1
    results_ws4['A' + str(results_ws4_row)].value = unique_muts[m]
    results_ws4['B' + str(results_ws4_row)].value = occurrences

results_ws4['B1'].value = j+1 # fill in the number of constructs processed. 
results_ws4['B2'].value = constructs_updated # fill in the number of construct names changed

cdb_file.save('/Users/tsanga/Documents/code/deepseq_update_cdb/testfiles/output/construct-db_updated_test.xlsx')
cdb_file.close()
results_wb.save('/Users/tsanga/Documents/code/deepseq_update_cdb/testfiles/output/Results_Summary.xlsx')
results_wb.close()
master_summary_file.close()

print('The operation has completed.')