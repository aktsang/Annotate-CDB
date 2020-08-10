#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:58:39 2020

@author: tsanga
"""

dna_mutlist = [] # list to track all non-SDM mutations encountered 
constructs_updated = 0

### File info

# location of master_summary
ms_file_loc = '/Users/tsanga/Documents/Arthur_GENIE_Work/Tn5/DeepSeqResults/20200807_DeepSequencing_full_run/master_summary_514.xlsx'


#location of CDB file
cdb_database_loc = '/Users/tsanga/Documents/code/deepseq_update_cdb/ConstructDB_Downloads/514_original.xlsx'

# save the updated CDB file as:
cdb_output_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Tn5/DeepSeqResults/20200807_DeepSequencing_full_run/output/ConstructDB_Category_514_annotated.xlsx'

#save the results_summary file as:
results_summary_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Tn5/DeepSeqResults/20200807_DeepSequencing_full_run/output/Results_Summary_514.xlsx'