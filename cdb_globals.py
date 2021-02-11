#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 16:58:39 2020

@author: tsanga
"""


### global declarations
dna_mutlist = [] # list to track all non-SDM mutations encountered 
constructs_updated = 0

### File info

# location of master_summary
ms_file_loc = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/master_summary_514.xlsx'
#location of CDB file
# cdb_database_loc = '/Users/tsanga/Documents/code/deepseq_update_cdb/ConstructDB_Downloads/514_original.xlsx'
cdb_database_loc ='/Users/tsanga/Documents/code/deepseq_update_cdb/ConstructDB_Downloads/514_F101X_is_L101X.xlsx'
# save the updated CDB file as:
cdb_output_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/output/ConstructDB_Category_514_annotated.xlsx'
#save the results_summary file as:
results_summary_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/output/Results_Summary_514.xlsx'

# # location of master_summary
# ms_file_loc = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/master_summary_514.xlsx'
# #location of CDB file
# cdb_database_loc = '/Users/tsanga/Documents/code/deepseq_update_cdb/ConstructDB_Downloads/514_sequential_numbering.xlsx'
# # save the updated CDB file as:
# cdb_output_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/output/ConstructDB_Category_514_sequential_annotated.xlsx'
# #save the results_summary file as:
# results_summary_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/output/Results_Summary_514_sequential.xlsx'


# # location of master_summary
# ms_file_loc = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/master_summary_421.xlsx'
# #location of CDB file
# cdb_database_loc = '/Users/tsanga/Documents/code/deepseq_update_cdb/ConstructDB_Downloads/421_original.xlsx'
# # save the updated CDB file as:
# cdb_output_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/output/ConstructDB_Category_421_annotated.xlsx'
# #save the results_summary file as:
# results_summary_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/output/Results_Summary_421.xlsx'

# # location of master_summary
# ms_file_loc = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/master_summary_376.xlsx'
# #location of CDB file
# cdb_database_loc = '/Users/tsanga/Documents/code/deepseq_update_cdb/ConstructDB_Downloads/376_original.xlsx'
# # save the updated CDB file as:
# cdb_output_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/output/ConstructDB_Category_376_annotated.xlsx'
# #save the results_summary file as:
# results_summary_file = '/Users/tsanga/Documents/Arthur_GENIE_Work/Deepseq/DeepSeqResults/20200807_DeepSequencing_full_run/output/Results_Summary_376.xlsx'
