#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 15:45:03 2020

@author: tsanga

Runs under main.py
"""

scaffold_protein = {
    
    "376" : "METTVRYEQGSELTKTSSSPTADEPTIKIDDGRDEGNEQDSCSNTIRRKISPFVMSFGFRVFGVVLIIVDIIVVIVDLAISEKKRGIREILEGVSLAIALFFLVDVLMRVFVEGFKNYFRSKLNTLDAVIVVGTLLINMTYSFSDLAAFNSHNVYITADKQKNGIKANFTVRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQTVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYGGTGGSASQGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATIGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGKYKTRAVVKFEGDTLVNRIELKGTDFKEDGNILGHKLEYNTDQMPQMVTLLRVLRIVILIRIFRLASQKKQLEVVT",
    "421" : "MADVETETGMIAQWIVFAIMAAAAIAFGVAVHFRPSELKSAYYINIAICTIAATAYYAMAVNYQDLTMNGERQVVYARYINWVLTTPLLLLDLIVMTKMGGVMISWVIGADIFMIVFGILGAFEDEHKFKWVYFIAGCVMQAVLTYGMYNATWKDDLKKSPEYHSSYVSLLVFLSILWVFYPVVWAFGSGSGVLSVDNEAILMGILDVLAKPLFGMGCLIAHETIFKIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISGEPTTKSRITSEGEYIPLDQIDINVFCYENEV",
    "476" : "MADVETETGMIAQWIVFAIMAAAAIAFGVAVHFRPSELKSAYYINIAICTIAATAYYAMAVNYQDLTMNGERQVVYARYINWVLTTPLLLLDLIVMTKMGGVMISWVIGADIFMIVFGILGAFEDEHKFKWVYFIAGCVMQAVLTYGMYNATWKDDLKKSPEYHSSYVSLLVFLSILWVFYPVVWAFGSGSGVLSVDNEAILMGILDVLAKPLFGMGCLIAHETIFKIGTGFPFDPHYVEVLGERMHYVDVGPRDGTPVLFLHGNPTSSYVWRNIIPHVAPTHRCIAPDLIGMGKSDKPDLGYFFDDHVRFMDAFIEALGLEEVVLVIHDWGSALGFHWAKRNPERVKGIAFMEFIRPIPTWDEWPEFARETFQAFRTTDVGRKLIIDQNVFIEGTLPMGVVRPLTEVEMDHYREPFLNPVDREPLWRFPNELPIAGEPANIVALVEEYMDWLHQSPVPKLLFWGTPGVLIPPAEAARLAKSLPNCKAVDIGPGLNLLQEDNPDLIGSEIARWLSTLEISGEPTTKSRITSEGEYIPLDQIDINVFCYENEVQSQPILNTKEMAPQSKPPEELEMSSMPSPVAPLPARTEGVIDMRSMSSIDSFISCATDFPEATRF",
    "514" : "METDTLLLWVLLLWVPGSTGDRSESINFVSWGGSTQDAQKQAWADPFSKASGITVVQDGPTDYGKLKAMVESGNVQWDVVDVEADFALRAAAEGLLEPLDFSVIQRDKIDPRFVSDHGVGSFLFSFVLGYNEGKLGASKPQDWTALFDTKTYPGKRALYKWPSPGVLELALLADGVPADKLYPLDLDRAFKKLDTIKKDIVWWGGGAQSQQLLASGEVSMGQFWNGRIHALQEDGAPVGVSWKQNLVMADILVVPKGTKNKAAAMKFLASASSAKGQDDFSALTAYAPVNIDSVQRLDLAQVRITADKQKNGIMANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITLGMDELYKGGTGGSMSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNWNANLAPNLPTAYVKDQITLDFAYWAKNGPAIATRWNEWLVKLQVDEQKLISEEDLNAVGQDTQEVIVVPHSLPFKVVVISAILALVVLTIISLIILIMLWQKKPR",
    }


scaffold_dna = {
    
    "376" : "ATGGAGACGACTGTGAGGTATGAACAGGGGTCAGAGTTAACTAAAACTTCGAGCTCTCCAACAGCAGATGAGCCCACGATAAAGATTGATGATGGTCGTGATGAGGGTAATGAACAAGACAGCTGTTCCAATACCATTAGGAGAAAAATTTCCCCGTTTGTGATGTCATTTGGATTCAGAGTATTTGGAGTTGTGCTTATCATTGTAGACATCATAGTGGTGATTGTGGATCTGGCCATCAGTGAGAAGAAAAGAGGCATTAGAGAGATTCTTGAAGGTGTTTCCCTGGCTATAGCACTCTTCTTCCTTGTTGATGTTCTCATGAGAGTGTTTGTTGAAGGCTTCAAGAACTATTTCCGGTCCAAACTGAATACTTTGGATGCAGTCATAGTAGTGGGCACTCTGCTAATTAATATGACCTACTCCTTCTCTGACCTTGCTGCCTTTAACAGCCATAACGTGTATATTACCGCGGATAAACAGAAAAACGGCATTAAAGCGAACTTTACCGTGCGCCATAACGTGGAAGATGGCAGCGTGCAGCTGGCGGATCATTATCAGCAGAACACCCCGATTGGCGATGGCCCGGTGCTGCTGCCGGATAACCATTATCTGAGCACCCAGACCGTGCTGAGCAAAGATCCGAACGAAAAACGCGATCACATGGTGCTGCTGGAATTTGTGACCGCAGCGGGCATTACACACGGCATGGATGAACTGTATGGCGGCACCGGCGGCAGCGCGAGCCAGGGCGAAGAACTGTTTACCGGCGTGGTGCCGATTCTGGTGGAACTGGATGGCGATGTGAACGGCCATAAATTTAGCGTGCGCGGCGAAGGCGAAGGCGATGCGACCATTGGCAAACTGACCCTGAAATTTATTTGCACCACCGGCAAACTACCGGTGCCGTGGCCGACCCTGGTGACCACCTTAACCTATGGCGTGCAGTGCTTTAGCCGCTATCCGGATCATATGAAACGCCATGATTTTTTTAAAAGCGCGATGCCGGAAGGCTATGTGCAGGAACGCACCATTAGCTTTAAAGATGATGGCAAATATAAAACCCGCGCGGTGGTGAAATTTGAAGGCGATACCCTGGTGAACCGCATTGAACTGAAAGGCACCGATTTTAAAGAAGATGGCAACATTCTGGGGCATAAACTGGAATATAACACAGATCAGATGCCGCAGATGGTTACTCTTCTTCGAGTTCTGAGAATTGTTATCTTAATAAGAATATTTCGCCTGGCTAGCCAGAAGAAACAACTTGAAGTGGTAACCTAA",
    "421" : "ATGGCTGACGTGGAAACCGAGACCGGCATGATTGCACAGTGGATTGTCTTTGCTATTATGGCTGCTGCTGCTATTGCTTTTGGAGTGGCTGTGCACTTTCGGCCTTCAGAGCTGAAGAGCGCATACTATATCAACATTGCCATCTGCACTATCGCCGCTACCGCTTACTATGCAATGGCCGTGAACTACCAGGACCTGACAATGAATGGTGAAAGGCAGGTGGTCTACGCAAGATATATTAACTGGGTGCTGACCACACCACTGCTCCTGCTCGATCTCATCGTCATGACCAAGATGGGCGGAGTGATGATTTCTTGGGTCATCGGCGCAGACATTTTCATGATCGTGTTTGGTATTCTGGGCGCCTTCGAGGATGAACACAAGTTCAAATGGGTGTACTTTATCGCTGGATGTGTGATGCAGGCAGTCCTGACATACGGGATGTATAACGCCACTTGGAAAGACGATCTGAAGAAAAGCCCCGAGTACCATAGCTCCTATGTCAGTCTGCTCGTCTTCCTGTCAATCCTCTGGGTGTTTTATCCTGTCGTGTGGGCTTTCGGGTCTGGTAGTGGCGTGCTGTCCGTCGACAATGAGGCCATTCTCATGGGAATCCTGGATGTGCTCGCTAAGCCACTGTTTGGAATGGGGTGCCTCATTGCCCATGAGACTATCTTCAAGATCGGTACTGGCTTTCCATTCGACCCCCATTATGTGGAAGTCCTGGGCGAGCGCATGCACTACGTCGATGTTGGTCCGCGCGATGGCACCCCTGTGCTGTTCCTGCACGGTAACCCGACCTCCTCCTACGTGTGGCGCAACATCATCCCGCATGTTGCACCGACCCATCGCTGCATTGCTCCAGACCTGATCGGTATGGGCAAATCCGACAAACCAGACCTGGGTTATTTCTTCGACGACCACGTCCGCTTCATGGATGCCTTCATCGAAGCCCTGGGTCTGGAAGAGGTCGTCCTGGTCATTCACGACTGGGGCTCCGCTCTGGGTTTCCACTGGGCCAAGCGCAATCCAGAGCGCGTCAAAGGTATTGCATTTATGGAGTTCATCCGCCCTATCCCGACCTGGGACGAATGGCCAGAATTTGCCCGCGAGACCTTCCAGGCCTTCCGCACCACCGACGTCGGCCGCAAGCTGATCATCGATCAGAACGTTTTTATCGAGGGTACGCTGCCGATGGGTGTCGTCCGCCCGCTGACTGAAGTCGAGATGGACCATTACCGCGAGCCGTTCCTGAATCCTGTTGACCGCGAGCCACTGTGGCGCTTCCCAAACGAGCTGCCAATCGCCGGTGAGCCAGCGAACATCGTCGCGCTGGTCGAAGAATACATGGACTGGCTGCACCAGTCCCCTGTCCCGAAGCTGCTGTTCTGGGGCACCCCAGGCGTTCTGATCCCACCGGCCGAAGCCGCTCGCCTGGCCAAAAGCCTGCCTAACTGCAAGGCTGTGGACATCGGCCCGGGTCTGAATCTGCTGCAAGAAGACAACCCGGACCTGATCGGCAGCGAGATCGCGCGCTGGCTGTCGACGCTCGAGATTTCCGGCGAGCCAACCACTAAGAGCAGGATCACCAGCGAGGGCGAGTACATCCCCCTGGACCAGATCGACATCAACGTGTTCTGCTACGAGAACGAGGTGTAA",
    "476" : "ATGGCTGACGTGGAAACCGAGACCGGCATGATTGCACAGTGGATTGTCTTTGCTATTATGGCTGCTGCTGCTATTGCTTTTGGAGTGGCTGTGCACTTTCGGCCTTCAGAGCTGAAGAGCGCATACTATATCAACATTGCCATCTGCACTATCGCCGCTACCGCTTACTATGCAATGGCCGTGAACTACCAGGACCTGACAATGAATGGTGAAAGGCAGGTGGTCTACGCAAGATATATTAACTGGGTGCTGACCACACCACTGCTCCTGCTCGATCTCATCGTCATGACCAAGATGGGCGGAGTGATGATTTCTTGGGTCATCGGCGCAGACATTTTCATGATCGTGTTTGGTATTCTGGGCGCCTTCGAGGATGAACACAAGTTCAAATGGGTGTACTTTATCGCTGGATGTGTGATGCAGGCAGTCCTGACATACGGGATGTATAACGCCACTTGGAAAGACGATCTGAAGAAAAGCCCCGAGTACCATAGCTCCTATGTCAGTCTGCTCGTCTTCCTGTCAATCCTCTGGGTGTTTTATCCTGTCGTGTGGGCTTTCGGGTCTGGTAGTGGCGTGCTGTCCGTCGACAATGAGGCCATTCTCATGGGAATCCTGGATGTGCTCGCTAAGCCACTGTTTGGAATGGGGTGCCTCATTGCCCATGAGACTATCTTCAAGATCGGTACTGGCTTTCCATTCGACCCCCATTATGTGGAAGTCCTGGGCGAGCGCATGCACTACGTCGATGTTGGTCCGCGCGATGGCACCCCTGTGCTGTTCCTGCACGGTAACCCGACCTCCTCCTACGTGTGGCGCAACATCATCCCGCATGTTGCACCGACCCATCGCTGCATTGCTCCAGACCTGATCGGTATGGGCAAATCCGACAAACCAGACCTGGGTTATTTCTTCGACGACCACGTCCGCTTCATGGATGCCTTCATCGAAGCCCTGGGTCTGGAAGAGGTCGTCCTGGTCATTCACGACTGGGGCTCCGCTCTGGGTTTCCACTGGGCCAAGCGCAATCCAGAGCGCGTCAAAGGTATTGCATTTATGGAGTTCATCCACCCTATCCCGACCTGGGACGAATGGCCAGAATTTGCCCGCGAGACCTTCCAGGCCTTCCGCACCACCGACGTCGGCCGCAAGCTGATCATCGATCAGAACGTTTTTATCGAGGGTACGCTGCCGATGGGTGTCGTCCGCCCGCTGACTGAAGTCGAGATGGACCATTACCGCGAGCCGTTCCTGAATCCTGTTGACCGCGAGCCACTGTGGCGCTTCCCAAACGAGCTGCCAATCGCCGGTGAGCCAGCGAACATCGTCGCGCTGGTCGAAGAATACATGGACTGGCTGCACCAGTCCCCTGTCCCGAAGCTGCTGTTCTGGGGCACCCCAGGCGTTCTGATCCCACCGGCCGAAGCCGCTCGCCTGGCCAAAAGCCTGCCTAACTGCAAGGCTGTGGACATCGGCCCGGGTCTGAATCTGCTGCAAGAAGACAACCCGGACCTGATCGGCAGCGAGATCGCGCGCTGGCTGTCGACGCTCGAGATTTCCGGCGAGCCAACCACTAAGAGCAGGATCACCAGCGAGGGCGAGTACATCCCCCTGGACCAGATCGACATCAACGTGTTCTGCTACGAGAACGAGGTGCAAAGTCAGCCTATCCTGAACACAAAGGAAATGGCTCCACAGTCTAAGCCTCCCGAAGAGCTTGAGATGTCCAGTATGCCAAGTCCCGTGGCTCCCCTCCCTGCCAGGACTGAAGGAGTGATTGACATGAGGAGTATGTCATCTATTGATAGCTTCATCTCTTGCGCAACAGATTTCCCCGAGGCTACTCGATTCTAA",
    "514" : "ATGGAGACAGACACACTCCTGCTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGACAGATCTGAAAGCATCAATTTCGTGAGCTGGGGCGGTAGCACCCAGGATGCGCAGAAGCAGGCCTGGGCCGACCCGTTCAGCAAGGCCAGCGGCATTACCGTGGTCCAGGATGGGCCCACCGACTACGGCAAACTCAAGGCCATGGTCGAAAGCGGCAACGTGCAGTGGGACGTGGTGGATGTAGAGGCCGACTTCGCCTTGCGCGCCGCCGCCGAAGGCCTGCTCGAACCCCTGGATTTCTCGGTGATCCAGCGCGACAAGATCGACCCGCGCTTCGTCTCCGACCATGGCGTGGGCTCGTTCCTGTTCTCTTTCGTGCTCGGCTACAACGAAGGCAAGCTCGGAGCCAGCAAGCCCCAGGACTGGACCGCGCTGTTCGACACCAAGACCTACCCCGGCAAACGCGCCCTCTACAAATGGCCGAGCCCTGGCGTGCTCGAACTGGCCCTGCTCGCCGACGGCGTACCGGCCGACAAGCTCTACCCGCTGGACCTGGACCGCGCCTTCAAGAAACTCGACACCATCAAGAAAGACATCGTCTGGTGGGGCGGCGGTGCACAGTCGCAGCAGCTGCTGGCCTCCGGCGAAGTCAGCATGGGCCAGTTCTGGAACGGTCGCATCCACGCCCTGCAAGAAGACGGCGCGCCAGTGGGCGTGAGCTGGAAGCAGAACCTGGTGATGGCCGACATCCTCGTCGTGCCTAAAGGCACGAAGAACAAAGCCGCCGCGATGAAGTTCCTGGCCAGTGCCAGCAGCGCCAAAGGCCAGGACGACTTTTCCGCCCTGACCGCCTATGCTCCAGTGAACATCGACAGTGTGCAGCGCCTCGACCTTGCGCAGGTCCGGATCACCGCCGACAAGCAGAAGAACGGCATCATGGCGAACTTCAAGATCCGCCACAACGTGGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGTGCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGGGCGGTACCGGAGGGAGCATGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGCGCGGCGAGGGCGAGGGCGATGCCACCAACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCAGCTTCAAGGACGACGGCACCTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTGGAACGCCAACCTGGCACCGAACCTGCCGACCGCCTACGTCAAGGATCAGATCACCCTCGATTTCGCCTACTGGGCCAAGAACGGTCCGGCCATTGCGACACGGTGGAATGAATGGCTAGTCAAACTGCAGGTCGACGAACAAAAACTCATCTCAGAAGAGGATCTGAATGCTGTGGGCCAGGACACGCAGGAGGTCATCGTGGTGCCACACTCCTTGCCCTTTAAGGTGGTGGTGATCTCAGCCATCCTGGCCCTGGTGGTGCTCACCATCATCTCCCTTATCATCCTCATCATGCTTTGGCAGAAGAAGCCACGTTAG",
    
    }


# borrowed from:
# https://towardsdatascience.com/starting-off-in-bioinformatics-turning-dna-sequences-into-protein-sequences-c771dc20b89f
# DNA codon dictionary
protein = {"TTT" : "F", "CTT" : "L", "ATT" : "I", "GTT" : "V",
           "TTC" : "F", "CTC" : "L", "ATC" : "I", "GTC" : "V",
           "TTA" : "L", "CTA" : "L", "ATA" : "I", "GTA" : "V",
           "TTG" : "L", "CTG" : "L", "ATG" : "M", "GTG" : "V",
           "TCT" : "S", "CCT" : "P", "ACT" : "T", "GCT" : "A",
           "TCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "TCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "TCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "TAT" : "Y", "CAT" : "H", "AAT" : "N", "GAT" : "D",
           "TAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "TAA" : "*", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "TAG" : "*", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "TGT" : "C", "CGT" : "R", "AGT" : "S", "GGT" : "G",
           "TGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "TGA" : "*", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "TGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }


# these mutations occur in so many constructs that they are excluded from CDB comments to reduce noise. 
# they don't seem to affect sensor function, since 376.13 (one of the best ASAP hits) has these. 
# the commas and spaces are important for matching. 
# used in the polymorphisms function (under sequence_funcs.py)

suppressed_snps = {
    
    "376" : ['snp-t87c, ', 'snp-t1203c, ', 'snp-c1204t, ', 'snp-t1206a, ', 'snp-c1207a, '],
    "421" : ['snp-t690c, ']
    
    }