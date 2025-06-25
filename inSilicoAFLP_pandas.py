#----------------------------------------------------------------------------------------------------------------------------------------
#IMPORTS
#----------------------------------------------------------------------------------------------------------------------------------------

import sys
import os
import re
import itertools as it
import numpy as np
import pandas as pd

#----------------------------------------------------------------------------------------------------------------------------------------
#FUNCTIONS
#----------------------------------------------------------------------------------------------------------------------------------------
"""
def job_queue(job_definition, enzyme_data, adaptor_data, primer_data)
    #This function gets the dataframes with the job definitions, and the data of restriction
    #enzymes, adaptors and primers and uses it to fill the jobs dataframe which records Enzymes,
    #Adaptors primers and the recognition and extensions equences required for each.

    #This function uses a temporary list of lists: From the job definition the digestion, adaptor
    #and primer mix is put into a list and the degenerated selective bases are extended into
    #clear seqeunces. Over these clear sequences is iterated to create the individual jobs 
    #in the job dataframe.

    temp_jobs = []

    for job_entry in job_definition:
        method = job_entry[0]
        digestion_mix = job_entry[1].split("+")
        adaptor_mix = job_entry[2].split("+")
        primer_mix = job_entry[3].split("+")
        job = [method, digestion_mix, adaptor_mix, primer_mix]
        job_list.append(job)
"""
#----------------------------------------------------------------------------------------------------------------------------------------
#SCRIPT
#----------------------------------------------------------------------------------------------------------------------------------------

##Pandas dataframes

###Enzyme data

enzyme_data = pd.DataFrame(
    {"Enzyme name"                  : ["PstI", "EcoRI", "MseI"],
     "Recognition/Cut sequence:"    : ["CTGCA/G", "G/AATTC", "T/TAA"],
     "Recognition sequence"         : ["CTGCAG", "GAATTC", "TTAA"],
     "Cut Index"                    : [5, 1, 1]}
)

###Adaptor data

adaptor_data = pd.DataFrame(
    {"Adaptor name"                 : ["PstI Adaptor 1", "EcoRI Adaptor 1", "MseI Adaptor 1"],
     "Enzyme name"                  : ["PstI", "EcoRI", "MseI"],
     "Adaptor sequence 1"           : ["CTCGTAGACTGCGTACATGCA", "CTCGTAGACTGCGTACC", "GACGATGAGTCCTGAG"],
     "Adaptor sequence 2"           : ["CATCTGACGCATGT", "AATTGGTACGCAGTCTAC", "TACTCAGGACTCAT"]}
)

###Primer data

primer_data = pd.DataFrame(
    {"Primer name"                  : ["PstI Primer 1", "EcoRI Primer 1", "Mse1 Primer 1"],
     "Enzyme name"                  : ["PstI", "EcoRI", "MseI"],
     "Primer core sequence"         : ["GACTGCGTACATGCAG", "GACTGCGTACCAATTC", "GATGAGTCCTGAGTAA"],
     "Selective bases"              : ["R", "Y", "ATG"]}
)

###Jobs
job_definition = pd.DataFrame(
    {"Method"                       : ["single", "double", "triple"],
     "Digestion Enzymes"            : ["PstI", "PstI+EcoRI", "PstI+EcoRI+MseI"],
     "Adaptors"                     : ["PstI Adaptor 1", "PstI Adaptor 1+EcoRI Adaptor 1", "PstI Adaptor 1+EcoRI Adaptor 1"],
     "Primers"                      : ["PstI Primer 1", "PstI Primer 1+EcoRI Primer 1", "PstI Primer 1+EcoRI Primer 1"]}
)

jobs = pd.DataFrame(
    {"Job Number"                   : [],
     "Method"                       : [],
     "Name Enzyme 1"                : [],
     "Name Enzyme 2"                : [],
     "Name Enzyme 3"                : [],
     "Name Adaptor 1"               : [],
     "Name Adaptor 2"               : [],
     "Name Primer 1"                : [],
     "P1: 5' Recognition Sequences" : [],
     "P1: 3' Recognition Sequences" : [],
     "P1: 5' Extention Sequence"    : [],
     "P1: 3' Extention Sequence"    : [],
     "Name Primer 2"                : [],
     "P2: 5' Recognition Sequences" : [],
     "P2: 3' Recognition Sequences" : [],
     "P2: 5' Extention Sequence"    : [],
     "P2: 3' Extention Sequence"    : []}
)

##Data
path_genome = "Support\\ECPlambda.fasta"


##Script

#Print all dataframes as test
print("The enzyme data frame")
print(enzyme_data)
print("The adaptor data frame")
print(adaptor_data)
print("The primer data frame")
print(primer_data)
print("The job definition data frame")
print(job_definition)
print("The jobs data frame")
print(jobs)