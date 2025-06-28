#----------------------------------------------------------------------------------------------------------------------------------------
#IMPORTS
#----------------------------------------------------------------------------------------------------------------------------------------

import sys
import os
import re
import itertools as it
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import math as m

#----------------------------------------------------------------------------------------------------------------------------------------
#FUNCTIONS
#----------------------------------------------------------------------------------------------------------------------------------------

def expand_selective_bases(sequence):
    #This function takes a seqeunce of nuclotides according to IUPAC numencalture
    #and returns a list of all seqeunces they represent. 

    IUPAC_decoder = {'A': ['A'],'C': ['C'],'G': ['G'],'T': ['T'],'U': ['T'],'R': ['A', 'G'],'Y': ['C', 'T'],'S': ['G', 'C'],'W': ['A', 'T'],'K': ['G', 'T'],'M': ['A', 'C'],'B': ['C', 'G', 'T'],'D': ['A', 'G', 'T'],'H': ['A', 'C', 'T'],'V': ['A', 'C', 'G'],'N': ['A', 'C', 'G', 'T']}
    options = [IUPAC_decoder[base.upper()] for base in sequence]
    sequence_list = [''.join(p) for p in it.product(*options)]

    return sequence_list

def five_prime_end(recognition_sequence, cut_index, selective_bases):
    #This function takes a recognition sequence and cut index of a Restriction enyzme and selective bases and creates the matching seqeunces for the 5' end
    #of fragments

    rec_site_fragment = recognition_sequence[cut_index:len(recognition_sequence)]
    five_prime_end = rec_site_fragment + selective_bases

    return five_prime_end

def three_prime_end(recognition_seqeunce, cut_index, selective_bases):
    #This function takes a recognition sequence and cut index of a restriction enyzme and selective bases and creates the matching seqeunces for the 3' end
    #of fragments
    
    rec_site_fragment = recognition_seqeunce[0:cut_index]
    complement_selective_bases = selective_bases.translate(str.maketrans("ATGC", "TACG"))
    three_prime_end = complement_selective_bases + rec_site_fragment

    return three_prime_end

def five_primer_extention(primer_sequence, cut_index):
    #This fucntion takes primer and recognition sequence and cut index to generate the 5' and 3' extention of the amplicon

    five_prime_extension = primer_sequence[0:(len(primer_sequence)-cut_index)]

    return five_prime_extension

def three_primer_extention(primer_sequence, cut_index):
    #This function takes primer and recognition sequence and cut index to generate the 5' and 3' extention of the amplicon

    reversed_primer = primer_sequence[::-1]
    reversed_complement = reversed_primer.translate(str.maketrans("ATGC", "TACG"))
    three_primer_extension = reversed_complement[cut_index:]

    return three_primer_extension

def job_queue(job_definition, enzyme_data, adaptor_data, primer_data):
    #This function gets the dataframes with the job definitions, and the data of restriction
    #enzymes, adaptors and primers and uses it to fill the jobs dataframe which records Enzymes,
    #Adaptors primers and the recognition and extensions equences required for each.

    #This function uses a temporary list of lists: From the job definition the digestion, adaptor
    #and primer mix is put into a list and the degenerated selective bases are extended into
    #clear seqeunces. Over these clear sequences is iterated to create the individual jobs 
    #in the job dataframe.

    temp_jobs = []

    for index, job_entry in job_definition.iterrows():
        Method = job_entry["Method"]
        Digestion_Enzymes = job_entry["Digestion Enzymes"].split("+")
        Adaptors = job_entry["Adaptors"].split("+")
        Primers = job_entry["Primers"].split("+")
        temp_job_entry = [Method, Digestion_Enzymes, Adaptors, Primers]
        temp_jobs.append(temp_job_entry)
    
    #decoding the selective bases
    for temp_job_entry in temp_jobs:
        for primer in temp_job_entry[3]:
            primer_name = primer
            #print(primer_name, type(primer))
            primer_index = primer_data[primer_data["Primer name"] == primer].index[0]
            #print(primer_index, type(primer_index))
            degen_sel_bases = primer_data.at[primer_index, "Selective bases"]
            #print(degen_sel_bases, type(degen_sel_bases))
            expanded_sel_bases = expand_selective_bases(degen_sel_bases)
            index = temp_job_entry[3].index(primer)
            temp_job_entry[3][index] = [primer_name, degen_sel_bases, expanded_sel_bases]
    print(temp_jobs)

    #Jobs data frame

    jobs = pd.DataFrame(
    {"Job Number"                   : [],
     "Method"                       : [],
     "Name Enzyme 1"                : [],
     "Name Enzyme 2"                : [],
     "Name Enzyme 3"                : [],
     "Name Adaptor 1"               : [],
     "Name Adaptor 2"               : [],
     "Name Primer 1"                : [],
     "P1: Selective Bases"           : [],
     "P1: 5' Recognition Sequences" : [],
     "P1: 3' Recognition Sequences" : [],
     "P1: 5' Extention Sequence"    : [],
     "P1: 3' Extention Sequence"    : [],
     "Name Primer 2"                : [],
     "P2: Selective Bases"          : [],
     "P2: 5' Recognition Sequences" : [],
     "P2: 3' Recognition Sequences" : [],
     "P2: 5' Extention Sequence"    : [],
     "P2: 3' Extention Sequence"    : []}
    )

    #Fill job dataframe

    for temp_job_entry in temp_jobs:
        if temp_job_entry[0] == "single":
            method = temp_job_entry[0]
            #print("Method of Job", method)
            enzyme = temp_job_entry[1][0]
            #print("Enzyme of Job", enzyme)
            name_adaptor = temp_job_entry[2][0]
            #print("adaptor of Job", name_adaptor)
            name_primer = temp_job_entry[3][0][0]
            #print("Primer of job", name_primer)
            for selective_base in temp_job_entry[3][0][2]:
                RE_index = enzyme_data[enzyme_data["Enzyme name"] == enzyme].index[0]
                #print("RE index", RE_index)
                recognition_sequence = enzyme_data.at[RE_index, "Recognition sequence"]
                #print("RE_seqeunce", recognition_sequence)
                cut_index = enzyme_data.at[RE_index, "Cut Index"]
                #print("cut index", cut_index)
                five_rec_seq = five_prime_end(recognition_sequence, cut_index, selective_base)
                #print("Five rec seq", five_rec_seq)
                three_rec_seq = three_prime_end(recognition_sequence, cut_index, selective_base)
                #print("Three rec seq", three_rec_seq)
                primer_index = primer_data[primer_data["Primer name"] == name_primer].index[0]
                primer_enzyme = primer_data.at[primer_index, "Enzyme name"]
                primer_sequence = primer_data.at[primer_index, "Primer core sequence"]
                five_extention_seq = five_primer_extention(primer_sequence, cut_index)
                #print("Five extention seq", five_extention_seq)
                three_extention_seq = three_primer_extention(primer_sequence, cut_index)
                #print(three_extention_seq)
                new_job = pd.DataFrame(
                                    {"Job Number"                  : ["none"],
                                    "Method"                       : [method],
                                    "Name Enzyme 1"                : [enzyme],
                                    "Name Enzyme 2"                : ["none"],
                                    "Name Enzyme 3"                : ["none"],
                                    "Name Adaptor 1"               : [name_adaptor],
                                    "Name Adaptor 2"               : ["none"],
                                    "Name Primer 1"                : [name_primer],
                                    "P1: Enzyme"                   : [primer_enzyme],
                                    "P1: Selective Bases"          : [selective_base],
                                    "P1: 5' Recognition Sequences" : [five_rec_seq],
                                    "P1: 3' Recognition Sequences" : [three_rec_seq],
                                    "P1: 5' Extention Sequence"    : [five_extention_seq],
                                    "P1: 3' Extention Sequence"    : [three_extention_seq],
                                    "Name Primer 2"                : ["none"],
                                    "P2: Enzyme"                   : [primer_enzyme],
                                    "P2: Selective Bases"          : ["none"],
                                    "P2: 5' Recognition Sequences" : ["none"],
                                    "P2: 3' Recognition Sequences" : ["none"],
                                    "P2: 5' Extention Sequence"    : ["none"],
                                    "P2: 3' Extention Sequence"    : ["none"],
                                    "Number of Overlaps"           : [0]}
                                )
                #print(new_job)
                jobs = pd.concat([jobs, new_job], ignore_index=True)
        elif temp_job_entry[0] == "double":
            #print(temp_job_entry)
            method = temp_job_entry[0]
            #print("Method of Job", method)
            enzyme_list = temp_job_entry[1]
            print("Enzymes of Job", enzyme)
            name_adaptor = temp_job_entry[2]
            #print("Adaptors of Job", name_adaptor)
            name_primer = [temp_job_entry[3][0][0], temp_job_entry[3][1][0]]
            #print("Primers of job", name_primer)
            #Create combinations of selective bases
            all_primers = []
            for primer in temp_job_entry[3]:
                enzyme = primer[0]
                for selective_base in primer[2]:
                    comb = (enzyme, selective_base)
                    all_primers.append(comb)
            #print("All primers")
            #print(all_primers)
            all_primer_pairs = list(it.combinations(all_primers, 2))
            #print("all_primer_pairs")
            #print(all_primer_pairs)
            filtered_combinations = []
            for entry in all_primer_pairs:
                if entry[0][0] != entry[1][0]:
                    filtered_combinations.append(entry)
            #print("filtered combinations")
            #print(filtered_combinations)
            for combination in filtered_combinations:
                coll_data = []
                for primer in combination:
                    #print(f"Primer: {primer}")
                    primer_name = primer[0]
                    selective_base = primer[1]
                    primer_index = primer_data[primer_data["Primer name"] == primer_name].index[0]
                    primer_enzyme = primer_data.at[primer_index, "Enzyme name"]
                    #print(f"primer enzyme: {primer_enzyme}")
                    RE_index = enzyme_data[enzyme_data["Enzyme name"] == primer_enzyme].index[0]
                    recognition_sequence = enzyme_data.at[RE_index, "Recognition sequence"]
                    cut_index = enzyme_data.at[RE_index, "Cut Index"]
                    five_rec_seq = five_prime_end(recognition_sequence, cut_index, selective_base)
                    three_rec_seq = three_prime_end(recognition_sequence, cut_index, selective_base)
                    primer_sequence = primer_data.at[primer_index, "Primer core sequence"]
                    primer_enzyme = primer_data.at[primer_index, "Enzyme name"]
                    five_extention_seq = five_primer_extention(primer_sequence, cut_index)
                    three_extention_seq = three_primer_extention(primer_sequence, cut_index)
                    new_entry = [primer_name, selective_base, five_rec_seq, three_rec_seq, five_extention_seq, three_extention_seq, primer_enzyme]
                    coll_data.append(new_entry)
                #print("Collected data")
                #print(coll_data)
                new_job = pd.DataFrame(
                                    {"Job Number"                  : ["none"],
                                    "Method"                       : [method],
                                    "Name Enzyme 1"                : [enzyme_list[0]],
                                    "Name Enzyme 2"                : [enzyme_list[1]],
                                    "Name Enzyme 3"                : ["none"],
                                    "Name Adaptor 1"               : [name_adaptor[0]],
                                    "Name Adaptor 2"               : [name_adaptor[1]],
                                    "Name Primer 1"                : [coll_data[0][0]],
                                    "P1: Enzyme"                   : [coll_data[0][6]],
                                    "P1: Selective Bases"          : [coll_data[0][1]],
                                    "P1: 5' Recognition Sequences" : [coll_data[0][2]],
                                    "P1: 3' Recognition Sequences" : [coll_data[0][3]],
                                    "P1: 5' Extention Sequence"    : [coll_data[0][4]],
                                    "P1: 3' Extention Sequence"    : [coll_data[0][5]],
                                    "Name Primer 2"                : [coll_data[1][0]],
                                    "P2: Enzyme"                   : [coll_data[1][6]],
                                    "P2: Selective Bases"          : [coll_data[1][1]],
                                    "P2: 5' Recognition Sequences" : [coll_data[1][2]],
                                    "P2: 3' Recognition Sequences" : [coll_data[1][3]],
                                    "P2: 5' Extention Sequence"    : [coll_data[1][4]],
                                    "P2: 3' Extention Sequence"    : [coll_data[1][5]],
                                    "Number of Overlaps"           : [0]}
                                )
                #print(new_job)
                jobs = pd.concat([jobs, new_job], ignore_index=True)
        elif temp_job_entry[0] == "triple":
            #print(temp_job_entry)
            method = temp_job_entry[0]
            #print("Method of Job", method)
            enzyme_list = temp_job_entry[1]
            #print("Enzymes of Job", enzyme)
            name_adaptor = temp_job_entry[2]
            #print("Adaptors of Job", name_adaptor)
            name_primer = [temp_job_entry[3][0][0], temp_job_entry[3][1][0]]
            #print("Primers of job", name_primer)
            #Create combinations of selective bases
            all_primers = []
            for primer in temp_job_entry[3]:
                enzyme = primer[0]
                for selective_base in primer[2]:
                    comb = (enzyme, selective_base)
                    all_primers.append(comb)
            #print("All primers")
            #print(all_primers)
            all_primer_pairs = list(it.combinations(all_primers, 2))
            #print("all_primer_pairs")
            #print(all_primer_pairs)
            filtered_combinations = []
            for entry in all_primer_pairs:
                if entry[0][0] != entry[1][0]:
                    filtered_combinations.append(entry)
            #print("filtered combinations")
            #print(filtered_combinations)
            for combination in filtered_combinations:
                coll_data = []
                for primer in combination:
                    #print(f"Primer: {primer}")
                    primer_name = primer[0]
                    selective_base = primer[1]
                    primer_index = primer_data[primer_data["Primer name"] == primer_name].index[0]
                    primer_enzyme = primer_data.at[primer_index, "Enzyme name"]
                    #print(f"primer enzyme: {primer_enzyme}")
                    RE_index = enzyme_data[enzyme_data["Enzyme name"] == primer_enzyme].index[0]
                    recognition_sequence = enzyme_data.at[RE_index, "Recognition sequence"]
                    cut_index = enzyme_data.at[RE_index, "Cut Index"]
                    five_rec_seq = five_prime_end(recognition_sequence, cut_index, selective_base)
                    three_rec_seq = three_prime_end(recognition_sequence, cut_index, selective_base)
                    primer_sequence = primer_data.at[primer_index, "Primer core sequence"]
                    primer_enzyme = primer_data.at[primer_index, "Enzyme name"]
                    five_extention_seq = five_primer_extention(primer_sequence, cut_index)
                    three_extention_seq = three_primer_extention(primer_sequence, cut_index)
                    new_entry = [primer_name, selective_base, five_rec_seq, three_rec_seq, five_extention_seq, three_extention_seq, primer_enzyme]
                    coll_data.append(new_entry)
                #print("Collected data")
                #print(coll_data)
                new_job = pd.DataFrame(
                                    {"Job Number"                  : ["none"],
                                    "Method"                       : [method],
                                    "Name Enzyme 1"                : [enzyme_list[0]],
                                    "Name Enzyme 2"                : [enzyme_list[1]],
                                    "Name Enzyme 3"                : [enzyme_list[2]],
                                    "Name Adaptor 1"               : [name_adaptor[0]],
                                    "Name Adaptor 2"               : [name_adaptor[1]],
                                    "Name Primer 1"                : [coll_data[0][0]],
                                    "P1: Enzyme"                   : [coll_data[0][6]],
                                    "P1: Selective Bases"          : [coll_data[0][1]],
                                    "P1: 5' Recognition Sequences" : [coll_data[0][2]],
                                    "P1: 3' Recognition Sequences" : [coll_data[0][3]],
                                    "P1: 5' Extention Sequence"    : [coll_data[0][4]],
                                    "P1: 3' Extention Sequence"    : [coll_data[0][5]],
                                    "Name Primer 2"                : [coll_data[1][0]],
                                    "P1: Enzyme"                   : [coll_data[1][6]],
                                    "P2: Selective Bases"          : [coll_data[1][1]],
                                    "P2: 5' Recognition Sequences" : [coll_data[1][2]],
                                    "P2: 3' Recognition Sequences" : [coll_data[1][3]],
                                    "P2: 5' Extention Sequence"    : [coll_data[1][4]],
                                    "P2: 3' Extention Sequence"    : [coll_data[1][5]],
                                    "Number of Overlaps"           : [0]}
                                )
                #print(new_job)
                jobs = pd.concat([jobs, new_job], ignore_index=True)

        jobs["Job Number"] = jobs.apply(lambda row: f"Job {row.name + 1}", axis = 1)
    return(jobs)

def sequence_fasta(filename):
    #This function extracts a nuclotide sequence from a fasta file.
    #It checks and creates errors when the provided file path is wrong, the file is not a fasta file, does not contain a DNA sequence, 
    #is empty or contains multiple sequences.

    if not os.path.isfile(filename):
        sys.exit(f"Error: {filename} does not exist")

    with open(filename, "r") as file:
        first_line = file.readline()
        if not first_line.startswith(">"):
            sys.exit(f"Error: {filename} is not a fasta file")
    
    sequence = ""

    with open(filename, "r") as file:
        next(file)
        for line in file:
            line = line.strip()
            line = line.upper()
            if line.startswith(">"):
                sys.exit(f"Error: {filename} contains multiple sequences.")
            elif not re.fullmatch("[ATGCN]*", line):
                sys.exit(f"Error: {filename} does not contain a nuclotide sequence.")
            sequence += line

    if len(sequence) == 0:
        sys.exit(f"Error: {filename} is empty")
    
    return sequence

def cut_sites(seqeunce, enzyme, enzyme_data):
    #This function gets a sequence, an enzyme and an enzyme dataframe and returns a dataframe of 
    #all the cut indices and the Enzyme

    enzyme_index = enzyme_data[enzyme_data["Enzyme name"] == enzyme].index[0]
    recognition_sequence = enzyme_data.at[enzyme_index, "Recognition sequence"]
    cut_index = enzyme_data.at[enzyme_index, "Cut Index"]

    cut_sites = []

    match_sequence = re.compile(recognition_sequence)

    for match in match_sequence.finditer(seqeunce):
        cut_pos = match.start() + cut_index
        cut_sites.append({"Enzyme name": enzyme, "Cut position": cut_pos})

    enzyme_cut_sites = pd.DataFrame(cut_sites)

    return enzyme_cut_sites


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

resolution = 0.1
upper_limit_electrophoresis = 11000
lower_limit_electrophoresis = 80

DNA_ladder = [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1200, 1000, 900, 800, 600, 500, 400, 300, 200, 100]

##Data
path_genome = "Support\\ECPlambda.fasta"

##Script

#Print all dataframes as test
"""
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
"""

#sys.stdout = open('output.txt', 'w')

jobs = job_queue(job_definition, enzyme_data, adaptor_data, primer_data)
#with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', None):
#    print(jobs)

#Get path to the genome
path_genome = "Support\\ECPlambda.fasta"

absolute_path = os.path.dirname(__file__)
relative_path_genome = path_genome
full_path_genome = os.path.join(absolute_path, relative_path_genome)

sequence = sequence_fasta(full_path_genome)

#Create folder for reports
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
parent_folder = Path(f"Report_{timestamp}")
plots_folder = parent_folder / "plots"
plots_folder.mkdir(parents = True, exist_ok = True)


#create the fingerprints
for index, job_entry in jobs.iterrows():
    #Create a dataframe with all cut sites in the sequence.
    #Gather the enzymes.
    enzyme_name_list = [job_entry["Name Enzyme 1"], job_entry["Name Enzyme 2"], job_entry["Name Enzyme 3"]]
    #print(enzyme_name_list)
    cut_sites_job = pd.DataFrame(
            {"Enzyme name"      :[],
             "Cut position"        :[]}
        )
    #Create dattaframe of cut sites
    for enzyme_entry in enzyme_name_list:
        if enzyme_entry == "none":
            continue
        else:
            new_cut_sites = cut_sites(sequence, enzyme_entry, enzyme_data)
            #print(new_cut_sites)
            cut_sites_job = pd.concat([cut_sites_job, new_cut_sites], ignore_index=True)
    cut_sites_job = cut_sites_job.sort_values("Cut position", ascending=True).reset_index(drop=True)
    #print(cut_sites_job)
    #Now create fragments dataframe.
    fragments = pd.DataFrame(
        {"Start position"       :[],
         "End position"         :[],
         "Start enzyme"         :[],
         "End enzyme"           :[],
         "5' Adaptor sequence"  :[],
         "3' Adaptor sequence"  :[],
         "Fragment sequence"    :[],
         "Length (wA)"          :[],
         "Upper length"         :[],
         "Lower length"         :[],
         "Overlap"              :[]}
    )
    for i in range(len(cut_sites_job)-1):
        #Get fragment sequence
        start_position = int(cut_sites_job.at[i, "Cut position"])
        end_position = int(cut_sites_job.at[i+1, "Cut position"])
        #print(f"Start position: {start_position}, End position: {end_position}")
        fragment_sequence = sequence[start_position:end_position]
        #Check if fragment will be amplified
        start_enzyme = cut_sites_job.at[i, "Enzyme name"]
        end_enzyme = cut_sites_job.at[i+1, "Enzyme name"]
        #Skip fragments which start or end with Enyzme 3 digestion
        enzyme_3 = jobs.at[index, "Name Enzyme 3"]
        if start_enzyme == enzyme_3 or end_enzyme == enzyme_3:
            continue
        else:
            #Get the start and stop sequences
            enzyme_1 = jobs.at[index, "P1: Enzyme"]
            enzyme_2 = jobs.at[index, "P2: Enzyme"]
            print("Enzyme 1 and Enzyme 2:")
            print(enzyme_1, type(enzyme_1), enzyme_2, type(enzyme_2))
            if start_enzyme == enzyme_1:
                start_match = jobs.at[index, "P1: 5' Recognition Sequences"]
            elif start_enzyme == enzyme_2:
                start_match = jobs.at[index, "P2: 5' Recognition Sequences"]
            if end_enzyme == enzyme_1:
                end_match = jobs.at[index, "P1: 3' Recognition Sequences"]
            elif start_enzyme == enzyme_2:
                end_match = jobs.at[index, "P2: 3' Recognition Sequences"]
            if fragment_sequence.startswith(start_match) and fragment_sequence.endswith(end_match):
                #Find extentions
                if start_enzyme == enzyme_1:
                    start_extension = jobs.at[index, "P1: 5' Extention Sequence"]
                elif start_enzyme == enzyme_2:
                    start_extension = jobs.at[index, "P2: 5' Extention Sequence"]
                if end_enzyme == enzyme_1:
                    end_extension = jobs.at[index, "P1: 3' Extention Sequence"]
                elif start_enzyme == enzyme_2:
                    end_extension = jobs.at[index, "P2: 3' Extention Sequence"]
                #Calculate length of fragment
                fragment_length = len(start_extension) + len(fragment_sequence) + len(end_extension)
                upper_length = fragment_length + resolution * fragment_length
                lower_length = fragment_length - resolution * fragment_length
                new_fragment = pd.DataFrame(
                    {"Start position"      :[start_position],
                    "End position"         :[end_position],
                    "Start enzyme"         :[start_enzyme],
                    "End enzyme"           :[end_enzyme],
                    "5' Adaptor sequence"  :[start_extension],
                    "3' Adaptor sequence"  :[end_extension],
                    "Fragment sequence"    :[fragment_sequence],
                    "Length (wA)"          :[fragment_length],
                    "Upper length"         :[upper_length],
                    "Lower length"         :[lower_length],
                    "Overlap"              :[0]}
                    )
                print(new_fragment)
                fragments = pd.concat([fragments, new_fragment], ignore_index=True)
    #Remove too large and too small fragments:
    fragments = fragments[fragments["Length (wA)"].between(lower_limit_electrophoresis, upper_limit_electrophoresis)]
    #look and mark the overlaps
    for i in range(len(fragments)-1):
        i_upper_length = fragments.at[i, "Upper length"]
        i1_lower_length = fragments.at[i, "Lower length"]
        if i_upper_length >= i1_lower_length:
            fragments.at[i, "Overlap"] = 1
            fragments.at[i+1, "Overlap"] = 1
    #Enter number of overlaps into job dataframe
    no_overlaps = fragments["Overlap"].sum()
    jobs.at[i, "Number of Overlaps"] = no_overlaps

    #Create graphic





    

