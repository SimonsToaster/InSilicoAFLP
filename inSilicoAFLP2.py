#----------------------------------------------------------------------------------------------------------------------------------------
#IMPORTS
#----------------------------------------------------------------------------------------------------------------------------------------

import sys
import os
import re
import itertools as it

#----------------------------------------------------------------------------------------------------------------------------------------
#FUNCTIONS
#----------------------------------------------------------------------------------------------------------------------------------------

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

def expand_selective_bases(sequence):
    #This function takes a seqeunce of nuclotides according to IUPAC numencalture
    #and returns a list of all seqeunces they represent. 

    IUPAC_decoder = {'A': ['A'],'C': ['C'],'G': ['G'],'T': ['T'],'U': ['T'],'R': ['A', 'G'],'Y': ['C', 'T'],'S': ['G', 'C'],'W': ['A', 'T'],'K': ['G', 'T'],'M': ['A', 'C'],'B': ['C', 'G', 'T'],'D': ['A', 'G', 'T'],'H': ['A', 'C', 'T'],'V': ['A', 'C', 'G'],'N': ['A', 'C', 'G', 'T']}
    options = [IUPAC_decoder[base.upper()] for base in sequence]
    sequence_list = [''.join(p) for p in it.product(*options)]

    return sequence_list

def five_prime_end(recognition_seqeunce, cut_index, selective_bases):
    #This function takes a recognition seqeunce and cut index of a Restriction enyzme and selective bases and creates the matching seqeunces for the 5' end
    #of fragments

    rec_site_fragment = recognition_seqeunce[cut_index:len(recognition_seqeunce)]
    five_prime_end = rec_site_fragment + selective_bases

    return five_prime_end

def three_prime_end(recognition_seqeunce, cut_index, selective_bases):
    #This function takes a recognition sequence and cut index of a restriction enyzme and selective bases and creates the matching seqeunces for the 3' end
    #of fragments
    
    rec_site_fragment = recognition_seqeunce[0:cut_index]
    complement_selective_bases = selective_bases.translate(str.maketrans("ATGC", "TACG"))
    three_prime_end = complement_selective_bases + rec_site_fragment

    return three_prime_end

def match_dictionary(Name_RE, recognition_sequence, cut_index, selective_bases):
    #This function takes a recognitin sequenc and cut site of a restriction enzyme and selective bases and returns a dictionary with the selective bases as keys and the 5' and 3' matching ends as values
    
    selective_bases_list = expand_selective_bases(selective_bases)

    match_dictionary = {(Name_RE, entry): [] for entry in selective_bases_list}

    for key in match_dictionary:
        five_match = five_prime_end(recognition_sequence, cut_index, key[1])
        three_match = three_prime_end(recognition_sequence, cut_index, key[1])
        match_dictionary[key].append(five_match)
        match_dictionary[key].append(three_match)

    return match_dictionary

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

def job_queue(jobs, enzyme_data, adaptor_data, primer_data):
    #This function takes the job list, and the dictionaries with enzyme, primer and adaptor data and creates a job list. Each job 
    #is equivalent to a new finger print

    job_list = []

    #Turns the unformated job list into a list of lists so the information on enzymes, adaptors and primers in each reaction can be
    #retrieved from the dictionaries.
    for job_entry in jobs:
        method = job_entry[0]
        digestion_mix = job_entry[1].split("+")
        adaptor_mix = job_entry[2].split("+")
        primer_mix = job_entry[3].split("+")
        job = [method, digestion_mix, adaptor_mix, primer_mix]
        job_list.append(job)

    #Replaces the names of enzymes, adaptors and primers with information from the dictionary
    for job_entry in job_list:
        digestion_mix = job_entry[1]
        for i in range(len(digestion_mix)):
            enzyme = digestion_mix[i]
            digestion_mix[i] = enzyme_data[enzyme]
        adaptor_mix = job_entry[2]
        for i in range(len(adaptor_mix)):
            adaptor = adaptor_mix[i]
            adaptor_mix[i] = adaptor_data[adaptor]
        primer_mix = job_entry[3]
        for i in range(len(primer_mix)):
            primer = primer_mix[i]
            primer_mix[i] = primer_data[primer]
    
    #Degenerated selective bases are decoded
    #for job_entry in job_list:
    #    primer_list = job_entry[2]
    #    print(primer_list)
    #    for primer_entry in primer_list:
    #        if len(primer_entry) == 4:
    #            decoded_selective_bases = expand_selective_bases(primer_entry[3])
    #            primer_entry.append(decoded_selective_bases)
    
    #Create 5', 3' match sequences for every decoded selective sequence
    for job_entry in job_list:
        primer_list = job_entry[3]
        for primer_entry in primer_list:
            if len(primer_entry) == 4:
                enzyme = primer_entry[1]
                enzyme_info = enzyme_data[enzyme]
                match_seqeunce_dictionary = match_dictionary(enzyme_info[0], enzyme_info[1], enzyme_info[2], primer_entry[3])
                primer_entry.append(match_seqeunce_dictionary)
    
    #create 5' and 3' extension for the amplicons

    for job_entry in job_list:
        primer_list = job_entry[3]
        for primer_entry in primer_list:
            if len(primer_entry) == 5:
                enzyme = primer_entry[1]
                enzyme_info = enzyme_data[enzyme]
                five_extension = five_primer_extention(primer_entry[2], enzyme_info[2])
                three_extension = three_primer_extention(primer_entry[2], enzyme_info[2])
                primer_entry.append(five_extension)
                primer_entry.append(three_extension)
    
    job_dict = {}

    for job_entry in job_list:
        method = job_entry[0]
        RE_info = job_entry[1]
        adaptor_info = job_entry[2]
        primer_info = job_entry[3]
        #gather the selctive bases together
        primer_list = job_entry[3]
        complete_selective_bases_dict = {}
        for primer_entry in primer_list:
            complete_selective_bases_dict.update(primer_entry[4])
        print("Complete selective bases dict")
        print(complete_selective_bases_dict)
        if method == "single":
            for entry in complete_selective_bases_dict:
                job = {}
                primer_info_added = primer_info[0][0:3]
                selctive_bases_info = complete_selective_bases_dict[entry]
                primer_info_added.append([selctive_bases_info])
                primer_info_added.append(primer_info[0][5:7])
                job_name = (method, entry[0], entry[1])
                job = {job_name: [RE_info, adaptor_info, primer_info_added]}
                job_dict.update(job)
        else:
            #All combinations need to be made first.
            all_primers = list(complete_selective_bases_dict.keys())
            print("All primers")
            print(all_primers)
            all_primer_pairs = list(it.combinations(all_primers, 2))
            print("All primer pairs")
            print(all_primer_pairs)
            filtered_combinations = []
            for entry in all_primer_pairs:
                if entry[0][0] != entry[1][0]:
                    filtered_combinations.append(entry)
            print("filtered combinations")
            #Now create jobs for every combination
            for entry in filtered_combinations:
                job = {}
                primer_info_added = primer_info[0][0:3]
                selctive_bases_info = []
                matches_enzyme_1 = complete_selective_bases_dict[entry[0]]
                matches_enzyme_2 = complete_selective_bases_dict[entry[1]]
                mixed_fragment_1 = [matches_enzyme_1[0], matches_enzyme_2[1]]
                mixed_fragment_2 = [matches_enzyme_2[0], matches_enzyme_1[1]]
                selctive_bases_info = [matches_enzyme_1, matches_enzyme_2, mixed_fragment_1, mixed_fragment_2]
                primer_info_added.append(selctive_bases_info)
                primer_info_added.append(primer_info[0][5:7])
                job_name = (method, entry[0], entry[1])
                job = {job_name: [RE_info, adaptor_info, primer_info_added]}
                job_dict.update(job)

    return job_dict

def digest(sequence_list, enzyme_entry):
    #Function takes a library of sequences, and a library of a restriction enzymes recognition site and cutting index
    #and returns a library of fragments which would result from this digestion.

    fragment_list = []

    for sequence in sequence_list:
        last_end_index = 0
        for match in re.finditer(enzyme_entry[1], sequence):
            split_index = match.start() + enzyme_entry[2]
            fragment_list.append(sequence[last_end_index:split_index])
            last_end_index = split_index
        fragment_list.append(sequence[last_end_index:])
    
    return fragment_list

def selection(sequence_list, selection_list):
    #Function takes a list of fragments and a list of an beginning and end sequence and filters fragments which begin an end with these seqeunces.
    #IMPORTANT: This means the selection list needs to consist of selective bases and the restriction overhang.
    #first entry is at 5' end, second entry is at 3'end

    filtered_sequence_list = []

    for sequence in sequence_list:
        if str(sequence).startswith(selection_list[0]) and str(sequence).endswith(selection_list[1]):
            filtered_sequence_list.append(sequence)
    
    return filtered_sequence_list

def extend_fragments(fragment_list, extension_list):
    extended_fragments_list = []
    return extended_fragments_list

#----------------------------------------------------------------------------------------------------------------------------------------
#SCRIPT
#----------------------------------------------------------------------------------------------------------------------------------------

##Data
path_genome = "Support\\ECPlambda.fasta"

enzyme_data = {"PstI": ["PstI", "CTGCAG", 5],
               "EcoRI": ["EcoRI", "GAATTC", 1],
               "MseI": ["MseI", "TTAA", 1]}

adaptor_data = {"PstI Adaptor 1": ["PstI Adaptor 1", "PstI", "CTCGTAGACTGCGTACATGCA", "CATCTGACGCATGT"],
                "EcoRI Adaptor 1": ["EcoRI Adaptor 1", "EcoRI", "CTCGTAGACTGCGTACC", "AATTGGTACGCAGTCTAC"],
                "MseI Adaptor 1":["MseI Adaptor 1", "MseI", "GACGATGAGTCCTGAG", "TACTCAGGACTCAT"]}

primer_data = {"PstI Primer 1": ["PstI Primer 1", "PstI", "GACTGCGTACATGCAG", "R"],
               "EcoRI Primer 1": ["EcoRI Primer 1", "EcoRI", "GACTGCGTACCAATTC", "Y"],
               "Mse1 Primer 1": ["Mse1 Primer 1", "MseI", "GATGAGTCCTGAGTAA", "ATG"]}

jobs = [["single", "PstI", "PstI Adaptor 1", "PstI Primer 1"],
        ["double", "PstI+EcoRI", "PstI Adaptor 1+EcoRI Adaptor 1", "PstI Primer 1+EcoRI Primer 1"],
        ["triple", "PstI+EcoRI+MseI", "PstI Adaptor 1+EcoRI Adaptor 1", "PstI Primer 1+EcoRI Primer 1"]]

##Script

job_dictionary = job_queue(jobs, enzyme_data, adaptor_data, primer_data)
print(job_dictionary)

#print(job_dictionary)

absolute_path = os.path.dirname(__file__)
relative_path_genome = path_genome
full_path_genome = os.path.join(absolute_path, relative_path_genome)

sequence = [sequence_fasta(full_path_genome)]

results_dict = {}

for job_entry in job_dictionary:
    job_data = job_dictionary[job_entry]
    enzyme_list = job_data[0]
    fragment_library = sequence
    for enzyme in enzyme_list:
        fragment_library = digest(fragment_library, enzyme)
    #print(f"Number of restrcition fragments{len(fragment_library)}")
    match_list = job_data[2][3]
    filtered_fragments = []
    for match_pair in match_list:
        print(match_pair)
        filtered = selection(fragment_library, match_pair)
        filtered_fragments.append(filtered)
    #print(f"Number of filtered fragments: {len(filtered_fragments)}")
    #attach extension to fragments

