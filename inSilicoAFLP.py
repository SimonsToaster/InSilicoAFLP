
import sys
import os
import re
import itertools as it


#### FUNCTIONS FUNCTIONS FUNCTIONS FUNCTIONS ####

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

def digest(sequence_list, enzyme_entry):
    #Function takes a library of sequences, and a library of a restriction enzymes recognition site and cutting index
    #and returns a library of fragments which would result from this digestion.

    fragment_list = []

    for seqequence in sequence_list:
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

def generate_triple_digests(enzyme_list):
    #This function gets a list of enzymes and uses it to generate all possible triple digests: 
    #All unique combinations of pairs, then adding a third which isnt already in the pair.

    pairs = list(it.combinations(enzyme_list, 2))

    triplets = []

    for pair in pairs:
        pair_list = list(pair)
        for enzyme in enzyme_list:
            if enzyme != pair[0] and enzyme != pair[1]:
                triplet = pair_list + [enzyme]
                triplets.append(triplet)

    return triplets

def convert_recognition_sequence(recognition_cut_sequence):
    #Function takes a recognition seqeunce with indicated cut position and returns 
    #a clean recognition seqeunce and the cut indices. 

    clean_recognition_sequence = recognition_cut_sequence.replace("/", "")
    
    match = re.search("/", recognition_cut_sequence)
    cut_index = match.start()

    return [clean_recognition_sequence, cut_index]

def expand_selective_bases(sequence):
    #This function takes a seqeunce of nuclotides according to IUPAC numencalture
    #and returns a list of all seqeunces they represent. 

    IUPAC_decoder = {'A': ['A'],'C': ['C'],'G': ['G'],'T': ['T'],'U': ['T'],'R': ['A', 'G'],'Y': ['C', 'T'],'S': ['G', 'C'],'W': ['A', 'T'],'K': ['G', 'T'],'M': ['A', 'C'],'B': ['C', 'G', 'T'],'D': ['A', 'G', 'T'],'H': ['A', 'C', 'T'],'V': ['A', 'C', 'G'],'N': ['A', 'C', 'G', 'T']}
    options = [IUPAC_decoder[base.upper()] for base in sequence]
    sequence_list = [''.join(p) for p in it.product(*options)]

    return sequence_list

def five_prime_ends(recognition_seqeunce, cut_index, selective_bases):
    #This function takes a recognition seqeunce and cut index of a Restriction enyzme and a list of selective bases and creates the matching seqeunces for the 5' end
    #of fragments

    five_prime_ends_list = []

    rec_site_fragment = recognition_seqeunce[cut_index:len(recognition_seqeunce)]

    if type(selective_bases) == str:
        selective_bases_list = [selective_bases]

    for sequence in selective_bases_list:
        five_prime_end = rec_site_fragment + sequence
        five_prime_ends_list.append(five_prime_end)

    return five_prime_ends_list

def three_prime_ends(recognition_seqeunce, cut_index, selective_bases):
    #This function takes a recognition sequence and cut index of a restriction enyzme and a list of selective bases and creates the matching seqeunces for the 3' end
    #of fragments

    three_prime_ends_list = []

    rec_site_fragment = recognition_seqeunce[0:cut_index]

    if type(selective_bases) == str:
        selective_bases_list = [selective_bases]

    for sequence in selective_bases_list:
        complement_selective_bases = sequence.translate(str.maketrans("ATGC", "TACG"))
        three_prime_end = complement_selective_bases + rec_site_fragment
        three_prime_ends_list.append(three_prime_end)

    return three_prime_ends_list

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

def match_dictionary(recognition_sequence, cut_index, selective_bases):
    #This function takes a recognitin sequenc and cut site of a restriction enzyme and selective bases and returns a dictionary with the selective bases as keys and the 5' and 3' matching ends as values
    
    selective_bases_list = expand_selective_bases(selective_bases)

    match_dictionary = {entry: [] for entry in selective_bases_list}

    for key in match_dictionary:
        five_match = five_prime_end(recognition_sequence, cut_index, key)
        three_match = three_prime_end(recognition_sequence, cut_index, key)
        match_dictionary[key].append(five_match)
        match_dictionary[key].append(three_match)

    return match_dictionary

def five_primer_extention(primer_seqeunce, cut_index):
    #This fucntion takes primer and recognition sequence and cut index to generate the 5' and 3' extention of the amplicon

    five_prime_extension = primer_seqeunce[0:(len(primer_seqeunce)-cut_index)]

    return five_prime_extension

def three_primer_extention(primer_seqeunce, cut_index):
    #This function takes primer and recognition sequence and cut index to generate the 5' and 3' extention of the amplicon

    reversed_primer = primer_seqeunce[::-1]
    reversed_complement = reversed_primer.translate(str.maketrans("ATGC", "TACG"))
    three_primer_extension = reversed_complement[cut_index:]

    return three_primer_extension

def create_enzyme_dictionary(enzyme_instruction_data):
    #This function takes the enzyme/primer data specified by the user in the instruction file and creates the enzyme dictionary with all required information for
    #creation of the job queue and calcualtion of the findgerprints

    #Create dictionary with enzyme names as keys
    enzyme_instruction_dictionary = {entry[0]: entry[1:] for entry in enzyme_instruction_data}

    for enzyme, data in enzyme_instruction_dictionary.items():
        #Convert recognition seqeunce into clear seqeunce and cut index
        sequence_and_cutindex = convert_recognition_sequence(data[0])
        data[0] = sequence_and_cutindex[0]
        data.insert(1, sequence_and_cutindex[1])

        #Expand the degenerated coding of selective bases into propper sequences and create
        #5' and 3' match seqeunces, stored in dictionary:
        match_dict = match_dictionary(data[0], data[1], data[-1])
        data[-1] = match_dict

        #Create the 5' and 3' primer extensions
        five_extension = five_primer_extention(data[6], data[1])
        data.insert(7, five_extension)
        three_extension = three_primer_extention(data[6], data[1])
        data.insert(8, three_extension)

    return enzyme_instruction_dictionary


def job_queue(enzyme_dictionary, mode):
    #This function takes the enzyme dictionary, the mode and creates a job dictionary. Key is the names of the enzymes in the digest, value is a library with the enzyme data and the matching seqeunces
    
    if mode == "single":
        sys.exit("not defined yet")
    elif mode == "double":
        sys.exit("not defined yet")
    elif mode == "triple":
        #create all digest combinations
        enzymes = enzyme_dictionary.keys()
        triple_digests = generate_triple_digests(enzymes)
        job_dict = {tuple(triplet): triplet for triplet in triple_digests}

        #Create all matching sequence combinations

        
    
    return job_dict

#### Script Script Script ####

#Constants and such

path_genome = "Support\\ECPlambda.fasta"

##GENERAL
mode = "optimize"
best = 20
method = "triple"

##ELECTROPHORESIS
upper_limit = 2000
lower_limit = 200
resolution = 0.1

##SPACE
enzyme_instruction_data = [["PstI", "CTGCA/G", "PstI Adaptor 1", "CTCGTAGACTGCGTACATGCA", "CATCTGACGCATGT",  "PstI Primer 1", "GACTGCGTACATGCAG", "ATY"],
                           ["EcoRI", "G/AATTC", "EcoRI Adaptor 1", "CTCGTAGACTGCGTACC", "AATTGGTACGCAGTCTAC", "EcoRI Primer 1", "GACTGCGTACCAATTC", "ATY"],
                           ["MseI", "T/TAA", "MseI Adaptor 1", "GACGATGAGTCCTGAG", "TACTCAGGACTCAT", "Mse1 Primer 1", "GATGAGTCCTGAGTAA", "ATY"]]

restriction_library = [["PstI", "CTGCAG", 5]]

selection_library = []

mode = "triple"

enyzme_data = {"PstI": ["PstI", "CTGCAG", 5],
               "EcoRI": ["EcoRI", "GAATTC", 1],
               "MseI": ["MseI", "TTAA", 1]}


adaptor_data = {"PstI Adaptor 1": ["PstI Adaptor 1", "PstI", "CTCGTAGACTGCGTACATGCA", "CATCTGACGCATGT"],
                "EcoRI Adaptor 1": ["EcoRI Adaptor 1", "EcoRI", "CTCGTAGACTGCGTACC", "AATTGGTACGCAGTCTAC"],
                "MseI Adaptor 1":["MseI Adaptor 1", "MseI", "GACGATGAGTCCTGAG", "TACTCAGGACTCAT"]}

primer_data = {"PstI Primer 1": ["PstI Primer 1", "PstI", "GACTGCGTACATGCAG", "ATY"],
               "EcoRI Primer 1": ["EcoRI Primer 1", "EcoRI", "GACTGCGTACCAATTC", "ATY"],
               "Mse1 Primer 1": ["Mse1 Primer 1", "MseI", "GATGAGTCCTGAGTAA", "ATY"]}

Jobs = [["PstI", "PstI Adaptor 1", "PstI Primer 1"],
        ["PstI+EcoRI" "PstI Adaptor 1+EcoRI Adaptor 1", "PstI Primer 1+EcoRI Primer 1"],
        ["PstI+EcoRI+MseI" "PstI Adaptor 1+EcoRI Adaptor 1", "PstI Primer 1+EcoRI Primer 1"]]

#The Script

##Read seqeunce
absolute_path = os.path.dirname(__file__)
relative_path_genome = path_genome
full_path_genome = os.path.join(absolute_path, relative_path_genome)

sequence = sequence_fasta(full_path_genome)

fragment_library = [sequence]

##Construct the job queue

enzyme_dictionary = create_enzyme_dictionary(enzyme_instruction_data)

jobs = job_queue(enzyme_dictionary, mode)

print(jobs)

##Work over the search space

for enzyme_combination in restriction_library:
    fragment_library = digest(fragment_library, enzyme_combination)
