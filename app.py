#1. Variables

recognition_site = "TTAA"
cut_index = 1

adaptor_seqence = "GACGATGAGTCCTGAG"

pcr_sequence = "TACTCAGGACTC"
select_bases = "A"

path_genome = 'Support\ECPlambda.fasta'
path_resultfile = "Results\savefile"

#2. Code
#2.1 Imports

import os
import re 

#2.2 Program 

with open(path_genome, "r") as file:
    next(file)
    genome = file.read()

indices = re.search(recognition_site, genome)
print (indices)

#with open(path_resultfile, "w") as file:
 #   file.write(genome)



