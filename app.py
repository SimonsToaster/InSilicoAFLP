#1. Variables

recognition_site = "TTAA"
cut_index = 1

adaptor_seqence = "GACGATGAGTCCTGAG"

pcr_sequence = "TACTCAGGACTC"
select_bases = "A"

path_genome = "Support\\ECPlambda.fasta"

#path_resultfile = "Results\savefile"

#2. Code
#2.1 Imports

import os
import re

#2.2 Program 

#2.2.1 Import genome from file

absolute_path = os.path.dirname(__file__)
relative_path_genome = path_genome
full_path_genome = os.path.join(absolute_path, relative_path_genome)

with open(full_path_genome, "r") as file:
    next(file)
    genome_raw = file.read()
    genome = genome_raw.replace('\n', '')

#2.2.2 Find leght of amplified fragments
#Find position of recognition sites in genome

index_recognition_sites = [m.start() for m in re.finditer(recognition_site, genome)]

print(index_recognition_sites)
print(len(index_recognition_sites))

#Correct for site of restriction

index_restriction_site = []
for index in index_recognition_sites:
    index_restriction_site.append(index+cut_index)

print(index_restriction_site)
print(len(index_restriction_site))

#Make fragments from restriction site list
#I'm a bit concerned whether this will work rerectly wit respect to the last fragment.

restriction_fragments = []
i = 0
j = 1
while j < len(index_restriction_site):
    restriction_fragments.append(genome[index_restriction_site[i]:index_restriction_site[j]])
    i += 1
    j += 1

print(restriction_fragments[1:10])
print(len(restriction_fragments))

#Select fragments which are amplified

selected_fragments = []
selector_sequence = recognition_site[cut_index:len(recognition_site)] + select_bases
print(selector_sequence)
for fragment in restriction_fragments:
    if fragment[0:len(selector_sequence)] == selector_sequence:
        selected_fragments.append(fragment)

print(selected_fragments[1:10])
print(len(selected_fragments))

#I should append adaptors but im lazy and want to now reduce it to lengths in numbers

lenghts_amplicons = []
for fragment in selected_fragments:
    lenghts_amplicons.append(len(fragment)+2*len(pcr_sequence))

print(lenghts_amplicons)
print(len(lenghts_amplicons))
