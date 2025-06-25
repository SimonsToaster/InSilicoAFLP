#1. Variables

recognition_site = "CTGCAG"
cut_index = 5

adaptor_seqence = "CTCGTAGACTGCGTACATGCA"

pcr_sequence = "TACTCAGGACTC"
select_bases = "G"

path_genome = "Support\\ECPlambda.fasta"

#Thermo Scientific GeneRuler DNA Ladder Mix
DNA_ladder = [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1200, 1000, 900, 800, 600, 500, 400, 300, 200, 100]

#path_resultfile = "Results\savefile"

#2. Code
#2.1 Imports

import os
import re
import matplotlib.pyplot as plt
import numpy as np

#2.2 Program 

#2.2.1 Import genome from file

absolute_path = os.path.dirname(__file__)
relative_path_genome = path_genome
full_path_genome = os.path.join(absolute_path, relative_path_genome)

with open(full_path_genome, "r") as file:
    next(file)
    genome_raw = file.read()
    genome = genome_raw.replace('\n', '').upper()

#2.2.2 Find length of amplified fragments
#Find position of recognition sites in genome

index_recognition_sites = [m.start() for m in re.finditer(recognition_site, genome)]

print(index_recognition_sites[1:10])
print("Amount of recognition site indices", len(index_recognition_sites))

#Correct for site of restriction

index_restriction_site = []
for index in index_recognition_sites:
    index_restriction_site.append(index+cut_index)

print(index_restriction_site[1:10])
print("Amount of corrected recognition site indices", len(index_restriction_site))

#Make fragments from restriction site list
#I'm a bit concerned whether this will work rerectly wit respect to the last fragment.

restriction_fragments = []
i = 0
j = 1
while j < len(index_restriction_site):
    restriction_fragments.append(genome[index_restriction_site[i]:index_restriction_site[j]])
    i += 1
    j += 1

print(restriction_fragments[1])
print("Number of fragments by restriction:", len(restriction_fragments))

#Select fragments which are amplified

selected_fragments = []
selector_sequence_5prime = recognition_site[cut_index:len(recognition_site)] + select_bases
print("Selector sequence at 5 end of sequence", selector_sequence_5prime)

complement_select_bases = select_bases.translate(str.maketrans("ATGC", "TACG"))
print("Complement_sequence", complement_select_bases)

selector_sequence_3prime = complement_select_bases + recognition_site[0:cut_index]
print("Selector sequence at 3 end of sequence", selector_sequence_3prime)

for fragment in restriction_fragments:
    if str(fragment).startswith(selector_sequence_5prime) and str(fragment).endswith(selector_sequence_3prime):
        selected_fragments.append(fragment)

#print(selected_fragments[1])
print("number of amplicons without size selection:", len(selected_fragments))


#I should append adaptors but im lazy and want to now reduce it to lengths in numbers
#In big genomes even with three selective bases thousands of fragments are created. 
#According to the Internet Taq struggles to amplify fragments of more than 5kB. Everythign above 6 will thus be removed
#Agarose gels have trouble resolving below a certain size. Atm everything below 100bp is cut off

lenghts_amplicons = []
for fragment in selected_fragments:
    potential_amplicon = len(fragment)+2*len(pcr_sequence)
    if potential_amplicon <= 5000 and potential_amplicon >= 100:
        lenghts_amplicons.append(potential_amplicon)

print(lenghts_amplicons)
print("Number of amplicons after size selection:", len(lenghts_amplicons))
print("Longest amplicon:", max(lenghts_amplicons))
print("Shortest amplicon:", min(lenghts_amplicons))

#2.3 Creating the plot
# Draw ladder
lane_ladder = []
i=0
while i < len(DNA_ladder):
    lane_ladder.append(1)
    i += 1

x = np.array(lane_ladder)
y = np.array(DNA_ladder)

plt.scatter(x,y, marker = 's')

for i, txt in enumerate(DNA_ladder):
    plt.annotate(txt, (x[i], y[i]), va ="center", xytext=(1.1,y[i]))

#Draw fingerprint
lane_fingerprint = []
i=0
while i < len(lenghts_amplicons):
    lane_fingerprint.append(2)
    i += 1

x = np.array(lane_fingerprint)
y = np.array(lenghts_amplicons)

plt.scatter(x,y, marker = 's')
for i, txt in enumerate(lenghts_amplicons):
    plt.annotate(txt, (x[i], y[i]), va ="center", xytext=(2.1,y[i]))


plt.xticks(np.arange(0, 5, 1))
plt.yscale("log")

plt.show()