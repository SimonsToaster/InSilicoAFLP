#----------------------------------------------------------------------------------------------------------------------------------------
# This script formats a .txt file of neb REs into something i can use in yaml files
#----------------------------------------------------------------------------------------------------------------------------------------

import re

pattern_four = re.compile(r'^[ATGC]/[ATGC]{3}$|^[ATGC]{2}/[ATGC]{2}$|^[ATGC]{3}/[ATGC]$')
pattern_six = re.compile(r'^[ATGC]/[ATGC]{5}$|^[ATGC]{2}/[ATGC]{4}$|^[ATGC]{3}/[ATGC]{3}$|^[ATGC]{4}/[ATGC]{2}$|^[ATGC]{5}/[ATGC]$')

with open("functions\helper scripts\RE_cat.txt", "r") as infile, open("functions\helper scripts\\fourBP.txt", "a") as outfile:
    for line in infile:
        content = line.split()
        name = content[1].strip()
        sequence = content[0].strip()
        if pattern_four.match(sequence):
            four_string = f'''
- enzyme_name: {name}
  recognition_sequence: {sequence}
'''
            outfile.write(four_string)

with open("functions\helper scripts\RE_cat.txt", "r") as infile, open("functions\helper scripts\\sixBP.txt", "a") as outfile:
    for line in infile:
        content = line.split()
        name = content[1].strip()
        sequence = content[0].strip()
        if pattern_six.match(sequence):
            six_string = f'''
- enzyme_name: {name}
  recognition_sequence: {sequence}
'''
            outfile.write(six_string)