# InSilicoAFLP
My attempt at creating a Program which allows simulating AFLP fingerprints from a sequenced genome, and, eventually, automatically optimize  the parameters

/Config
Contains configuration data for the script:
primers.yaml: Data on primers: Name, associated restriction enzyme, core sequence, selective bases
enzymes.yaml: Data on restriction enzymes. Name, recognition seqeunce with cut position
adaptors.yaml: Data on adaptors: Name, associated restriction enzyme, 5'-3' sequence, 3-5' sequence
electrophoresis.yaml: Data on electrophoretic system, name, resolution, upper cut-off, lower cut-off, ladder
jobs.yaml: Data on fingerprint sot be generated. Method, restriction enzymes, primers, adaptors, electrophoretic system.

template: Template used to generate the report

/Data
Contains seqeunce data from which fingerprint will be generated. 
