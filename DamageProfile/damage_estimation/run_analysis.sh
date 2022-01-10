#/usr/bin/bash

# This wrapper runs:
# 1. get_reference_genomes.sh to download every reference genomic fasta file for each given organism ( here we will use
# genbank assembly database ) get_damage_profile.sh scripts.
# 2. get_damage_profile.sh to have a view on damage profiles of each relevant microrganism in the microbiota community
# 3. makePMDtable.py to have the results reported in a table

bash get_reference_genomes.sh

bash get_damage_profile.sh

python makePMDtable.py > damage_patterns.tsv
