#!/usr/bin/bash

# This is the wrapper for the genome assembly workflow.

# Run metaspades
cd final_assembly
bash run_metaspades.sh

# Run remapping
bash remap_contigs.sh

# Run metaBAT
bash run_metabat.sh

cd ..


# FOR USERS: it is suggested to manually analyze the good quality bins and undergo the below kraken2 analysis.
# Then the Good quality genome assemblies will be taken in consideration

# Run kraken2
cd taxonomy_identification
bash runKraken2.sh
cd ..
