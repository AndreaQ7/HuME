#!/usr/bin/bash


# This script trims and filters the reads and runs spades.py on all of them in order to generate 
# a contig file for each sample.

# Create result folders
mkdir results
mkdir results/trimmed
mkdir results/metaspadesOut

# debugging command
#set -uex
# Set variables for the file calling inside while loop
var1="_R1_001.fastq.gz"
var2="_R2_001.fastq.gz"
rawreads="/path/to/raw_reads"
list_samples="/path/to/basename_list.txt"


# Enter the while loop and run AdapterRemoval with no --collapse flag. It will generate 2 trimmed files (forward and reverse) for each sample
source activate preprocessing-2020.3
while read -r i; do
  # Create a folder for each sample
  mkdir results/trimmed/$i
  echo $i$var1
  echo $i$var2
  AdapterRemoval --basename trimmed/$i/$i --file1 $rawreads/$i$var1 --file2 $rawreads/$i$var2
done < $list_samples


# Assign trimmed folder to a path variable
path="trimmed"

# Run spades.py in --meta mode (metaspades) in order to get the create assemble genomes for each sample
for i in $(ls -A $path); do
  cat $path/$i$var2 | head -2
  spades.py --meta -1 $path/$i/$i$var1 -2 $path/$i/$i$var2 -o results/metaspadesOut/$i
done
