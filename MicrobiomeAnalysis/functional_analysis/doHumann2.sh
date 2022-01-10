#!/usr/bin/bash

# Set debugging mode
set -uex

# Set variable for list of samples (change this to user defined directory)
sample_list="/path/to/sample_list"

while read file; do
  mkdir Human2Out/$file
  echo "$file"
  humann2 -i /mnt/hume/TrimmedDentalCalculus/$file/$file.collapsed.gz --metaphlan MetaPhlan2/metaphlan2.py -o Human2Out/$file/ --nucleotide-database ../ChocophlanDB/chocophlan/ --protein-database ../UniRef90EC/uniref/ --metaphlan-options '--mpa_pkl /mnt/hume/Gabriel/bowtie2DB/mpa_v20_m200/mpa_v20_m200.pkl --bowtie2db /mnt/hume/Gabriel/bowtie2DB/mpa_v20_m200/mpa_v20_m200' --nucleotide-database ../ChocophlanDB/chocophlan/ --protein-database ../UniRef90EC/uniref/ --threads 20 --pathways-database /home/ginnocenti/humann_legacy-0.99b/data/keggc
done < ${sample_list}


