#!/usr/bin/bash

# This script runs Bracken on all samples

# Set debugging mode and create list variable for files
set -uex
list_files="/path/to/list_samples.txt"


# Create a folder for each desired taxonomy level abundance estimation
mkdir brackenOut/SpeciesAbundance
mkdir brackenOut/GenusAbundance
mkdir brackenOut/PhylaAbundance

# Loop across the samples
while read line; do
  echo "estimating abundances with Bracken:"
  echo "file ${line}.kraken.report"
  bracken -d ../KrakenDB -i krakenReport/$line.kraken.report -o brackenOut/SpeciesAbundance/$line.bracken -r 50 -l S -t 10
  #bracken -d ../KrakenDB -i krakenReport/$line.kraken.report -o brackenOut/GenusAbundance/$line.genus.bracken -r 50 -l G -t 10
  #bracken -d ../KrakenDB -i krakenReport/$line.kraken.report -o brackenOut/PhylaAbundance/$line.phylum.bracken -r 50 -l P -t 10
done < $list_files
