#!/usr/bin/bash

# This script Normalize the counts made by humann2 for each sample.

# Set debugging mode
set -uex

# Set variable for list of samples (change this to user defined directory)
sample_list="/path/to/sample_list"

# Enter the loop and normalize each count for pathcoverage, pathabundance and genefamilies files
while read l; do
  echo "normalizing $l"
  humann2_renorm_table --input Human2Out/$l/${l}_pathcoverage.tsv --output NormalizedPathCoverage/${l}_pathcoverage_norm.tsv --units cpm
  humann2_renorm_table --input Human2Out/$l/${l}_pathabundance.tsv --output NormalizedPathAbundances/${l}_pathabundance_norm.tsv --units cpm
  humann2_renorm_table --input Human2Out/$l/${l}_genefamilies.tsv --output NormalizedGeneFamilies/${l}_genefamilies_norm.tsv --units cpm
done < $sample_list
