#!/usr/bin/bash

# Set debugging mode
set -uex

# Create list variable with all samples
list_files="/path/to/list_of_samples"


while read l; do
  kreport2mpa.py -r krakenReport/$l.kraken_bracken_species.report -o kraken2mpa/$l.txt --display-header
done < $list_files
