#!/usr/bin/bash


# Once determined the good quality bins with CheckM, kraken2 will be rerun on these assemblies to assess if they correspond
# to a given species

# Create Good_quality_Bins folder and copy these assembly files according to CheckM estimation

# FOR USERS: if necessary, activate the environment in which you have kraken2 binaries

# Run kraken on all the files in Good_quality_Bins
mkdir kraken
mkdir krakenReport
pathkrakendb="/path/to/krakenDB"
for i in Good_quality_Bins/*; do
  echo "$i"
  file=$(basename $i)
  echo $file
  for bin in $i/*; do
    echo $bin
    namebin=$(basename $bin)
    mkdir kraken/$file
    mkdir krakenReport/$file
    kraken2 --db $pathkrakendb --threads 10 --output kraken/$file/$namebin.kraken --report krakenReport/$file/$namebin.kraken.report $bin
  done
done

