#!/usr/bin bash

# Script for running kraken over the trimmend and collapsed read files.

# Variable for list of samples txt file
list_files="/path/to/list_of_files.txt"

while read line; do
  echo "Parsing $line"
  kraken2 --db ../KrakenDB /mnt/hume/TrimmedDentalCalculus/$line/$line.collapsed --threads 15 --output krakenOut/$line.kraken --report krakenReport/$line.kraken.report
done < $list_files
