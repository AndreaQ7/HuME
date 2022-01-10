#!/usr/bin/bash

# This script uses metaBAT binaries to bin each contig to its predicted species
# Create metaBAT results folder 
mkdir metaBAT

# FOR USER: Please download metabat in your home before running the script


# Create variables with Path:
# 1. To a user defined list with sample names
# 2. To metaspades out
# 3. To bam remapping result  
list_samples="/path/to/list_samples.txt"
pathmetaspades="results/metaspades"
pathbam="results/bamresults"

# Loop through each file in metaspades output
while read -r namesample; do
  namesample=$(basename $namesample)
  echo $namesample
  mkdir metaBAT/$namesample
  cd metaBAT/$namesample
  location=$(pwd)
  echo "Into "${location}
  bash ~/metabat/bin/runMetaBat.sh -t 0 $pathmetaspades/$namesample/contigs.fasta $pathbam/$namesample/$namesample.sorted.bam
  cd ../../
done < $list_files


# Now, to simplify the next steps, filter only the samples which binned file was created, and take only the bin files 
# Create a new folder in which to store results
mkdir results/metaBatOut

# Loop into metaBAT output folder and extract only the created bins.
# Output in -> results/metaBatOut
pathout="results/metaBatOut

for file in metaBAT/*; do
  for f in $file/*; do
    namefile=$(basename $file)
    if [ -d $f ]; then
      echo $f
      namebins=$(ls $f | grep -E '\.fa$' | grep -v $namefile | grep -v 'corrected')
      echo $namebins
      if [ -z "$namebins" ]; then
        :
      else
        ls $f | grep -E '\.fa$' | grep -v $namefile | grep -v 'corrected' > $file/binList.txt
        #cat $file/binList.txt | parallel --will-cite "mv "$f"/{} "$f"/"$namefile".{}"
        mkdir $pathout/$namefile
        mv $f/* $pathout/$namefile/
        # rm $file/binList.txt   Uncomment this if you want to remove binList report.
      fi
    fi
  done
done
