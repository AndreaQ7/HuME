#!/usr/bin/bash

# This script downloads all the reference genomes of the relevant microrganisms and put them in their homologous folder.

# Run in debugging mode
set -uex

# Create the variable with the list of microrganisms
list_microrganims="/path/to/list.txt"
suffix="_genomic.fna.gz"

# Enter the while loop and use the downloaded assembly summary to get the links to the genomes of each microrganism
# Download the assembly summary by:
wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

while read l; do
  l2=${l//_/ }
  oral_sp="sp. oral taxon"
  if [[ "$l2" == *"$oral_sp"* ]]; then
    awk -v l2="$l2" -F '\t' '$8~l2 {print $0; exit}' assembly_summary_genbank.txt >> strainInfo.txt
    ftp=$(awk -v l2="$l2" -F '\t' '$8~l2 {print $20; exit}' assembly_summary_genbank.txt)
    file=$(echo $ftp| cut -d'/' -f 10)
    fasta=$file$suffix
  else  
    awk -v l2="$l2" -F '\t' '$8~l2 && $5=="representative genome" {print $0; exit}' assembly_summary_genbank.txt >> strainInfo.txt
    ftp=$(awk -v l2="$l2" -F '\t' '$8~l2 && $5=="representative genome" {print $20; exit}' assembly_summary_genbank.txt)
    file=$(echo $ftp| cut -d'/' -f 10)
    fasta=$file$suffix
  fi
  if test -z $ftp
  then
    echo "No ftp found for this species"
    echo ${l} >> species_not_found.txt
  else
    mkdir $l
    cd $l
    wget $ftp/$fasta
    cd ../
  fi
done < $list_microrganisms

# Please note: some of the reference genomes are not indicated in the assembly_summary_genbank.txt as "representative genome" so 
# they need to be downloaded manually. For this purpose, a file "species_not_found.txt will be generated in the previous step,
# while a file with all the downloaded strain info will be available as strainInfo.txt
