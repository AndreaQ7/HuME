#!/usr/bin/bash


# This script runs bowtie2 and samtools to remap the trimmed and quality filtered reads on the new generated contigs.fasta files
# by spades.py in --meta mode.

# Create results folder for bam files
mkdir results/bamresults/

# Create list variable
list_files="/path/to/trimmed_files.txt"

# Iterate the mapping workflow on each sample
while read namesample; do
  # create a folder for each sample
  echo $namesample
  cat results/metaspadesOut/$namesample/contigs.fasta | head -2
  mkdir results/bamresults/$namesample

  # Run bowtie build on contigs to index them and then map the reads (forward and reverse) against contigs
  bowtie2-build results/metaspadesOut/$namesample/contigs.fasta results/bamresults/$namesample/$namesample
  bowtie2 --sensitive-local -p 10 -x results/bamresults/$namesample/$namesample -1 results/trimmed/$namesample/$namesample.pair1.truncated.fastq -2 results/trimmed/$namesample/$namesample.pair2.truncated.fastq -S results/bamresults/$namesample/$namesample.sam
  samtools faidx results/metaspadesOut/$namesample/contigs.fasta
  # samtools import results/metaspadesOut/$namesample/contigs.fasta.fai results/bamresults/$namesample/$namesample.sam results/bamresults/$namesample/$namesample.bam
  # Generate the bam file from sam, and then delete the bigger sam file.
  samtools view -S -b results/bamresults/$namesample/$namesample.sam > results/bamresults/$namesample/$namesample.bam -@ 10
  rm results/bamresults/$namesample/$namesample.sam

  # Sort and index the mapping file
  samtools sort -@ 10 results/bamresults/$namesample/$namesample.bam -o results/bamresults/$namesample/$namesample.sorted.bam
  samtools index -@ 10 results/bamresults/$namesample/$namesample.sorted.bam
  samtools idxstats results/bamresults/$namesample/$namesample.sorted.bam > results/bamresults/$namesample/idxstats.txt
done < $list_files

