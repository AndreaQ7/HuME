# Example using sample bin 1 fasta from Tr-USD177_S2_L001
# Align reads to contigs
bowtie2-build Tr-USD177_S2_L001.bin.1.fa Tr-USD177_S2_L001.bin.1
bowtie2 -t -x  Tr-USD177_S2_L001.bin.1 -U  /mnt/hume/TrimmedDentalCalculus/Fer-T75_S25_001/Fer-T75_S25_001.collapsed.gz -p 30 -S Tr-USD177_S2_L001.bin.1.sam
	
# Output only aligned reads in bam format
samtools view -bSF4 Tr-USD177_S2_L001.bin.1.sam > Tr-USD177_S2_L001.bin.1.bam -@ 30

# Sort bam file
samtools sort Tr-USD177_S2_L001.bin.1.bam -o Tr-USD177_S2_L001.bin.1.sorted.bam

# Index BAM file
samtools index Tr-USD177_S2_L001.bin.1.sorted.bam

# Count number of reads aligned to each contig
samtools idxstats Tr-USD177_S2_L001.bin.1.sorted.bam > Tr177_num_reads_per_contig.txt

# Based on third column, remove rows that have # of reads <1000 
awk -F"\t" '$3>999' Tr177_num_reads_per_contig.txt | cut -f1 > Tr177_contigs.txt


while read contigs;
do
	contig=$(echo "$contigs")
samtools view -H Tr-USD177_S2_L001.bin.1.sorted.bam | grep -P -v "^@SQ" > Tr-USD177_S2_L001.bin.1_${contig}.header
samtools view -H Tr-USD177_S2_L001.bin.1.sorted.bam | grep -P "^@SQ\tSN:${contigs}\t" >> Tr-USD177_S2_L001.bin.1_${contig}.header
samtools view -bh Tr-USD177_S2_L001.bin.1.sorted.bam ${contig} > Tr-USD177_S2_L001.bin.1_${contig}.bam
samtools view -h Tr-USD177_S2_L001.bin.1_${contig}.bam | python2 ../../../script/PMDtools/pmdtools.0.60.py  --threshold 1 --header | samtools view -Sb -> pmd1_${contigs}_filtered.bam
samtools view "pmd1_${contigs}_filtered.bam" |  python2 ../../../script/PMDtools/pmdtools.0.60.py  --platypus --requiremapq 30 > "PMD_${contigs}_temp.txt"
done <  Tr177_contigs.txt


