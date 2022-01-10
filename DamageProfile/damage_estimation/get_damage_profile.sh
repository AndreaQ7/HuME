#!/usr/bin/bash	


# This is the script that will be used to compute the damage profiles on each sample for the most relevant species found with
# Kraken2.
# Tools used in this part are bowtie2, samtools, bedtools, 2 commands from the suite GATK and PMDtools for damage profiling.

# Activate the environment in which Mapdamage binary is present
# FOR USERS: comment out this line or change the name of the environment with yours
source activate bam-filter-2020.3

# Create bamresults folder
mkdir bamresults
bampath="/path/to/bamresults"
fastqPath="/path/to/fastq_trimmed"

# This for loop iterates on all folders created within previous script
for file in $(ls .);
do
  
  # Activate bowtie and samtools environment and unzip the reference files saved in each organism folder
  source activate bam-filter-2020.3
  if [ -d "$file" ]; then
    echo $file;
    gz="fna.gz"
    fasta=$(ls $file)
    if [[ "$fasta" == *"$gz"* ]]; then
      gunzip $file/$fasta;
      ref2=$fasta;
    else
      echo "no file zipped"
      ref2=$(ls $file);
    fi
    echo $ref2;
    
    # Reference sequence preparation: the first step is the fasta indexing using the Borrow-Wheeler transform. 4 index files will be produced
    # The -a flag stands for algoritm and the is algoritm is efficient for genomes shorter than 2Gb, such as bacterial genomes
    bwa index -a is $file/*.fna;
    
    # Here we create a dictionary file for the reference genome, that will be then used from GATK RealignerTargetCreator tool
    picard CreateSequenceDictionary R=$file/$ref2 O=$file/${ref2}.dict;
    
    # Also samtools Index will be run to make sequence regions accessible from other tools (GATK).
    samtools faidx $file/*.fna;
    
    # Create the folder of the given bacteria inside bamresults folder
    mkdir $bampath/$file;
    
    # Create a folder in each bacteria file to store the damage calculations and also a path variable
    mkdir ${file}/damageCalculation;
    damagepath="/path/to/${file}/damageCalculation";
    
    # The alignment algoritm is that provided from bwa aln that use a end-to-end alingment and it's optimal for short reads, such as 
    # ancient ones. Flag -n is the maximum edit distance and we set it to be very stringent (we do not want that similar reads of other
    # species also match our reference), the -o is the max number of gap opens (increased over the default value to 2, tolerating more gaps)
    # finally -l is the seed length, that is turned off by indicating a huge number; -e is set to -1 to disable long gaps.
    for sample in $(ls $fastqPath); do
      # Use bwa aln to create the sam index file
      bwa aln $file/*.fna $fastqPath$sample/${sample}.collapsed.gz -n 0.1 -l 1000 -t 5 > $bampath/$file/${sample}.$file.sai;
      
      # Do the alignment with bwa samse to obtain the sam file
      bwa samse $file/*.fna $bampath/$file/${sample}.$file.sai $fastqPath$sample/${sample}.collapsed.gz -f $bampath/$file/${sample}.$file.sam;

      # Convert the sam file in bam using samtools
      samtools view -Sb $bampath/$file/${sample}.$file.sam > $bampath/$file/${sample}.$file.bam;

      # Removing samfile for some space
      rm $bampath/$file/${sample}.$file.sam
      
      # Filtering step 1: selecting reads with alignment quality > 30
      samtools view -q 30 -b $bampath/$file/${sample}.$file.bam > $bampath/$file/${sample}.${file}_qfilter30.bam;

      # Filtering step 2 : selecting reads that are best hits with no more than 1 match (no alternative matching reads)
      # bamtools filter -tag X0:1 -in $bampath$file/${sample}.${file}_qfilter20.bam -out $bampath$file/${file}_X0.bam;
      # Sort and index the bam file in order to make it accessible to other softwares
      samtools sort $bampath/$file/${sample}.${file}_qfilter30.bam -o $bampath/$file/${sample}.${file}.sort.bam;
      echo "end sort ${file}";
      samtools index $bampath/$file/${sample}.${file}.sort.bam;

      # Create the Read groups indexes, to distinguish samples and remove duplicates using Picard tools. Then index with samtools
      picard AddOrReplaceReadGroups INPUT=$bampath/$file/${sample}.${file}.sort.bam OUTPUT=$bampath/$file/${sample}.${file}.RG.bam RGID=rg_id RGLB=lib_id RGPL=platform RGPU=plat_unit RGSM=sam_id VALIDATION_STRINGENCY=LENIENT;
      samtools index $bampath/$file/${sample}.${file}.RG.bam;
      echo "end index1";
      picard MarkDuplicates I=$bampath/$file/${sample}.${file}.RG.bam O=$bampath/$file/${sample}.${file}.DR.bam M=output_metrics.txt REMOVE_DUPLICATES=True VALIDATION_STRINGENCY=LENIENT &>  $bampath/$file/${sample}.${file}_picardlogFile.log;
      samtools index $bampath/$file/${sample}.${file}.DR.bam;
      echo "end index 2";
      source activate alignment-2020.3;
      
      # We used from GATK suite the RealignertargeCreator and IndelRealigner tools:
      # This is used to realign some parts of the reads whose mismatches are actually no SNPs but accidental misalignment
      # For eny documentation: https://gatk.broadinstitute.org/hc/en-us
      java -jar $GATK -T RealignerTargetCreator -R $file/$ref2 -I $bampath/$file/${sample}.${file}.DR.bam -o $bampath/$file/${sample}.${file}_target.intervals --filter_mismatching_base_and_quals;
      java -jar $GATK -T IndelRealigner -R $file/$ref2 -I $bampath/$file/${sample}.${file}.DR.bam -targetIntervals $bampath/$file/${sample}.${file}_target.intervals -o $bampath/$file/${sample}.${file}.final.bam --filter_bases_not_stored;
      echo "end java ${sample}.${file}";
      # Then sort and index the new realigned files
      samtools sort $bampath/$file/${sample}.${file}.final.bam -o $bampath/$file/${sample}.${file}.final.sort.bam;
      samtools index $bampath/$file/${sample}.${file}.final.sort.bam;
      samtools flagstat $bampath/$file/${sample}.${file}.final.sort.bam > $bampath/$file/${sample}.${file}_flagstat.txt;
      echo "end ${sample}.${file}";

      # Activate new environment to use PMDtools
      # PMDtools will be used to infer the frequency of G to A and C to T misincorporations.
      source activate adna-2020.3;
      samtools view -h $bampath/$file/${sample}.${file}.final.sort.bam | python2 /path/to/Scripts/PMDtools/pmdtools.0.60.py --threshold 1 --header --requiremapq 30 | samtools view -Sb > /home/ginnocenti/DentalCalculus/MAPDAMAGE/PMDtools/$file/damageCalculation/pmd1_${sample}.${file}_filteredX0.bam;
      bedtools bamtobed -i /home/ginnocenti/DentalCalculus/MAPDAMAGE/PMDtools/$file/damageCalculation/pmd1_${sample}.${file}_filteredX0.bam -tag NM | awk -F '\t' '{print $4,$5}' > /home/ginnocenti/DentalCalculus/MAPDAMAGE/PMDtools/$file/damageCalculation/pmd1_${sample}.${file}.NM.bed # Edit distance in quinta colonna quinta colonna > ricavare sulla edit distance una frequency table.
      samtools view ${damagepath}/*_filteredX0.bam | python2 /path/to/Scripts/PMDtools/pmdtools.0.60.py --platypus --requirebaseq 30 > ${damagepath}/PMD_${sample}.${file}_temp.txt;  
      echo "selezione fatta per ${sample}.${file}";
      # Plot with  r script from PMDtools
      R CMD BATCH path/to/Scripts/PMDtools/plotPMD.v2.R;
      mv PMD_plot.frag.pdf ${damagepath}/PMD.plot.${sample}.${file}.pdf;
      
      # Now sort the final bam
      samtools sort ${damagepath}/pmd1_${sample}.${file}_filteredX0.bam -o ${damagepath}/${file}_filtered.final.sortX0.bam;
      
      # CHECK THESE STEPS to know if is necessary to remove them
      mv PMD_temp.txt $damagepath
      mv plotPMD.v2.Rout $damagepath
      mv output_metrics.txt
      source activate bam-filter-2020.3;
    done
  else 
    echo "txt file";
  fi
done
