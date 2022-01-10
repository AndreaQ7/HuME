#!/usr/bin/bash


# CheckM is then used to asses the quality of obtained bins, resulting in good quality assemblies through the evaluation
# of completness and contamination


# The lineage workflow runs checkm on
checkm lineage_wf \
	-f checkm_output.txt \
	-t 30 \
	-x fa \
	"${input_folder}" checkm_output
	
checkm qa \
	-o 2 \
	-t 30 \
	-f checkm_stats \
	--tab_table checkm_output/lineage.ms \
	checkm_output
	
#calculate coverage per contig
checkm coverage \
	"${sample_name}_final.contigs.fa.metabat-bins0" \
	-x fa \
	-t 35 \
	"${sample_name}_coverage.tsv" \
	"${sample_name}.sorted.bam"


