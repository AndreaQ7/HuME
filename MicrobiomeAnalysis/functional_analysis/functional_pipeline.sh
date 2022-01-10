#!/usr/bin/bash

# Functional Analysis will be made using Humann2 toolkit. ( https://huttenhower.sph.harvard.edu/humann2/ ) and using the 
# library UniRef90EC

mkdir UniRef90EC
# Use humann and download humann2 databases
humann2_databases --download uniref uniref90_ec_filtered_diamond UniRef90EC/
humann2_database --download chocophlan full ChocophlanDB/

# Run humann2 on previous analyzed samples (56) of dental calculus using a while loop on all fastq trimmed files
# and providing the taxonomic profiles and then inferring the functional traits of the microbiota community per each sample
bash doHumann2.sh


# We get 3 output per sample: ${sample}_pathcoverage.tsv, ${sample}_pathabundance.tsv, ${sample}_genefamilies.tsv
# Now want to normalize the counts of each file for each sample: so we run this script based on humann2_renorm_table using cpm (Counts
# per million)
bash normalizeHumann2.sh

# Now join normalized tables
humann2_join_tables -i NormalizedGeneFamilies/ --file_name genefamilies_norm.tsv -o Joined_Gene_families.txt
humann2_join_tables -i NormalizedPathAbundances/ --file_name pathabundance_norm.tsv -o Joined_Path_abundance.txt
humann2_join_tables -i NormalizedPathCoverage/ --file_name pathcoverage_norm.tsv -o Joined_Path_coverage.txt


# Download utility mappings from humann2 databases to get a KEGG id output
humann2_databases --download utility_mapping full UniRef90EC/

# Rename Uniref90 to KEGG Orthology IDs
humann2_regroup_table -i Joined_Gene_families.txt -g uniref90_rxn -o ${output_name} -c ${path/to/utility_mapping/map_ko_uniref90.txt.gz}
