# This is the pipeline that will be used over Dental Calculus Samples using Kraken, and that is composed by 2 main steps:
# - Building of a new custom Kraken DB
# - Mapping of Dental Calculus Adapter trimmed reads against out Custom DB
# The Database will be composed by all Bacteria and Archea Genomes and Human mitocondrial sequences that could be useful for further
# investigation.

# FOR USERS: this is not a ready-to-use wrapper but needs some changes the "/path/to/*" with your actually available paths in your machine.
# Plus it needs some manually installing as specified in commented lines.

# Building of Kraken2 Database in custom mode
mkdir KrakenDB
# Activate kraken environment
source activate metaprocessing{}
kraken2-build --download-taxonomy --threads 10 --db KrakenDB
kraken2-build --download-library bacteria --threads 20 --db KrakenDB
kraken2-build --download-library archaea --threads 20 --db KrakenDB

# Now add the human mithocondrion sequence retrieved from ncbi db at https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1?report=fasta 
# with the accession id NC_012920.1 (Homo sapians itochondrion, complete sequence). Download it and put it in /mnt/hume/Gabriel/mtDNA_human
kraken2-build --add-to-library /path/to/mitochondria_human.fasta --threads 4 --db KrakenDB

# Build the final Database
kraken2-build --build --max-db-size 10000000000 --db KrakenDB

# Use the doKraken.sh script
cd DentalCalculus
mkdir krakenOut
mkdir krakenReport
bash doKraken.sh

# Download bracken
git clone https://github.com/jenniferlu717/Bracken.git
# install
cd Bracken
bash install_bracken.sh

# Insert bracken scripts in your path
PATH=$PATH:/path/to/Bracken/

# Use bracken-build to assign each kmer in Kraken database to a taxonomic unit, which will be used for computing taxonomic 
# abundance on every level
# Using kmers of length 35 and reads of 50 bp
bracken-build -d KrakenDB -t 10 -k 35 -l 50

# These file will be created inside KrakenDB
# --- database50mers.kmer_distrib 
# --- database50mers.kraken
# --- hash.k2d
# --- opts.k2d
# --- seqid2taxid.map
# --- taxo.k2d

# Use Bracken: I used reads of 50bp and kmers of 35 bp (Kraken2 default)
mkdir brackenOut
bash doBracken.sh
# FOR USERS: change script in order to have Species, Genus or Phylum level abundance calculation: default will create all these and
# store output in brackenOut/PhylaAbundance, brackenOut/SpeciesAbundance, brackenOut/GenusAbundance

# FOR USERS: Then git clone Kraken Tools in home and add both bracken bin directory and KrakenTools bin directory to my env in .bashrc file

# Finally, convert kraken reports to 
bash doKreport2mpa.sh
# And merge it with 
combine_mpa.py -i kraken2mpa/*.txt -o CombinedKraken2mpa.txt


##### Repeat the pipeline also for SRR samples from HMP #####
# Created 2 new folders for kraken outputs: ~/DentalCalculus/KrakenSRRreport and ~/DentalCalculus/krakenSRROut
