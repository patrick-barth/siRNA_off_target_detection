#!/bin/bash

# This script downloads all genomes and if available, all annotations and CDS which are deposited in https://docs.google.com/spreadsheets/d/1RIiyyyTIbz6HRRbVrbnf56rNVvgsLUk3FSJ4xjt5c2k/edit#gid=0
#  These data are then processed into a database that can be used to find potentiial siRNA Off-Target sites

# Variables
DATE=$(date '+%Y-%m-%d_%H-%M-%S')
TMP_DIR="./tmp/"
OUTPUT_DIR="./insect_genomes_database_$DATE)/"
SPREADSHEET_LINK="https://docs.google.com/spreadsheets/d/1RIiyyyTIbz6HRRbVrbnf56rNVvgsLUk3FSJ4xjt5c2k/export?format=csv"
ID_COLUMN="assembly_ID"




main() {
	mkdir $TMP_DIR
	mkdir $OUTPUT_DIR
	# Download spreadsheet as a CSV file
	curl -L "$SPREADSHEET_LINK" > "$TMP_DIR"/insect_genomes_overview.csv
	# Extract all assembly IDs from the corresponding column
	GENOME_IDS=$(awk -v id_col="$ID_COLUMN" -F ',' 'NR>1 {printf $col" "} NR==1 {for(i=1; i<=NF; i++) if($i == id_col) col=i}' "$TMP_DIR"/insect_genomes_overview.csv)
	# Download all genomes
	~/software/RefSeq_downloader/datasets download genome accession $GENOME_IDS --include gff3,cds,genome,seq-report --dehydrated --filename "$TMP_DIR"/insect_genomes.zip
	# unzip downloaded data 
	unzip "$TMP_DIR"/insect_genomes.zip -d "$TMP_DIR"
	# Rehydrate datasets (needed due to download of large datasets)
	~/software/RefSeq_downloader/datasets rehydrate --directory "$TMP_DIR"
	# Get important files
	mkdir "$OUTPUT_DIR"/raw_data
	mkdir "$OUTPUT_DIR"/metadata
	cp "$TMP_DIR"/insect_genomes.zip "$OUTPUT_DIR"/raw_data
	cp "$TMP_DIR"/insect_genomes_overview.csv "$OUTPUT_DIR"/metadata
	cp "$TMP_DIR"/ncbi_dataset/fetch.txt "$OUTPUT_DIR"/metadata
	# Get md5sums
	touch "$OUTPUT_DIR"/metadata/md5sum.tsv
	md5sum "$OUTPUT_DIR"/metadata/insect_genomes_overview.csv > "$OUTPUT_DIR"/metadata/md5sum.tsv
	md5sum "$OUTPUT_DIR"/metadata/insect_genomes.zip > "$OUTPUT_DIR"/metadata/md5sum.tsv
	md5sum "$TMP_DIR"/ncbi_dataset/data/GC{A,F}_*/*.{1_genomic.fna,gff} > "$OUTPUT_DIR"/metadata/md5sum.tsv
}

# Start main function
main
