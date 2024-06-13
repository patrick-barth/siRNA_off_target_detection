#!/bin/bash

# This script downloads all genomes and if available, all annotations and CDS which are deposited in https://docs.google.com/spreadsheets/d/1RIiyyyTIbz6HRRbVrbnf56rNVvgsLUk3FSJ4xjt5c2k/edit#gid=0
#  These data are then processed into a database that can be used to find potentiial siRNA Off-Target sites

# Variables
DATE=$(date '+%Y-%m-%d_%H-%M-%S')
WORK_DIR=$(pwd)
TMP_DIR="${WORK_DIR}/tmp"
OUTPUT_DIR="${WORK_DIR}/insect_genomes_database_$DATE/"
#SPREADSHEET_LINK="https://docs.google.com/spreadsheets/d/1RIiyyyTIbz6HRRbVrbnf56rNVvgsLUk3FSJ4xjt5c2k/export?format=csv"
SPREADSHEET_LINK="https://docs.google.com/spreadsheets/d/1YwQjEsfl1H3kDUFMwPjYQLgqeUrrrY0MBp0Q5n4NKhs/export?format=csv"
ID_COLUMN="assembly_ID"
GROUP_COLUMN="group"
PATH_TO_DATASETS="/homes/pbarth/software/RefSeq_downloader"

# TODO: Add hash to tmp directory
# TODO: Add tests if directories already exist
# TODO: Add argument to potentially not delete the tmp-directory
# TODO: put unwanted tool outputs to /dev/null



download_data() {
	mkdir ${TMP_DIR}
	# Download spreadsheet as a CSV file
	curl -L ${SPREADSHEET_LINK} > ${TMP_DIR}/insect_genomes_overview.csv
	# Get all groups from the corresponding column
	DB_GROUPS=$(awk -v group_col="$GROUP_COLUMN" -F ',' '
		NR>1 {
			printf $col" "
		} 
		NR==1 {
			for(i=1; i<=NF; i++) {
				if($i == group_col) col=i
			}
		}
		' "$TMP_DIR"/insect_genomes_overview.csv)
	# Make group names unique (only present once in array)
	declare -a DB_GROUPS_UNIQUE=()
	for element in ${DB_GROUPS[@]}; do
		if [[ ! ${DB_GROUPS_UNIQUE[@]} =~ ${element} ]]; then
			DB_GROUPS_UNIQUE+=($element)
		fi
	done
	echo "Groups found: ${DB_GROUPS_UNIQUE[@]}"
	# Go through every group, download the genomes
	for current_group in ${DB_GROUPS_UNIQUE[@]}; do
		echo "Processing insect group ${current_group}"
		# Extract all assembly IDs from the corresponding column which also belong to the group extracted previously
		GENOME_IDS=$(awk -v id_col="${ID_COLUMN}" -v group_col="${GROUP_COLUMN}" -v current_group="${current_group}" -F ',' '
			NR>1 {
				if($group_index == current_group) printf $id_index" "
			} 
			NR==1 {
				for(i=1; i<=NF; i++) {
					if($i == id_col) id_index=i
					if($i == group_col) group_index=i
				}
			}' ${TMP_DIR}/insect_genomes_overview.csv)
		
		# Generate directory within tmp to download genomes
		DIR_GROUP_TMP=${TMP_DIR}/${current_group}
		mkdir ${DIR_GROUP_TMP}
		# Download genomes 
		${PATH_TO_DATASETS}/datasets download genome accession $GENOME_IDS --include gff3,cds,genome,seq-report --dehydrated --filename ${DIR_GROUP_TMP}/insect_genomes.zip
		# Unzip downloaded genomes
		unzip ${DIR_GROUP_TMP}/insect_genomes.zip -d ${DIR_GROUP_TMP}
		# Rehydrate downloaded genomes (necessary for bulk download to restore all files)
		${PATH_TO_DATASETS}/datasets rehydrate --directory ${DIR_GROUP_TMP}
		
		# Generate file to collect all reference genomes
		PATH_COLLECTED_FASTA="${DIR_GROUP_TMP}/genomes_collected_${current_group}.fna"
		touch ${PATH_COLLECTED_FASTA}

		# Go through every genome
		for GENOME_ID_TMP in ${GENOME_IDS[@]}; do
			# Get path of current accession
			DIR_CURRENT_ID="${DIR_GROUP_TMP}/ncbi_dataset/data/${GENOME_ID_TMP}"
			# Get path of genome file
			GENOME_FILE=$(ls ${DIR_CURRENT_ID}/${GENOME_ID_TMP}*.fna 2>/dev/null )
			if [ -f "$GENOME_FILE" ]; then
				echo "Modifying genome ${GENOME_ID_TMP}"
				# Modify genome file to add the accession to the first string of the header. This way they can be easily assigned in later steps
				awk -v accession=${GENOME_ID_TMP} '{
					if (/^>/) {
						sub(/ /, "_" accession " ", $0)
					}; print $0
				}' ${GENOME_FILE} > temp_file.fa && mv temp_file.fa ${GENOME_FILE}
				# Copy genome to collect fasta file
				cat ${GENOME_FILE} >> ${PATH_COLLECTED_FASTA}
				# Check if annotation is present
				if [ -f "${DIR_CURRENT_ID}/genomic.gff" ]; then
					echo "Annotation found for ${GENOME_ID_TMP}"
					mv ${DIR_CURRENT_ID}/genomic.gff ${DIR_CURRENT_ID}/${GENOME_ID_TMP}.gff
					# TODO: Dir needs to be adapted. Maybe put whole file generation to DB_build
					echo -e "${GENOME_ID_TMP}\tTRUE\t$(realpath ${DIR_CURRENT_ID}/${GENOME_ID_TMP}.gff)" >> ${PATH_ANNOTATION_INFO}
				else
					echo "No annotation found for ${GENOME_ID_TMP}"
					echo -e "${GENOME_ID_TMP}\tFALSE\tNA" >> ${PATH_ANNOTATION_INFO}
				fi
			else
				echo "Genome file belonging to accession ${GENOME_ID_TMP} seems to missing"
			fi
		done
		echo -e "\n"
	done
}
	
	
	
	# Get important files
	##mkdir "$OUTPUT_DIR"/raw_data
	##mkdir "$OUTPUT_DIR"/metadata
	##cp "$TMP_DIR"/insect_genomes.zip "$OUTPUT_DIR"/raw_data
	##cp "$TMP_DIR"/insect_genomes_overview.csv "$OUTPUT_DIR"/metadata
	##cp "$TMP_DIR"/ncbi_dataset/fetch.txt "$OUTPUT_DIR"/metadata
	# Get md5sums
	##touch "$OUTPUT_DIR"/metadata/md5sum.tsv
	##md5sum "$OUTPUT_DIR"/metadata/insect_genomes_overview.csv >> "$OUTPUT_DIR"/metadata/md5sum.tsv
	##md5sum "$OUTPUT_DIR"/raw_data/insect_genomes.zip >> "$OUTPUT_DIR"/metadata/md5sum.tsv
	##md5sum "$TMP_DIR"/ncbi_dataset/data/GC{A,F}_*/*.{*_genomic.fna,gff} >> "$OUTPUT_DIR"/metadata/md5sum.tsv


generate_indexes() {
	mkdir ${OUTPUT_DIR}
	# Get list of all groups in tmp dir as they were downloaded during download data 
		# (list entries | filter for dirs | reverse | get first column (dir name) | reverse again | turn newline to space in order to get an array)
	GROUPS_LIST=$(ls -l ${TMP_DIR} | grep "^d" | rev | cut -d ' ' -f 1 | rev | tr '\n' ' ')
	# Go through all accessions of the current group
	for CURRENT_GROUP in ${GROUPS_LIST[@]}; do
		TMP_DIR_CURRENT_GROUP="${TMP_DIR}/${CURRENT_GROUP}"
		OUTPUT_DIR_CURRENT_GROUP="${OUTPUT_DIR}/${CURRENT_GROUP}"
		GENOMES_LIST=$(ls -l ${TMP_DIR_CURRENT_GROUP} | grep "^d" | rev | cut -d ' ' -f 1 | rev | tr '\n' ' ')
		
		ANNOTATION_DIR="${OUTPUT_DIR_CURRENT_GROUP}/annotations"
		CDS_DIR="${OUTPUT_DIR_CURRENT_GROUP}/cds"
		DB_DIR="${OUTPUT_DIR_CURRENT_GROUP}/db"
		mkdir ${OUTPUT_DIR_CURRENT_GROUP}
		mkdir ${ANNOTATION_DIR}
		mkdir ${CDS_DIR}
		mkdir ${DB_DIR}

		# Generate and prepare file to keep information about the presence of annotations
		PATH_ANNOTATION_INFO="${OUTPUT_DIR_CURRENT_GROUP}/info_annotation_${CURRENT_GROUP}.tsv"
		touch ${PATH_ANNOTATION_INFO}
		echo -e "accession\tannotation\tpath" > ${PATH_ANNOTATION_INFO}

		# Go through every accession
		for CURRENT_ACCESSION in ${GENOMES_LIST[@]}; do
			# TODO: Continue here:
				# Copy part from download function regarding annotations
				# Copy CDS and rename
				# Copy collected fasta
				# Get md5 checksum
				# Build DB
		done
		
	done


	# Build bowtie1 index (used for short reads)
	mkdir "$OUTPUT_DIR"/index_bowtie
	bowtie-build --threads 1 $GENOMES_LIST "$OUTPUT_DIR"/index_bowtie/insect_genomes
}



# Start functions
##download_data
generate_indexes
