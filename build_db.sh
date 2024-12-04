#!/bin/bash

# This script downloads all genomes and if available, all annotations and CDS which are deposited in https://docs.google.com/spreadsheets/d/1RIiyyyTIbz6HRRbVrbnf56rNVvgsLUk3FSJ4xjt5c2k/edit#gid=0
#  These data are then processed into a database that can be used to find potentiial siRNA Off-Target sites

# Variables
DATE=$(date '+%Y-%m-%d_%H-%M-%S')
WORK_DIR=$(pwd)
TMP_DIR="${WORK_DIR}/tmp"
TMP_DIR_GROUPS="${TMP_DIR}/groups"
LOG_FILE="${TMP_DIR}/logs.txt"
OUTPUT_DIR="${WORK_DIR}/insect_genomes_database_$DATE/"
##SPREADSHEET_LINK="https://docs.google.com/spreadsheets/d/1RIiyyyTIbz6HRRbVrbnf56rNVvgsLUk3FSJ4xjt5c2k/export?format=csv"
##SPREADSHEET_LINK="https://docs.google.com/spreadsheets/d/1YwQjEsfl1H3kDUFMwPjYQLgqeUrrrY0MBp0Q5n4NKhs/export?format=csv"
SPREADSHEET_LINK_NCBI="https://docs.google.com/spreadsheets/d/12LIbdqf7D1Lq7pppJMQw4u2hk9N2CxS5CapoZsJn_us/export?format=csv"
SPREADSHEET_LINK_INSECT_BASE="https://docs.google.com/spreadsheets/d/12LIbdqf7D1Lq7pppJMQw4u2hk9N2CxS5CapoZsJn_us/export?format=csv&gid=601092575"
ID_COLUMN="assembly_ID"
SPECIES_COLUMN="Species name"
GROUP_COLUMN="group"
PATH_TO_DATASETS="/homes/pbarth/software/RefSeq_downloader"
INSECTBASE_LINK_START="http://v2.insect-genome.com/api/Download/..-01_data-01_species-"
INSECTBASE_LINK_END_GENOME=".genome.fa.tar.bz2"
INSECTBASE_LINK_END_ANNOTATION=".gff3"
INSECTBASE_LINK_END_CDS=".cds.fa"

# TODO: Add hash to tmp directory
# TODO: Add tests if directories already exist
# TODO: Add argument to potentially not delete the tmp-directory
# TODO: put unwanted tool outputs to /dev/null

PATH_TO_INFORMATION="${TMP_DIR}/infos.tsv"

mkdir ${TMP_DIR}
mkdir ${TMP_DIR_GROUPS}
touch ${PATH_TO_INFORMATION}
touch ${LOG_FILE}
echo -e "accession\tgroup\torigin\tannotation\tannotation_path\tCDS\tCDS_path" > ${PATH_TO_INFORMATION}

download_data_ncbi() {

	echo "Downloading data from NCBI" | tee -a ${LOG_FILE}
	
	DOWNLOAD_DIR_TMP="${TMP_DIR}/NCBI"
	mkdir ${DOWNLOAD_DIR_TMP}
	# Download spreadsheet as a CSV file
	curl -L ${SPREADSHEET_LINK_NCBI} > ${TMP_DIR}/insect_genomes_overview_ncbi.csv
	# Get all groups from the corresponding column
	DB_GROUPS=$(awk -v group_col="${GROUP_COLUMN}" -F ',' '
		NR>1 {
			printf $col" "
		} 
		NR==1 {
			for(i=1; i<=NF; i++) {
				if($i == group_col) col=i
			}
		}
		' "${TMP_DIR}"/insect_genomes_overview_ncbi.csv)
	# Make group names unique (only present once in array)
	declare -a DB_GROUPS_UNIQUE=()
	for ELEMENT in ${DB_GROUPS[@]}; do
		if [[ ! ${DB_GROUPS_UNIQUE[@]} =~ ${ELEMENT} ]]; then
			DB_GROUPS_UNIQUE+=(${ELEMENT})
		fi
	done
	echo "Groups found: ${DB_GROUPS_UNIQUE[@]}" | tee -a ${LOG_FILE}
	# Go through every group, download the genomes
	for current_group in ${DB_GROUPS_UNIQUE[@]}; do
		echo "Processing insect group ${current_group}" | tee -a ${LOG_FILE}
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
			}' ${TMP_DIR}/insect_genomes_overview_ncbi.csv)
		
		# Generate directory within tmp to download genomes
		DIR_GROUP_TMP=${TMP_DIR_GROUPS}/${current_group}
		if [ ! -d "${DIR_GROUP_TMP}" ]; then
			mkdir ${DIR_GROUP_TMP}
		fi

		# Download genomes 
		${PATH_TO_DATASETS}/datasets download genome accession ${GENOME_IDS} --include gff3,cds,genome,seq-report --dehydrated --filename ${DOWNLOAD_DIR_TMP}/insect_genomes.zip
		# Unzip downloaded genomes
		unzip ${DOWNLOAD_DIR_TMP}/insect_genomes.zip -d ${DOWNLOAD_DIR_TMP}
		# Rehydrate downloaded genomes (necessary for bulk download to restore all files)
		${PATH_TO_DATASETS}/datasets rehydrate --directory ${DOWNLOAD_DIR_TMP}
		
		# Generate file to collect all reference genomes
		PATH_COLLECTED_FASTA="${DIR_GROUP_TMP}/genomes_collected_${current_group}.fna"
		if [ ! -f "${PATH_COLLECTED_FASTA}" ]; then
			touch ${PATH_COLLECTED_FASTA}
		fi

		# Go through every genome
		for GENOME_ID_TMP in ${GENOME_IDS[@]}; do
			# Set variables for info file
			ANNOTATION_FOUND='FALSE'
			OUTPUT_ANNOTATION_RELATIVE='NA'
			CDS_FOUND='FALSE'
			OUTPUT_CDS_RELATIVE='NA'

			# Generate variable for output dir files with changed names
			TARGET_DIR="${DIR_GROUP_TMP}/${GENOME_ID_TMP}"
			TARGET_GENOME="${DIR_GROUP_TMP}/${GENOME_ID_TMP}/${GENOME_ID_TMP}.genome.fa"
			TARGET_ANNOTATION="${DIR_GROUP_TMP}/${GENOME_ID_TMP}/${GENOME_ID_TMP}.gff3"
			TARGET_CDS="${DIR_GROUP_TMP}/${GENOME_ID_TMP}/${GENOME_ID_TMP}.cds.fa"
			mkdir ${TARGET_DIR}
			# Get path of current accession
			DIR_CURRENT_ID="${DOWNLOAD_DIR_TMP}/ncbi_dataset/data/${GENOME_ID_TMP}"
			# Get path of genome file
			GENOME_FILE=$(ls ${DIR_CURRENT_ID}/${GENOME_ID_TMP}*.fna 2>/dev/null )
			
			if [ -f "${GENOME_FILE}" ]; then
				echo "Modifying genome ${GENOME_ID_TMP}" | tee -a ${LOG_FILE}
				# Modify genome file to add the accession to the first string of the header. This way they can be easily assigned in later steps
				awk -v accession=${GENOME_ID_TMP} '{
					if (/^>/) {
						sub(/ /, ":" accession " ", $0)
					}; print $0
				}' ${GENOME_FILE} > temp_file.fa && mv temp_file.fa ${TARGET_GENOME}
				# Copy genome to collect fasta file
				cat ${TARGET_GENOME} >> ${PATH_COLLECTED_FASTA}
			else
				echo "Genome file belonging to accession ${GENOME_ID_TMP} seems to missing" | tee -a ${LOG_FILE}
			fi

			# Check if an annotation file is present. If yes rename it
			CURRENT_ANNOTATION_FILE="${DIR_CURRENT_ID}/genomic.gff"
			if [ -f "${CURRENT_ANNOTATION_FILE}" ]; then
				# Check if there was an error with downloading the annotation file (first line is <!DOCTYPE html>)
				if grep -q "^<!DOCTYPE html>" "${CURRENT_ANNOTATION_FILE}"; then
					echo "Error in annotation file of ${GENOME_ID_TMP}: File begins with: <!DOCTYPE html>" | tee -a ${LOG_FILE}
					echo "File will be removed" | tee -a ${LOG_FILE}
					rm ${CURRENT_ANNOTATION_FILE}
				else
					echo "Annotation found for ${GENOME_ID_TMP}" | tee -a ${LOG_FILE}
					cp ${CURRENT_ANNOTATION_FILE} ${TARGET_ANNOTATION}
					ANNOTATION_FOUND='TRUE'
					OUTPUT_ANNOTATION_RELATIVE="./${current_group}/annotations/${GENOME_ID_TMP}.gff3"
				fi
			else
				echo "No annotation found for ${GENOME_ID_TMP}" | tee -a ${LOG_FILE}
			fi

			# Check if CDS file is present. If yes, rename it
			CURRENT_CDS_FILE="${DIR_CURRENT_ID}/cds_from_genomic.fna"
			if [ -f "${CURRENT_CDS_FILE}" ]; then
				# Check if there was an error with downloading the CDS file (first line is <!DOCTYPE html>)
				if grep -q "^<!DOCTYPE html>" "${CURRENT_CDS_FILE}"; then
					echo "Error in CDS file of ${GENOME_ID_TMP}: File begins with: <!DOCTYPE html>" | tee -a ${LOG_FILE}
					echo "File will be removed" | tee -a ${LOG_FILE}
					rm ${CURRENT_CDS_FILE}
				else
					echo "CDS found for ${GENOME_ID_TMP}" | tee -a ${LOG_FILE}
					cp ${CURRENT_CDS_FILE} ${TARGET_CDS}
					CDS_FOUND='TRUE'
					OUTPUT_CDS_RELATIVE="./${current_group}/cds/${GENOME_ID_TMP}.cds.fna"
				fi
			else
				echo "No CDS found for ${GENOME_ID_TMP}" | tee -a ${LOG_FILE}
			fi

			echo -e "${GENOME_ID_TMP}\t${current_group}\tNCBI\t${ANNOTATION_FOUND}\t${OUTPUT_ANNOTATION_RELATIVE}\t${CDS_FOUND}\t${OUTPUT_CDS_RELATIVE}" >> ${PATH_TO_INFORMATION}
		done
		echo -e "\n" | tee -a ${LOG_FILE}
	done
}
	
download_data_insectbase(){
	echo "Downloading data from insectBase" | tee -a ${LOG_FILE}
	# Download spreadsheet as a CSV file
	curl -L ${SPREADSHEET_LINK_INSECT_BASE} > ${TMP_DIR}/insect_genomes_overview_insectbase.csv
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
		' "${TMP_DIR}"/insect_genomes_overview_insectbase.csv)
	# Make group names unique (only present once in array)
	declare -a DB_GROUPS_UNIQUE=()
	for element in ${DB_GROUPS[@]}; do
		if [[ ! ${DB_GROUPS_UNIQUE[@]} =~ ${element} ]]; then
			DB_GROUPS_UNIQUE+=($element)
		fi
	done
	echo "Groups found: ${DB_GROUPS_UNIQUE[@]}" | tee -a ${LOG_FILE}
	# Go through every group, download the genomes
	for current_group in ${DB_GROUPS_UNIQUE[@]}; do
		echo "Processing insect group ${current_group}" | tee -a ${LOG_FILE}
		# Extract all species names from the corresponding column which also belong to the group extracted previously and substitute spaces with underscores (_)
		GENOME_IDS=$(awk -v species_col="${SPECIES_COLUMN}" -v group_col="${GROUP_COLUMN}" -v current_group="${current_group}" -F ',' '
			NR>1 {
				if($group_index == current_group) {
					gsub(/ /, "_", $species_index)
					printf $species_index" "
				}
			} 
			NR==1 {
				for(i=1; i<=NF; i++) {
					if($i == species_col) species_index=i
					if($i == group_col) group_index=i
				}
			}' ${TMP_DIR}/insect_genomes_overview_insectbase.csv)

		# Generate directory within tmp to download genomes
		DIR_GROUP_TMP="${TMP_DIR_GROUPS}/${current_group}"
		if [ ! -d "${DIR_GROUP_TMP}" ]; then
			mkdir ${DIR_GROUP_TMP}
		fi
		PATH_COLLECTED_FASTA="${DIR_GROUP_TMP}/genomes_collected_${current_group}.fna"

		for SPECIES in ${GENOME_IDS[@]}; do
			# Set variables for info file
			ANNOTATION_FOUND='FALSE'
			OUTPUT_ANNOTATION_RELATIVE='NA'
			CDS_FOUND='FALSE'
			OUTPUT_CDS_RELATIVE='NA'

			DIR_CURRENT_SPECIES="${DIR_GROUP_TMP}/${SPECIES}"

			FILE_GENOME_DOWNLOAD="${DIR_CURRENT_SPECIES}/${SPECIES}${INSECTBASE_LINK_END_GENOME}"
			GENOME_FILE="${DIR_CURRENT_SPECIES}/${SPECIES}.genome.fa"
			ANNOTATION_FILE="${DIR_CURRENT_SPECIES}/${SPECIES}${INSECTBASE_LINK_END_ANNOTATION}"
			CDS_FILE="${DIR_CURRENT_SPECIES}/${SPECIES}${INSECTBASE_LINK_END_CDS}"
			mkdir ${DIR_CURRENT_SPECIES}
			
			# Download genome, annotation and CDS files in parallel
			echo "Downloading files for ${SPECIES}"
			curl -L ${INSECTBASE_LINK_START}${SPECIES}-${SPECIES}${INSECTBASE_LINK_END_GENOME} -o ${FILE_GENOME_DOWNLOAD} &
			curl -L ${INSECTBASE_LINK_START}${SPECIES}-${SPECIES}${INSECTBASE_LINK_END_CDS} -o ${CDS_FILE} &
			curl -L ${INSECTBASE_LINK_START}${SPECIES}-${SPECIES}${INSECTBASE_LINK_END_ANNOTATION} -o ${ANNOTATION_FILE}

			# Waits for all downloads to be finished
			wait
			# unpack genome
			tar -xjf ${FILE_GENOME_DOWNLOAD} -C ${DIR_CURRENT_SPECIES}

			# Check if tar file contained a directory witht he genome instead of just the genome file.
			#  If TRUE, then the file is moved to the current species directory
			if [ -d ${DIR_CURRENT_SPECIES}/${SPECIES} ]; then
				mv "${DIR_CURRENT_SPECIES}/${SPECIES}/${SPECIES}.genome.fa" ${GENOME_FILE}
				rm -r ${DIR_CURRENT_SPECIES}/${SPECIES}
			fi

			# Go through every genome and check that the file actually exists
			if [ -f "${GENOME_FILE}" ]; then
				echo "Modifying genome ${SPECIES}" | tee -a ${LOG_FILE}
				# Modify genome file to add accession to first string of the header
				awk -v accession=${SPECIES} '{
					if (/^>/) {
						sub(/[ \t]/, ":" accession " ", $0)
					}; print $0
				}' ${GENOME_FILE} > temp_file.fa && mv temp_file.fa ${GENOME_FILE}
				# Copy genome to collect fasta file
				cat ${GENOME_FILE} >> ${PATH_COLLECTED_FASTA}
			else
				echo "Genome file belonging to ${SPECIES} seems to missing" | tee -a ${LOG_FILE}
			fi

			if [ -f "${ANNOTATION_FILE}" ]; then
				# Check if there was an error with downloading the CDS file (first line is <!DOCTYPE html>)
				if grep -q "^<!DOCTYPE html>" "${ANNOTATION_FILE}"; then
					echo "Error in annotation file of ${SPECIES}: File begins with: <!DOCTYPE html>" | tee -a ${LOG_FILE}
					echo "File will be removed" | tee -a ${LOG_FILE}
					rm ${ANNOTATION_FILE}
				else
					echo "Annotation found for ${SPECIES}"
					ANNOTATION_FOUND='TRUE'
					OUTPUT_ANNOTATION_RELATIVE="./${current_group}/annotations/${SPECIES}.gff3"
				fi
			else
				echo "No annotation found for ${SPECIES}" | tee -a ${LOG_FILE}
			fi

			if [ -f "${CDS_FILE}" ]; then
				# Check if there was an error with downloading the CDS file (first line is <!DOCTYPE html>)
				if grep -q "^<!DOCTYPE html>" "${CDS_FILE}"; then
					echo "Error in CDS file of ${SPECIES}: File begins with: <!DOCTYPE html>" | tee -a ${LOG_FILE}
					echo "File will be removed" | tee -a ${LOG_FILE}
					rm ${CDS_FILE}
				else
					echo "CDS found for ${SPECIES}" | tee -a ${LOG_FILE}
					CDS_FOUND='TRUE'
					OUTPUT_CDS_RELATIVE="./${current_group}/cds/${SPECIES}.cds.fna"
				fi
			else
				echo "No CDS found for ${SPECIES}" | tee -a ${LOG_FILE}
			fi

			echo -e "${SPECIES}\t${current_group}\tinsect_base\t${ANNOTATION_FOUND}\t${OUTPUT_ANNOTATION_RELATIVE}\t${CDS_FOUND}\t${OUTPUT_CDS_RELATIVE}" >> ${PATH_TO_INFORMATION}
		done
		echo -e "\n"
	done	
}	



generate_indexes() {
	# Get list of all groups in tmp dir as they were downloaded during download data 
		# (list entries | filter for dirs | reverse | get first column (dir name) | reverse again | turn newline to space in order to get an array)
	GROUPS_LIST=$(ls -l ${TMP_DIR_GROUPS} | grep "^d" | rev | cut -d ' ' -f 1 | rev | tr '\n' ' ')

	mkdir ${OUTPUT_DIR}
	cp ${PATH_TO_INFORMATION} ${OUTPUT_DIR}

	# Generate file that contains paths to DBS
	PATH_DB_INFO="${OUTPUT_DIR}/info_db_dirs.tsv"
	touch ${PATH_DB_INFO}
	echo -e "group\tpath" > ${PATH_DB_INFO}

	# Go through all accessions of the current group
	for CURRENT_GROUP in ${GROUPS_LIST[@]}; do
		TMP_DIR_CURRENT_GROUP="${TMP_DIR_GROUPS}/${CURRENT_GROUP}"
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

		# Generate and prepare file to keep information about the presence of a CDS file. Might be redundant with annotation file
		PATH_CDS_INFO="${OUTPUT_DIR_CURRENT_GROUP}/info_CDS_${CURRENT_GROUP}.tsv"
		touch ${PATH_CDS_INFO}
		echo -e "accession\tCDS\tpath" > ${PATH_CDS_INFO}

		# Go through every accession
		for CURRENT_ACCESSION in ${GENOMES_LIST[@]}; do
			TMP_DIR_CURRENT_ACCESSION="${TMP_DIR_CURRENT_GROUP}/${CURRENT_ACCESSION}"

			GENOME_FILE="${TMP_DIR_CURRENT_ACCESSION}/${CURRENT_ACCESSION}.genome.fa"
			ANNOTATION_FILE="${TMP_DIR_CURRENT_ACCESSION}/${CURRENT_ACCESSION}.gff3"
			CDS_FILE="${TMP_DIR_CURRENT_ACCESSION}/${CURRENT_ACCESSION}.cds.fa"

			# Check if annotation is present. If yes it will be renamed into the output dir
			# TODO: Check if file is correct (when file is missing a HTML report is downloaded)
			if [ -f "${ANNOTATION_FILE}" ]; then
				TARGET_FILE="${ANNOTATION_DIR}/${CURRENT_ACCESSION}.gff3"
				cp ${ANNOTATION_FILE} ${TARGET_FILE}

				echo -e "${CURRENT_ACCESSION}\tTRUE\t./${CURRENT_GROUP}/annotations/${CURRENT_ACCESSION}.gff3" >> ${PATH_ANNOTATION_INFO}
			else
				echo -e "${CURRENT_ACCESSION}\tFALSE\tNA" >> ${PATH_ANNOTATION_INFO}
			fi

			# Do the same thing as above for the CDS file
			if [ -f "${CDS_FILE}" ]; then
				TARGET_FILE="${CDS_DIR}/${CURRENT_ACCESSION}.cds.fa"
				cp ${CDS_FILE} ${TARGET_FILE}

				echo -e "${CURRENT_ACCESSION}\tTRUE\t./${CURRENT_GROUP}/cds/${CURRENT_ACCESSION}.cds.fa" >> ${PATH_CDS_INFO}
			else
				echo -e "${CURRENT_ACCESSION}\tFALSE\tNA" >> ${PATH_CDS_INFO}
			fi

				
		done
		
		# Copy fasta file containing all genome chromosomes of the current group
		FILENAME_COLLECTED_GENOMES="genomes_collected_${CURRENT_GROUP}.fna"
		FILE_COLLECTED_GENOMES="${TMP_DIR_GROUPS}/${CURRENT_GROUP}/${FILENAME_COLLECTED_GENOMES}"
		OUTPUT_FILE_COLLECTED_GENOMES="${OUTPUT_DIR_CURRENT_GROUP}/${FILENAME_COLLECTED_GENOMES}"
		if [ -f "${FILE_COLLECTED_GENOMES}" ]; then
			echo "Copying file with collected genomes to output dir: ${FILENAME_COLLECTED_GENOMES}" | tee -a ${LOG_FILE}
			cp ${FILE_COLLECTED_GENOMES} ${OUTPUT_FILE_COLLECTED_GENOMES}
		else
			echo -e "Something went wrong with the collected genomes file of ${CURRENT_GROUP}.\nThis should be investigated manually." | tee -a ${LOG_FILE}
		fi

		# Build BLASTN DB
		echo "Generating database for ${CURRENT_GROUP}" | tee -a ${LOG_FILE}
		DB_BASENAME="${DB_DIR}/${CURRENT_GROUP}"
		makeblastdb -in ${OUTPUT_FILE_COLLECTED_GENOMES} -out ${DB_BASENAME} -title ${CURRENT_GROUP} -dbtype nucl

		echo -e "${CURRENT_GROUP}\t${DB_BASENAME}" >> ${PATH_DB_INFO}
		
	done
}

	#TODO: Get md5 checksums
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



# Start functions
download_data_ncbi
download_data_insectbase
generate_indexes
