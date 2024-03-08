#!/usr/bin/env python

#############################################
# Receives a longer RNA/DNA sequence and    #
#  detects potentially functional siRNAs    #
#  that could be derived from this sequence #
#############################################

import argparse
from typing import Dict, List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



parser = argparse.ArgumentParser()
parser.add_argument('--sequence', 	'-s', type=str,	help='Sequence to derive siRNAs from')
parser.add_argument('--rules',      '-r', type=str, default='mammal', help='Rules to determine siRNAs')
parser.add_argument('--output', 	'-o', type=str,	help='Output file')
args = parser.parse_args()

########################
########################
###    parameters    ###
########################
########################
add_overhang_to_output 			= True
add_passenger_strand_to_output 	= True

siRNA_len: 				int 	= 21
siRNA_overhang_len:		int		= 2
seed_start:				int 	= 2
seed_end:				int		= 7

rules_allowed = ['mammal']

mammal_start_nuc_guide: List 	= ['A','U']
mammal_start_nuc_pass:	List 	= ['C','G']
mammal_seed: Dict = {'nuc': 	['A','U'],
					 'count': 	4,
					 'rule': 	'equal or more'}

##########################################################################################################################################

###################
###################
### main script ###
###################
###################


def main(sequence,rules,output):
	#Check if arguments are correct
	check_parameters(rules)
	
	collect_records = []

	sequence: Dict = SeqIO.parse(sequence,"fasta")
	for entry in sequence:
		# Capitalize all nucleotides
		entry.seq = entry.seq.upper()

		# All T are transformed to U
		if 'T' in entry.seq:
			print('Warning: ' + entry.id + ' contained one or several instances of the nucleotide T which have been switched to U')
			entry.seq = entry.seq.transcribe()

		# Repeat for both strands
		for strand in ['forward','reverse']:
			# If the reverse strand is being checked the reverse complement is used as sequence
			sequence = entry.seq if strand == 'forward' else entry.seq.reverse_complement()
			# Retrieve all positions that have a valid start nucleotide for siRNA guiding strands
			start_positions: List = determine_siRNA_starts(sequence,rules)
			# Check for all positions if the starting nucleotide for the siRNA passenger strand is correct
			start_positions: List = validate_siRNA_pas_start(start_positions,sequence,rules)
			# Check if seed region is correct
			start_positions: List = validate_seed_region(start_positions,sequence,rules)
			# Extract siRNAs
			
			for idx, pos in enumerate(start_positions):
				# If overhang is supposed to be added to output Ns are added to the start (Ns represent overhang)
				seq = 'N' * siRNA_overhang_len if add_overhang_to_output else ''
				# Extract actual sequence
				seq += sequence[pos:pos+siRNA_len]
				# Generate record and add it to list
				record = SeqRecord(
					seq,
					id="guide_siRNA_" + strand + '_' + entry.id + '_' + str(idx),
					name="guide_siRNA_" + strand + '_' + entry.id + '_' + str(idx),
					description='Guide strand from siRNA Nr. %s from "%s" %s strand position %s' % (idx,entry.id,strand,pos + 1)
				)
				collect_records.append(record)
				# Same for passenger strand if it's supposed to be provided in the output
				if add_passenger_strand_to_output:
					# Empty variables just to avoid potential errors
					del seq
					del record
					# Ad Passenger strand is the complement of the guide strand the complement is used here
					seq = sequence[pos-siRNA_overhang_len:pos+siRNA_len-siRNA_overhang_len].complement()
					seq += 'N' * siRNA_overhang_len if add_overhang_to_output else ''
					record = SeqRecord(
						seq,
						id="passenger_siRNA_" + strand + '_' + entry.id + '_' + str(idx),
						name="passenger_siRNA_" + strand + '_' + entry.id + '_' + str(idx),
						description='Passenger strand from siRNA Nr. %s from "%s" %s strand position %s' % (idx,entry.id,strand,pos + 1)
					)
					collect_records.append(record)
	SeqIO.write(collect_records, output, "fasta")


#######################
#######################
###    functions    ###
#######################
#######################
# Print error message and exit script	
def errx(message):
    print(message)
    exit(1)

#Check if parameters are correct
def check_parameters(rules):
	if rules not in rules_allowed:
		errx('Chosen rule \'' + rules + '\' not in allowed list of rules (' + ', '.join(rules_allowed) + ')')

# Determines start positions of the guide strand of potential siRNAs
def determine_siRNA_starts(sequence,rules):
	start_positions: List = []
	# For mammals siRNAs start with an A or U at position 1 of the guide strand
	if rules == 'mammal': 
		# Goes through every Nucleotide
		for nuc in mammal_start_nuc_guide:
			# Checks if current nucleotide is matching the start nucleotides
			#  Additionally check if the guide and passenger sequence would go out of sequence length
			start_positions.append([idx for idx, item in enumerate(sequence) if nuc == item and idx + siRNA_len + siRNA_overhang_len - 2 < len(sequence) and idx - 2 >= 0])
	# Resolve List of Lists to get just a single list and sort it
	start_positions = sum(start_positions,[])
	start_positions.sort()

	return(start_positions)

# Checks if passenger strand of potential siRNAs starts with correct nucleotide
def validate_siRNA_pas_start(positions,sequence,rules):
	if(rules == 'mammal'):
		for idx, pos in enumerate(positions):
			pas_start_nuc = Seq(sequence[pos + 18]).complement()
			# Go through every potential start and check if the passenger strand start matches the viable nucleotides. Remove the position if not
			if pas_start_nuc not in mammal_start_nuc_pass:
				del positions[idx]
	return(positions)

# Check seed region of guide strand for rule
def validate_seed_region(positions,sequence,rules):
	for idx,pos in enumerate(positions):
		# Extract the seed region
		seed_seq = sequence[pos + seed_start - 1 : pos + seed_end]
		seed_rules = ''
		# Get correct rules (Is supposed to be expanded if other domains of life have different rules for siRNAs)
		if(rules == 'mammal'):
			seed_rules = mammal_seed
		# Count amount of target nucleotides in seed region
		nuc_count = 0
		for nuc in seed_rules['nuc']:
			nuc_count += seed_seq.count(nuc)
		# Check if seed rule applies
		if seed_rules['rule'] == 'equal or more':
			if not nuc_count >= seed_rules['count']:
				del positions[idx]

		
	return(positions)


##########################
### starts main script ###
##########################
main(sequence=args.sequence,
	rules=args.rules,
	output=args.output)