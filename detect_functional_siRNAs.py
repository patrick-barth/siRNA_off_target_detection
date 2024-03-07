#!/usr/bin/env python

#############################################
# Receives a longer RNA/DNA sequence and    #
#  detects potentially functional siRNAs    #
#  that could be derived from this sequence #
#############################################

import argparse
import os
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
mammal_start_nuc_guide: List 	= ['A','T','U']
mammal_start_nuc_pass:	List 	= ['C','G']
mammal_seed: Dict = {'nuc': 	['A','T','U'],
					 'count': 	4,
					 'rule': 	'equal or more'}

#TODO: Add check if args.rules is one of allowed terms (currently only mammal)

##########################################################################################################################################

###################
###################
### main script ###
###################
###################


def main(sequence,rules,output):
	sequence: Dict = SeqIO.parse(sequence,"fasta")
	for entry in sequence:
		entry.seq = entry.seq.upper() #TODO: Change all T to U and write a short message to stdout
		# Retrieve all positions that have a valid start nucleotide for siRNA guiding strands
		start_positions: List = determine_siRNA_starts(entry.seq,rules)
		#TODO: Sort for starting positions
		# Check for all positions if the starting nucleotide for the siRNA passenger strand is correct
		start_positions: List = validate_siRNA_pas_start(start_positions,entry.seq,rules)
		# Check if seed region is correct
		start_positions: List = validate_seed_region(start_positions,entry.seq,rules)
		# Extract siRNAs
		collect_records = []
		for idx, pos in enumerate(start_positions):
			seq = 'N' * siRNA_overhang_len if add_overhang_to_output else ''
			seq += entry.seq[pos:pos+siRNA_len]
			record = SeqRecord(
				seq,
				id="siRNA" + str(idx),
				name="guide_siRNA" + str(idx),
				description='Guide strand from siRNA Nr. "%s" from "%s" position "%s"' % (idx,entry.id,pos + 1)
			)
			collect_records.append(record)
			#TODO: add passenger strand in if statement

		SeqIO.write(collect_records, output, "fasta")


#######################
#######################
###    functions    ###
#######################
#######################
# Determines start positions of the guide strand of potential siRNAs
def determine_siRNA_starts(sequence,rules):
	start_positions: List = []
	# For mammals siRNAs start with an A or U at position 1 of the guide strand
	if rules == 'mammal': 
		for nuc in mammal_start_nuc_guide:
			start_positions.append([idx for idx, item in enumerate(sequence) if nuc == item and idx + siRNA_len + siRNA_overhang_len - 2 < len(sequence) and idx - 2 >= 0])
	return(sum(start_positions,[]))

# Checks if passenger strand of potential siRNAs starts with correct nucleotide
def validate_siRNA_pas_start(positions,sequence,rules):
	if(rules == 'mammal'):
		for idx, pos in enumerate(positions):
			pas_start_nuc = Seq(sequence[pos + 18]).complement()
			if pas_start_nuc not in mammal_start_nuc_pass:
				del positions[idx]
	return(positions)

# Check seed region of guide strand for rule
def validate_seed_region(positions,sequence,rules):
	for idx,pos in enumerate(positions):
		seed_seq = sequence[pos + seed_start - 1 : pos + seed_end]
		seed_rules = ''
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