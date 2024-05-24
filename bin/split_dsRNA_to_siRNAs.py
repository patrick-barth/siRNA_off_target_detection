#!/usr/bin/env python

#############################################
# Receives a longer RNA/DNA sequence, splits#
#  it into fragments of specific lengths    #
#  and provides information on how probable #
#  they are to appear as a functional siRNA #
#############################################

import argparse
from typing import Dict, List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



parser = argparse.ArgumentParser()
parser.add_argument('--sequence', 	'-s', type=str,	help='Sequence to derive siRNAs from')
parser.add_argument('--rule',      '-r', type=str, default='mammal', help='Rules to determine siRNAs')
parser.add_argument('--output_si', 	'-o', type=str,	help='Output file')
parser.add_argument('--output_score','-u', type=str,help='output file for information on how well the siRNA is')
parser.add_argument('--overhang_length','-a', type=int, default=2, help='Length of 5\' overhang ')
args = parser.parse_args()

########################
########################
###    parameters    ###
########################
########################
add_overhang_to_output 			= False

siRNA_len: 				int 	= 21
seed_start:				int 	= 2
seed_end:				int		= 7

rules_allowed = ['mammal']

rules = {'mammal':{
	'start_guide': 	{
		'nuc':['A','T'],
		'score':0.2},
	'start_pass' :	{
        'nuc':['C','T'],
        'score':0.1},
	'seed_region': {
		'pos_start': 	2,
		'pos_end':		7,
		'rules': [
			{'nuc': 	['A','T'],
			'count': 	4,
			'rule': 	'equal or more',
			'value':	0.8}
		]
	}
}}

##########################################################################################################################################

###################
###################
### main script ###
###################
###################

def main(sequence,rule,output_si,overhang_len):
	#Check if arguments are correct
	check_parameters(rule)

	# If no overhang length is provided it needs to be assumed that the length is 0
	if not overhang_len:
		overhang_len = 0

	sequence: Dict = SeqIO.parse(sequence,"fasta")
	collected_potential_siRNAs = []
	for entry in sequence:
		# Capitalize all nucleotides
		entry.seq = entry.seq.upper()
		
    # All U are transformed to T
		if 'U' in entry.seq:
			print('Warning: ' + entry.id + ' contained one or several instances of the nucleotide T which have been switched to U')
			entry.seq = entry.seq.back_transcribe()

		# extract all potential siRNAs from the forward as well as the reverse strand
		potential_siRNAs = split_sequences( entry, siRNA_len, overhang_len, add_overhang_to_output, rules[rule] )
		
		for siRNA in potential_siRNAs:

			# Checks the probability of an siRNA being valid

			# Write fasta entries for potential siRNAs
			collected_potential_siRNAs.append( SeqRecord(siRNA['seq_forward'],
									id=siRNA['origin_id'] + '_forward_' + str(siRNA['count']),
									description='siRNA generated from ' + siRNA['origin_id'] + '\'s forward strand position ' + siRNA['pos_forward'] + ", siRNA score: " + str(siRNA['si_score_forward'])))
			collected_potential_siRNAs.append( SeqRecord(siRNA['seq_reverse'],
									id=siRNA['origin_id'] + '_reverse_' + str(siRNA['count']),
									description='siRNA generated from ' + siRNA['origin_id'] + '\'s reverse strand position ' + siRNA['pos_reverse'] + ", siRNA score: " + str(siRNA['si_score_reverse'])))
			
	SeqIO.write(collected_potential_siRNAs,output_si,'fasta')


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
def check_parameters(rule):
	if rule not in rules_allowed:
		errx('Chosen rule \'' + rule + '\' not in allowed list of rules (' + ', '.join(rules_allowed) + ')')

def split_sequences(sequence_entry,si_len,overhang_len,add_overhang,rule):
	seq = sequence_entry.seq
	id	= sequence_entry.id

	# Start an empty array that will be filled with entries
	collected_sequences = []
	count = 0
	
	# Go through every possible position, disregarding the once that can not be due to missing overhang
	for pos in range(overhang_len,len(seq)-si_len + 1):
		# Overhang is always on 5'end. Ns will be added to show the overlap with the the other strand
		seq_forward			= 'N' * overhang_len + seq[pos:pos+si_len] if add_overhang else seq[pos:pos+si_len]
		seq_reverse			= 'N' * overhang_len + seq[pos - overhang_len:pos + si_len - overhang_len].reverse_complement() if add_overhang else seq[pos - overhang_len:pos + si_len - overhang_len].reverse_complement()
		forward_positions	= str(pos) + '-' + str(pos + si_len)
		reverse_positions	= str(pos - overhang_len) + '-' + str(pos + si_len - overhang_len)
		si_score_forward = validate_siRNA(seq_forward,overhang_len,rule)
		si_score_reverse = validate_siRNA(seq_reverse,overhang_len,rule)

		collected_sequences.append({	'origin_id':id,
										'count':count,
										'seq_forward': seq_forward,
										'seq_reverse': seq_reverse,
										'pos_forward':forward_positions,
										'pos_reverse':reverse_positions,
										'si_score_forward':round(si_score_forward,2),
										'si_score_reverse':round(si_score_reverse,2)})

		count += 1
	return(collected_sequences)
		
def validate_siRNA(seq,overhang_len,rules):
	score_accumulative = 0
	# If Ns are added due to overhang they are removed here
	seq = str(seq).replace('N','')

	# Check the start nucleotide of the guide strand
	if 'start_guide' in rules.keys():
		if 'nuc' in rules['start_guide'] and 'score' in rules['start_guide']:
			if seq[0] in rules['start_guide']['nuc']:
				score_accumulative += rules['start_guide']['score']
	# Check the start nucleotide of the passenger strand
	if 'start_pass' in rules.keys():
		if 'nuc' in rules['start_pass'] and 'score' in rules['start_pass']:
			# If the resulting siRNA has an overhang then the start position of the passenger strand needs to be adjusted accordingly
			if Seq(seq).reverse_complement()[overhang_len] in rules['start_pass']['nuc']:
				score_accumulative += rules['start_pass']['score']
	# Check seed region for rules
	if 'seed_region' in rules.keys():
		if 'pos_start' in rules['seed_region'] and 'pos_end' in rules['seed_region'] and 'rules' in rules['seed_region']:
			# Extract seed sequence
			seed_seq = seq[rules['seed_region']['pos_start']-1:rules['seed_region']['pos_end']]
			for rule in rules['seed_region']['rules']:
				if rule['rule'] == 'equal_or_more':
					count_occurrences = 0
					for nuc in rule['nuc']:
						count_occurrences += seed_seq.count(nuc)
					if count_occurrences >= rule['count']:
						score_accumulative += rule['value']

	return(score_accumulative)

##########################
### starts main script ###
##########################
main(sequence=args.sequence,
	rule=args.rule,
	output_si=args.output_si,
	overhang_len=args.overhang_length)