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



parser = argparse.ArgumentParser()
parser.add_argument('--sequence', 	'-s', type=str,	help='Sequence to derive siRNAs from')
parser.add_argument('--rules',      '-r', type=str, default='mammal', help='Rules to determine siRNAs')
args = parser.parse_args()

########################
########################
###    parameters    ###
########################
########################

mammal_start_nuc_guide: List = ['A','T','U']

##########################################################################################################################################

###################
###################
### main script ###
###################
###################


def main(sequence,rules):
	sequence: Dict = SeqIO.parse(sequence,"fasta")
	for entry in sequence:
		start_positions: List = determine_siRNA_starts(entry.seq,rules)
		print(start_positions)
		print(entry.seq)

		
			


#######################
#######################
###    functions    ###
#######################
#######################

def determine_siRNA_starts(sequence,rules):
	start_positions: List = []
	# For mammals siRNAs start with an A or U at position 1 of the guide strand
	if rules == 'mammal': 
		start_positions = [idx for idx, item in enumerate(sequence.lower()) if item in mammal_start_nuc_guide] #TODO: fix position getting

	return(start_positions)



##########################
### starts main script ###
##########################
main(sequence=args.sequence,
	rules=args.rules)