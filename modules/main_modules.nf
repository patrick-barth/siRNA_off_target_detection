process split_dsRNA {
	tag {dsRNA.baseName}

	input:
	path(dsRNA)

	output:
	path("${dsRNA.simpleName}_derived_siRNAs.fa"),	emit: siRNAs
	path("${task.process}.version.txt"), 			emit: version

	"""
	split_dsRNA_to_siRNAs.py \
		--sequence ${dsRNA} \
		--rule mammal \
		--output_si ${dsRNA.simpleName}_derived_siRNAs.fa \
		--overhang_length 2
	
	echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" > ${task.process}.version.txt
	"""
}

process align_siRNAs {
	tag {siRNAs.baseName}

	input:
	path(siRNAs)
	path(db)

	output:
	path("${siRNAs.simpleName}_derived_siRNAs.fa"),	emit: siRNAs
	path("${task.process}.version.txt"), 			emit: version

	"""
	bowtie \
		-x ${db} \
		-f ${siRNAs} \
		--all \
		--threads ${task.cpus} \
		--seed 0 \
		--suppress 6 \
		${siRNAs.simpleName}_alignments.bam
	
	echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" > ${task.process}.version.txt
	"""
}