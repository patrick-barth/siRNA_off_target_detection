process split_dsRNA {
	tag {dsRNA.baseName}
	publishDir "${params.output_dir}/generated_siRNAs", mode: 'copy', pattern: "${dsRNA.simpleName}_derived_siRNAs.fa"

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
	publishDir "${params.output_dir}/raw_alignments", mode: 'copy', pattern: "${siRNAs.simpleName}_alignments.bam"

	input:
	path(siRNAs)
	val(db)

	output:
	path("${siRNAs.simpleName}_alignments.bam"),	emit: siRNAs
	path("${task.process}.version.txt"), 			emit: version

	"""
	bowtie \
		-f \
		${db} \
		${siRNAs} \
		--all \
		--seedlen 8 \ 
		--seedmms 2 \
		--threads ${task.cpus} \
		--seed 0 \
		--suppress 6 \
		${siRNAs.simpleName}_alignments.bam > tmp_report.txt
	
	echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" > ${task.process}.version.txt
	"""
}