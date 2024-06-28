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
	publishDir "${params.output_dir}/raw_alignments", mode: 'copy', pattern: "hits_${group}.tsv"

	input:
	path(siRNAs)
	set val(group)
	path(db_base_dir)

	output:
	path("hits_${group}.tsv"),	emit: aligned_siRNAs
	path("${task.process}.version.txt"), 	emit: version

	"""
	blastn \
		-query ${siRNAs} \
		-db ${db_base_dir}/${group}/db/${group} \
		-out hits_${group}.tsv \ 
		-outfmt "6 qseqid sseqid nident mismatch gaps sstrand length qstart qend qseq sstart send sseq evalue" \
		-task blastn-short \
		-word_size 6 \
		-perc_identity 60 \
		-reward 3 \
		-penalty -4 \
		-ungapped \
		-num_threads ${task.cpus}
	
	echo -e "${task.process}\tpython\t\$(python --version | rev | cut -d' ' -f1 | rev)" > ${task.process}.version.txt
	"""
}