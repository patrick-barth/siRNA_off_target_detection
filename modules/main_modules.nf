process build_index {

	input:
	path(dsRNA)

	output:
	path("${ref}.*"),                       emit: dsRNAs
	path("${task.process}.version.txt"), 	emit: version

	"""
	python
	
	echo -e "${task.process}\tbowtie2\t\$(bowtie2-build --version | head -1 | rev | cut -f1 -d' ' | rev)" > ${task.process}.version.txt
	"""
}