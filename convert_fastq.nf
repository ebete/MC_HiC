#!/usr/bin/env nextflow


/* Copyright (c) 2019 Thom Griffioen
 * MIT License
 *
 * Create FASTA files from FASTQ
 */

 // default parameters
params.input = "*.fq.gz"
params.enzyme = "DpnII"
params.optional_args = "-s 50 -l 2000"
params.read_min = "1500"
params.read_max = "8000"
params.output_dir = "./"
params.script_dir = "/home/thom/PycharmProjects/McHiC"


// Queue channels
raw_files = Channel
	.fromPath(params.input, checkIfExists: true)
	.map { file -> tuple(file.simpleName, file) }


process filterFastq {
	label "multicore"
	memory "500MB"
	tag "${dataset}"

	input:
	set dataset, file(fq_file) from raw_files

	output:
	set dataset, file("${dataset}") into filtered_fq1, filtered_fq2

	script:
"""
pigz -cd ${fq_file} \
| bioawk -c fastx 'length(\$seq)<=${params.read_max} && length(\$seq)>=${params.read_min} {print "@" \$name " " \$comment "\\n" \$seq "\\n+\\n" \$qual;}' \
| pigz -p ${task.cpus} > "${dataset}"
"""
}


process createNonSplit {
	cpus 1
	memory "250MB"
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(fq_file) from filtered_fq1

	output:
	file "${dataset}.fa.gz"
	
	script:
	"""python3 "${params.script_dir}/fastq_subsample_and_merge.py" ${params.optional_args} "${dataset}.fa.gz" "${fq_file}" """
}

process createDigested {
	cpus 1
	memory "250MB"
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(fq_file) from filtered_fq2

	output:
	file "${dataset}_digested.fa.gz"
	
	script:
	"""python3 "${params.script_dir}/fastq_subsample_and_merge.py" -e "${params.enzyme}" ${params.optional_args} "${dataset}_digested.fa.gz" "${fq_file}" """
}

// vim: noet:ai:colorcolumn=0
