#!/usr/bin/env nextflow


/* Copyright (c) 2019 Thom Griffioen
 * MIT License
 *
 * Create FASTA files from FASTQ
 */

 // default parameters
params.input = "*.fq.gz"
params.output_dir = "./"
params.script_dir = "/home/thom/PycharmProjects/McHiC"


// Queue channels
Channel
	.fromPath(params.input, checkIfExists: true)
	.map { file -> tuple(file.simpleName, file) }
	.into { raw_files1; raw_files2 }


process createNonSplit {
	cpus 1
	memory "250MB"
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(fq_file) from raw_files1

	output:
	file "${dataset}.fa.gz"
	
	script:
	"""python3 "${params.script_dir}/fastq_subsample_and_merge.py" -s 50 "${dataset}.fa.gz" "${fq_file}" """
}

process createDigested {
	cpus 1
	memory "250MB"
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(fq_file) from raw_files2

	output:
	file "${dataset}_digested.fa.gz"
	
	script:
	"""python3 "${params.script_dir}/fastq_subsample_and_merge.py" -e DpnII -s 50 "${dataset}_digested.fa.gz" "${fq_file}" """
}

// vim: noet:ai:colorcolumn=0
