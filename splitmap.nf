#!/usr/bin/env nextflow


/* Copyright (c) 2019 Thom Griffioen
 * MIT License
 *
 * Create FASTA files from FASTQ
 */

 // default parameters
params.input = "*.fa.gz"
params.output_dir = "./"
params.script_dir = "/home/thom/PycharmProjects/McHiC"


// Queue channels
raw_files = Channel
	.fromPath(params.input, checkIfExists: true)
	.map { file -> tuple(file.simpleName, file) }
extend_len = [ 50, 100, 150, 200 ]


// Align digested reads
process mapFragments {
	label "multicore"
	memory { 2.GB * task.cpus }
	tag "${dataset}"

	input:
	set dataset, file(fa_file) from raw_files

	output:
	set dataset, file("aligned.bam"), file("${fa_file}") into alignment_file

	script:
"""
bwa mem -x ont2d -t ${task.cpus} -k 10 -q "/data0/thom/mm9/mm9.fa" "${fa_file}" \
| samtools view -q 1 -F 260 -b \
| samtools sort -o "aligned.bam"
"""
}

// Extract all chimeric reads
process extractChimeric {
	cpus 1
	memory "500MB"
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(alignment), file(raw_reads) from alignment_file

	output:
	set dataset, file("${dataset}.bam"), file("${raw_reads}") into chimeric_file

	script:
"""
python3 "${params.script_dir}/extract_chimeric.py" -m 3 "${alignment}" "${dataset}.bam"
"""
}

// Extract the reads with a mappable fragment
process extractMappable {
	label "multicore"
	memory "500MB"
	tag "${dataset}"

	input:
	set dataset, file(chimeric), file(raw_reads) from chimeric_file

	output:
	set dataset, file("mappable.fa.gz"), file("unmappable.fa.gz"), file("${chimeric}") into reads_excerpt

	script:
"""
python3 "${params.script_dir}/get_mapped_reads.py" ${chimeric} ${raw_reads} -m >(pigz -p ${task.cpus} > mappable.fa.gz) -u >(pigz -p ${task.cpus} > unmappable.fa.gz)
"""
}

// Run SplitMap read alignment
process makeSplitmapReads {
	label "multicore"
	memory { 2.GB * task.cpus }
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(mappable), file(unmappable), file(chimeric) from reads_excerpt
	each extension from extend_len

	output:
	file "${dataset}_${extension}.csv" into splitmap_perf

	script:
"""
python3 "${params.script_dir}/splitmap.py" -m ${extension} "${chimeric}" "${mappable}" \
| pigz -p ${task.cpus} \
> "splitmap.fa.gz"

bwa mem -x ont2d -t ${task.cpus} -k 10 -q "/data0/thom/mm9/mm9.fa" "splitmap.fa.gz" \
| samtools view -b -q 1 \
| samtools sort \
> "splitmap.bam"

python3 "${params.script_dir}/splitmap_merge.py" "${chimeric}" "splitmap.bam" \
> "${dataset}_${extension}.csv"
"""
}

// vim: noet:ai:colorcolumn=0
