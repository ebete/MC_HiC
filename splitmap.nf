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
params.bwa = "-k 14 -A 1 -B 4 -O 6 -E 1 -w 50 -T 5 -q"
params.ref = "/data0/thom/mm9/mm9.fa"
params.mq = 1

// Queue channels
raw_files = Channel
	.fromPath(params.input, checkIfExists: true)
	.map { file -> tuple(file.simpleName, file) }
extend_len = [ 50 ]


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
bwa mem -t ${task.cpus} ${params.bwa} "${params.ref}" "${fa_file}" \
	| samtools view -q ${params.mq} -F 260 -b \
	| samtools sort -o "aligned.bam"
"""
}

// Extract all chimeric reads
process extractChimeric {
	cpus 1
	memory "5GB"
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true, pattern: "*.bam"

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
	memory "1GB"
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true, pattern: "*mappable.fa.gz"

	input:
	set dataset, file(chimeric), file(raw_reads) from chimeric_file

	output:
	set dataset, file("${dataset}_mappable.fa.gz"), file("${dataset}_unmappable.fa.gz"), file("${chimeric}") into reads_excerpt1, reads_excerpt2

	script:
"""
python3 "${params.script_dir}/get_mapped_reads.py" "${chimeric}" "${raw_reads}" -m >(pigz -p ${task.cpus} > "${dataset}_mappable.fa.gz") -u >(pigz -p ${task.cpus} > "${dataset}_unmappable.fa.gz")
"""
}

// Run SplitMap read alignment
process makeSplitmapReads {
	label "multicore"
	memory { 2.GB * task.cpus }
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(mappable), file(unmappable), file(chimeric) from reads_excerpt1
	each extension from extend_len

	output:
	file "${dataset}_${extension}_splitmap.csv" into splitmap_perf
	file "${dataset}_${extension}_splitmap.bam" into splitmap_bam
	file "${dataset}_${extension}_splitmap.fa.gz" into splitmap_fasta

	script:
"""
python3 "${params.script_dir}/splitmap.py" -m ${extension} "${chimeric}" "${mappable}" \
	| pigz -p ${task.cpus} \
	> "${dataset}_${extension}_splitmap.fa.gz"

bwa mem -t ${task.cpus} ${params.bwa} "${params.ref}" "${dataset}_${extension}_splitmap.fa.gz" \
	| samtools view -b \
	| samtools sort -o "${dataset}_${extension}_splitmap.bam"

python3 "${params.script_dir}/splitmap_merge.py" "${chimeric}" "${dataset}_${extension}_splitmap.bam" \
	> "${dataset}_${extension}_splitmap.csv"
"""
}

// Run MergeMap read alignment
process makeMergemapReads {
	label "multicore"
	memory { 2.GB * task.cpus }
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(mappable), file(unmappable), file(chimeric) from reads_excerpt2

	output:
	file "${dataset}_mergemap.csv" into mergemap_perf
	file "${dataset}_mergemap.bam" into mergemap_bam
	file "${dataset}_mergemap.fa.gz" into mergemap_fasta

	script:
"""
python3 "${params.script_dir}/mergemap.py" "${chimeric}" "${mappable}" \
	| pigz -p ${task.cpus} \
	> "${dataset}_mergemap.fa.gz"

bwa mem -t ${task.cpus} ${params.bwa} "${params.ref}" "${dataset}_mergemap.fa.gz" \
	| samtools view -b \
	| samtools sort -o "${dataset}_mergemap.bam"

samtools view "${dataset}_mergemap.bam" \
	| awk 'BEGIN{FS="\\t"; OFS="\\t"; print "read_id","mapq";} {print \$1,\$5}' \
	> "${dataset}_mergemap.csv"
"""
}

// vim: noet:ai:colorcolumn=0
