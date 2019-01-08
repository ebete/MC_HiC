#!/usr/bin/env nextflow

/* Copyright (c) 2019 Thom Griffioen
 * MIT License
 *
 * Run the Multi-Contact Hi-C pipeline
 */

// Queue channels
raw_files = Channel
	.fromPath(params.input, checkIfExists: true)
	.map { file -> tuple(file.simpleName, file) }


// Fragment the reads
process createFragments {
	cpus 1
	memory "500MB"
	tag "${dataset}"

	input:
	set dataset, file(fq_file) from raw_files

	output:
	set dataset, file('fragments.fa.gz') into read_fragments
	
	script:
	"""
	python3 ${params.script_dir}/fastq_subsample_and_merge.py ${params.fragments} fragments.fa.gz ${fq_file}
	"""
}

// Align to reference and export alignments
process alignReads {
	label "multicore"
	memory { 4.GB * task.cpus }
	tag "${dataset}"
	publishDir "result", mode: "rellink", overwrite: true

	input:
	set dataset, reads from read_fragments

	output:
	set dataset, file("${dataset}.bam"), file("${dataset}.bam.bai") into alignment_file
	val dataset into alignment_done

	script:
	"""
	bwa mem -x ont2d -t ${task.cpus} ${params.bwa} ${params.ref} ${reads} \
	| samtools view -b \
	| samtools sort > ${dataset}.bam

	samtools index ${dataset}.bam ${dataset}.bam.bai
	"""
}

// Filter and merge alignments
process filterAndMerge {
	cpus 1
	memory "500MB"
	tag "${dataset}"

	input:
	set dataset, file(alignment), file(index) from alignment_file

	output:
	set dataset, file("positions.csv") into positions_file

	script:
	"""
	python3 ${params.script_dir}/read_map_freq.py ${params.filter} -c positions.csv ${alignment}
	"""
}

// Map reads to restriction sites
process mapToRestriction {
	cpus 1
	memory "2GB"
	tag "${dataset}"
	publishDir "result", mode: "rellink", overwrite: true

	input:
	set dataset, file(positions) from positions_file

	output:
	set dataset, file("${dataset}.csv") into interaction_file
	val dataset into interactions_done

	script:
	"""
	python3 ${params.script_dir}/sam_to_restriction_site.py ${params.restriction} ${params.ref}.gz.DpnII ${positions} ${dataset}.csv
	"""
}


// Alignment output channel
alignment_done.subscribe onNext: {
	println "Aligned ${it}"
}, onComplete: {
	println "Alignment done!"
}

interactions_done.subscribe onNext: {
	println "Mapped interactions of ${it}"
}, onComplete: {
	println "Interaction mapping done!"
}

// vim: noet:ai
