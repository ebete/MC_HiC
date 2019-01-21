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
	memory "250MB"
	tag "${dataset}"

	input:
	set dataset, file(fq_file) from raw_files

	output:
	set dataset, file("fragments.fa.gz") into read_fragments
	
	script:
	"""python3 "${params.script_dir}/fastq_subsample_and_merge.py" ${params.fragments} fragments.fa.gz "${fq_file}" """
}

// Align to reference and export alignments
process alignReads {
	label "multicore"
	memory { 2.GB * task.cpus }
	tag "${dataset}"

	input:
	set dataset, reads from read_fragments

	output:
	set dataset, file("aligned.bam"), file("aligned.bam.bai") into alignment_file
	val dataset into alignment_done

	script:
"""
bwa mem -t ${task.cpus} ${params.bwa} "${params.ref}" "${reads}" \
| samtools view -q 1 -F 260 -u \
| samtools sort -l 6 -o "aligned.bam"

samtools index "aligned.bam" "aligned.bam.bai"
"""
}

// Filter based on CIGAR string
process cigarFilter {
	cpus 1
	memory "100MB"
	tag "${dataset}"

	input:
	set dataset, file(alignment), file(index) from alignment_file

	output:
	set dataset, file("cigar.bam"), file("cigar.bam.bai") into cigarfiltered_file

	script:
"""
python3 "${params.script_dir}/cigar_parse.py" ${params.cigar} "${alignment}" "cigar.bam"

samtools index "cigar.bam" "cigar.bam.bai"
"""
}

// Extract all chimeric reads
process extractChimeric {
	cpus 1
	memory "500MB"
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(alignment), file(index) from cigarfiltered_file

	output:
	set dataset, file("${dataset}.bam"), file("${dataset}.bam.bai") into chimeric_file

	script:
"""
python3 "${params.script_dir}/extract_chimeric.py" -m 3 "${alignment}" "${dataset}.bam"

samtools index "${dataset}.bam" "${dataset}.bam.bai"
"""
}

// Filter and merge alignments
process filterAndMerge {
	cpus 1
	memory "5GB"
	tag "${dataset}"

	input:
	set dataset, file(alignment), file(index) from chimeric_file

	output:
	set dataset, file("positions.csv") into positions_file

	script:
	"""python3 "${params.script_dir}/read_map_freq.py" ${params.filter} -c "positions.csv" ${alignment}"""
}

// Map reads to restriction sites
process mapToRestriction {
	cpus 1
	memory "6GB"
	tag "${dataset}"
	publishDir "${params.output_dir}", mode: "copy", overwrite: true

	input:
	set dataset, file(positions) from positions_file

	output:
	set dataset, file("${dataset}.csv") into interaction_file
	val dataset into interactions_done

	script:
	"""python3 "${params.script_dir}/sam_to_restriction_site.py" ${params.restriction} "${params.ref}.gz.DpnII" "${positions}" "${dataset}.csv" """
}

// Generate statistics
process interactionStatistics {
	cpus 1
	memory "50MB"
	tag "${dataset}"

	input:
	set dataset, file(interactions) from interaction_file

	output:
	stdout interaction_stats

	script:
"""
#!/usr/bin/env python3

fname = "${interactions}"
with open(fname, "r") as f:
	next(f)
	trans = 0
	total = 0
	for line in f:
		cols = line.strip().split(';')
		if cols[0] != cols[2]:
			trans += 1
		total += 1
	print(f"[{fname}]:")
	print(f"\\x1b[97;41minteractions\\x1b[0m: {total}")
	print(f"\\x1b[97;41mtrans\\x1b[0m: {trans} ({trans/total*100:.1f}%)")
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

interaction_stats.subscribe onNext: {
	println "${it}"
}

// vim: noet:ai:colorcolumn=0
