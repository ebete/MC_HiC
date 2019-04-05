# Multi-Contact 4C Project
This directory contains most of the scripts used during my internship at the Hubrecht Institute, Utrecht.
This readme contains a quick description of what each script does.


## [align_multiconfig.py](align_multiconfig.py)
Execute all alignment configurations in [config.csv](config.csv).


## [bowtie_scoring.R](bowtie_scoring.R)
Plot the `bowtie2` scoring functions.


## [calculate_2mer_freq.py](calculate_2mer_freq.py)
Convert the 2-mer count generated by Jellyfish to a CSV format.


## [cigar_parse.py](cigar_parse.py)
Write the alignments that pass the length/matches cutoff to a new SAM file.


## [collapse_coverage.py](collapse_coverage.py)
Convert the `samtools depth` output to a range-based list of covered regions.


## [compare_merged_with_original.py](compare_merged_with_original.py)
Print statistics of alignments that occur in both SAM files to a CSV format.


## [convert_fastq.nf](convert_fastq.nf)
Nextflow workflow that consumes raw FASTQ files and creates filtered/digested fragments in FASTA format.


## [count_aligned_bases.py](count_aligned_bases.py)
Sums the number of nucleotides in all reads that were used in alignment.


## [create_heatmap.R](create_heatmap.R)
Plot a heatmap from the CSV file generated by [fragments_to_matrix.py](fragments_to_matrix.py).


## [create_reference.sh](create_reference.sh)
Generate a `bowtie2`, `bwa`, and `last` database from a genome.


## [demux_sam.py](demux_sam.py)
Demultiplex a SAM file. This will separate reads based on the source file.


## [distance_to_restriction.py](distance_to_restriction.py)
Compute the distance between alignment start/ends to the closest restriction site.


## [extendmap.py](extendmap.py)
ExtendMap script


## [extendmap_merge.py](extendmap_merge.py)
...


## [extendmap_mergemap.nf](extendmap_mergemap.nf)
Nextflow pipeline for executing [ExtendMap](extendmap.py)/[MergeMap](mergemap.py).


## [extract_chimeric.py](extract_chimeric.py)
Extract all reads that have at least a certain number of alignments (split-reads).


## [extract_primers_from_reads.sh](extract_primers_from_reads.sh)
Write the 50 bp tail-ends of reads to a new FASTA file.


## [fastq_subsample_and_merge.py](fastq_subsample_and_merge.py)
Convert a FASTQ to FASTA and perform digestion/filtering.


## [fragments_to_matrix.py](fragments_to_matrix.py)
Get the number of fragments per read from the [read_map_freq.py](read_map_freq.py) CSV and print it in a matrix format.


## [get_cut_region.py](get_cut_region.py)
Gets the viewpoint region from a reference genome in FASTA format based on the primer positions given.


## [get_depth.sh](get_depth.sh)
Writes the coverage of a SAM file to a BED formatted file that can be loaded in IGV.


## [get_mapped_reads.py](get_mapped_reads.py)
Split a FASTA file into a FASTA with the aligned and one with the unaligned reads.


## [get_mq_dist.py](get_mq_dist.py)
Write the MAPQ values of the reads to a CSV format.


## [get_read_mapping.py](get_read_mapping.py)
Print the aligned positions of a single read.


## [homopolymer_freq.py](homopolymer_freq.py)
Count the number of homopolymers in a FASTA file and print it in CSV format.


## [mapping_dist.R](mapping_dist.R)
Script that I used to store a bunch of random plots.


## [mergemap.py](mergemap.py)
MergeMap script.
Generate new fragments from two unmapped ones.
All adjoined unmapped fragments will be concatenated to create a single, larger fragment.


## [nextflow.config](nextflow.config)
Nextflow workflow configuration file for [pipeline.nf](pipeline.nf).


## [normalise_aln_score.py](normalise_aln_score.py)
Calculate an alignment score for each alignment that has been normalised based on alignment length.
_Not really used._


## [pipeline.nf](pipeline.nf)
Runs part of the MC-4C filtering steps.
Configured using [nextflow.config](nextflow.config).


## [plot_cutsite_dist.R](plot_cutsite_dist.R)
Plot statistics about restriction enzymes in a genome.


## [plot_interactions.R](plot_interactions.R)
Plot found interactions in a 2D plot.


## [plot_kmer_deviation.R](plot_kmer_deviation.R)
Plot the base association rates of different reads compared to the background rates.


## [plot_length_dist.R](plot_length_dist.R)
Plot the alignment length distribution of unmapped and mapped fragments superimposed.


## [plot_mq_dist.R](plot_mq_dist.R)
Plot the MAPQ values distribution of alignments.


## [print_mapped_len.py](print_mapped_len.py)
Print the unclipped length of all alignments.


## [read_map_freq.py](read_map_freq.py)
Get the alignment positions of the reads and write it to CSV.
Read fragments are merged into a single fragment if they are too close together.


## [read_mapping_fraction.py](read_mapping_fraction.py)
Calculate various statistics about the read fragment alignment contiguity (e.g. current and next fragment both mapped _i_ times out of _j_).


## [restriction_frequency.py](restriction_frequency.py)
Create a dictionary containing all positions of restriction sites on a reference genome.
This will be pickled to a file.


## [run_2mer.sh](run_2mer.sh)
Calculate 2-mer rates of a FASTA file and write it to CSV.
Uses Jellyfish.


## [sam_to_restriction_site.py](sam_to_restriction_site.py)
Create an interaction matrix from a SAM file.


## [seed_mismatch.R](seed_mismatch.R)
Plot estimates of mismatches in k-mers in Oxford Nanopore reads.


## [subsample_fasta.py](subsample_fasta.py)
Splits a FASTA file into n FASTA files, randomly distributing the records in a balanced fashion.


## [utils.py](utils.py)
Python file containing helper functions.


## [view_enzymecuts.py](view_enzymecuts.py)
Highlights all DpnII sites with at most a single point mutation.
