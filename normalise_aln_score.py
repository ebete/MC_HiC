#!/usr/bin/env python3

import argparse
import logging

import pysam

import utils


def get_normalised_alignment_score(sam_input):
    """
    Calculate an alignment score that has been normalised based on alignment length.

    :param sam_input: SAM file to read.
    """
    logging.info("Reading file %s ...", sam_input)

    print("qname", "mapq", "norm_aln_score", sep="\t")
    with pysam.AlignmentFile(sam_input, "r") as sam:
        for read in sam:
            total_length, aligned_length = unclipped_length(read.cigar)
            aln_score = read.get_tag("AS")
            print(read.qname, read.mapq, "{:.2f}".format(aln_score / aligned_length), sep="\t")


def unclipped_length(cigar_tuple):
    unclipped = 0
    total = 0
    for k, v in cigar_tuple:
        total += v
        if k in (4, 5):  # skip clipped
            continue
        unclipped += v
    return total, unclipped


if __name__ == "__main__":
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="SAM", action="store", type=str)
    args = parser.parse_args()

    get_normalised_alignment_score(args.input_sam)

    logging.shutdown()
