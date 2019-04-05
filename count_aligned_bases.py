#!/usr/bin/env python3

import argparse
import logging
import sys

import pysam

import utils


def count_matches(sam_input):
    """
    Get the total number of nucleotides in all reads that were used in alignment.

    :param sam_input: SAM file to take the alignments from.
    :return: Total number of bases used in aligning the fragments.
    """
    logging.info("Counting aligned bases in %s ...", sam_input.name)

    total_bases = 0
    with pysam.AlignmentFile(sam_input, "r") as sam:
        for read in sam:
            total_bases += aligned_bases(read.cigar)
    return total_bases


def aligned_bases(cigar_tuple):
    """
    Get the total number of operations that consume the query sequence.

    :param cigar_tuple: CIGAR tuple from a pysam alignment.
    :return: Number of query-consuming operations.
    """
    aligned = 0
    for k, v in cigar_tuple:
        if k in {0, 1, 7, 8}:  # query consumed
            aligned += v
    return aligned


if __name__ == "__main__":
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="SAM", nargs="?", type=argparse.FileType("r"),
                        default=sys.stdin)
    args = parser.parse_args()

    try:
        print(count_matches(args.input_sam))
    finally:
        args.input_sam.close()
        logging.shutdown()
