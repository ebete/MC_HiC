#!/usr/bin/env python3

import argparse
import logging

import pysam

import utils


def parse_cigar(sam_input, sam_output, region, unclipped_cutoff, match_cutoff):
    """
    Write the alignments that pass the length/matches cutoff to a new SAM file.

    :param sam_input: Input SAM file.
    :param sam_output: Output filtered SAM file.
    :param region: Limit to alignments within this region.
    :param unclipped_cutoff: Minimum unclipped length of alignments.
    :param match_cutoff: Minimum fraction of matches in a read.
    """
    logging.info("Reading %s and writing matching records to %s ...", sam_input, sam_output)

    with pysam.AlignmentFile(sam_input, "r") as sam, \
        pysam.AlignmentFile(sam_output, "wb", template=sam) as fout:
        for read in sam.fetch(region=region):
            align_len, matches = get_matches(read.cigar)
            seq_len, unclipped_len = unclipped_length(read.cigar)
            try:
                if unclipped_len < unclipped_cutoff:
                    continue
                if (matches / align_len) < match_cutoff:
                    continue
                fout.write(read)
            except ZeroDivisionError:
                # logging.warning("%s: Skipping alignment of length zero", read.qname)
                pass


def get_matches(cigar_tuple):
    """
    Count the number of matches in a CIGAR tuple.

    :param cigar_tuple: CIGAR tuple from a pysam alignment.
    :return:
    """
    total = 0
    matches = 0
    for k, v in cigar_tuple:
        code = cigar_decoder.get(k, "-")
        if code in "HS":  # skip clipped
            continue
        total += v
        if code in "M=X":  # (mis)matches
            matches += v
    return total, matches


def unclipped_length(cigar_tuple):
    """
    Number of bases in the CIGAR that are unclipped.

    :param cigar_tuple: CIGAR tuple from a pysam alignment.
    :return:
    """
    unclipped = 0
    total = 0
    for k, v in cigar_tuple:
        code = cigar_decoder.get(k, "-")
        total += v
        if code in "HS":  # skip clipped
            continue
        unclipped += v
    return total, unclipped


if __name__ == "__main__":
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("output_sam", help="Output BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("-r", "--region", help="Limit read to specific region.", metavar="REGION",
                        action="store", type=str, default=".")
    parser.add_argument("-u", "--unclipped", help="Minimum length of unclipped part of the alignment.",
                        metavar="CUTOFF",
                        action="store", type=float, default=0)
    parser.add_argument("-m", "--matching", help="Fraction of the alignment that has to be matches.", metavar="CUTOFF",
                        action="store", type=float, default=0)
    args = parser.parse_args()

    parse_cigar(args.input_sam, args.output_sam, args.region, args.unclipped, args.matching)

    logging.shutdown()
