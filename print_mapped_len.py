#!/usr/bin/env python3

import argparse
import logging
import sys

import pysam

cigar_decoder = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "=",
    8: "X"
}


def parse_cigar(sam_input):
    with pysam.AlignmentFile(sam_input) as sam:
        for read in sam.fetch():
            seq_len, unclipped_len = unclipped_length(read.cigar)
            print(unclipped_len)


def unclipped_length(cigar_tuple):
    unclipped = 0
    total = 0
    for k, v in cigar_tuple:
        total += v
        if cigar_decoder.get(k, "-") in "HS":  # skip clipped
            continue
        unclipped += v
    return total, unclipped


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="SAM", action="store", nargs="?",
                        type=argparse.FileType("r"), default=sys.stdin)
    args = parser.parse_args()

    parse_cigar(args.input_sam)
    args.input_sam.close()

    logging.shutdown()
