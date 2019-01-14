#!/usr/bin/env python3

import argparse
import logging

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


def parse_cigar(sam_input, sam_output, region, unclipped_cutoff, match_cutoff):
    logging.info("Reading %s and writing matching records to %s ...", sam_input, sam_output)

    with pysam.AlignmentFile(sam_input, "r") as sam, \
        pysam.AlignmentFile(sam_output, "wb", template=sam) as fout:
        for read in sam.fetch(region=region):
            align_len, matches = get_matches(read.cigar)
            seq_len, unclipped_len = unclipped_length(read.cigar)
            try:
                if (unclipped_len / seq_len) < unclipped_cutoff:
                    continue
                if (matches / align_len) < match_cutoff:
                    continue
                fout.write(read)
            except ZeroDivisionError:
                # logging.warning("%s: Skipping alignment of length zero", read.qname)
                pass


def to_cigar_string(cigar_tuple):
    cigar = ""
    for k, v in cigar_tuple:
        cigar += str(v) + cigar_decoder.get(k, " ")
    return cigar


def get_matches(cigar_tuple):
    total = 0
    matches = 0
    for k, v in cigar_tuple:
        code = cigar_decoder.get(k, " ")
        if code in "HS":  # skip clipped
            continue
        total += v
        if code in "M=X":  # (mis)matches
            matches += v
    return total, matches


def unclipped_length(cigar_tuple):
    unclipped = 0
    total = 0
    for k, v in cigar_tuple:
        code = cigar_decoder.get(k, " ")
        total += v
        if code in "HS":  # skip clipped
            continue
        unclipped += v
    return total, unclipped


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("output_sam", help="Output BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("-r", "--region", help="Limit read to specific region.", metavar="REGION",
                        action="store", type=str, default=".")
    parser.add_argument("-u", "--unclipped", help="Fraction of alignment that has to be unclipped.", metavar="CUTOFF",
                        action="store", type=float, default=0)
    parser.add_argument("-m", "--matching", help="Fraction of the alignment that has to be matches.", metavar="CUTOFF",
                        action="store", type=float, default=0)
    args = parser.parse_args()

    parse_cigar(args.input_sam, args.output_sam, args.region, args.unclipped, args.matching)

    logging.shutdown()
