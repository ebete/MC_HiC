#!/usr/bin/env python3

import argparse
import logging
import re

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


def find_read_mappings(read_id, sam_file):
    alignments = []
    read_seq = ""
    with pysam.AlignmentFile(sam_file, "r") as sam:
        for read in sam.fetch():
            if read.flag & 0x4 != 0:
                logging.debug("Skipping %s (unaligned)", read.qname)
                continue

            read_metadata = {}
            for x in str(read.qname).split(";"):
                read_metadata[x.split(":")[0]] = x.split(":")[1]
            rdid = "{}_{}".format(read_metadata["Fq.Id"], read_metadata["Rd.Id"])
            if rdid != read_id:
                continue

            alignments.append(read)
            if read.flag & 0x800 == 0:
                read_seq = read.seq
    return alignments, read_seq


def cigar_to_read_coverage(cigar_tuple):
    cov = ""
    for k, v in cigar_tuple:
        if k in (0, 1, 7, 8):  # aligned
            cov += cigar_decoder.get(k, "#") * v
        elif k in (4, 5):  # clipped
            cov += " " * v
    return cov


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s]: %(message)s")

    # get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("read_id", help="Read ID to look for. It is formatted like {Fq.Id}_{Rd.Id}.", metavar="ID",
                        action="store", type=str)
    parser.add_argument("input_sam", help="SAM/BAM file.", metavar="SAM", action="store", type=str)
    args = parser.parse_args()

    aln, seq = find_read_mappings(args.read_id, args.input_sam)
    seq = re.sub(r"(GATC)", "\033[1;97;41m\\g<0>\033[0m", seq)
    seq = re.sub(r"(GAT[^C])|(GA[^T]C)|(G[^A]TC)|([^G]ATC)", "\033[0;97;106m\\g<0>\033[0m", seq)

    print(">{}".format(args.read_id))
    print(seq)
    for x in aln:
        coverage = cigar_to_read_coverage(x.cigar)
        print(coverage)
