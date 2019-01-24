#!/usr/bin/env python3

import argparse
import logging

import pysam


def get_read_coverage(samfile):
    read_coverage = dict()
    read_length = dict()
    logging.info("Reading SAM file %s ...", samfile)
    with pysam.AlignmentFile(samfile, "r") as sam:
        for read in sam.fetch():
            if read.is_unmapped:
                logging.debug("Skipping %s (unaligned)", read.qname)
                continue
            read_metadata = read_header_to_dict(read.qname)
            rdid = "{}_{}".format(read_metadata["Fq.Id"], read_metadata["Rd.Id"])

            frag_start = int(read_metadata["Fr.St"])
            frag_end = int(read_metadata["Fr.En"])

            read_coverage.setdefault(rdid, list()).append((frag_start, frag_end))

            if rdid not in read_length:
                read_length[rdid] = int(read_metadata["Rd.Ln"])

    logging.info("Writing coverage per read ...")
    print("read_id", "coverage", sep=";")
    for read_id, ranges in read_coverage.items():
        # pretty inefficient, but it works
        coverage = set()
        for s, e in ranges:
            coverage |= set(range(s, e))

        covered = len(coverage) / read_length[read_id]
        print(read_id, "{:.3f}".format(covered), sep=";")


def read_header_to_dict(header):
    read_metadata = {}
    for x in str(header).split(";"):
        read_metadata[x.split(":")[0]] = x.split(":")[1]
    return read_metadata


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s]: %(message)s")

    # get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="SAM/BAM file.", metavar="SAM", action="store", type=str)
    args = parser.parse_args()

    get_read_coverage(args.input_sam)
