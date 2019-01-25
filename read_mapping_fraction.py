#!/usr/bin/env python3

import argparse
import logging

import pysam


def get_read_coverage(samfile):
    read_coverage = dict()
    read_length = dict()
    read_contiguity = dict()

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

            read_contiguity.setdefault(rdid, list()).append(int(read_metadata["Fr.Id"]))

            if rdid not in read_length:
                read_length[rdid] = int(read_metadata["Rd.Ln"])

    get_read_contiguity(read_contiguity)
    write_coverage_per_read(read_coverage, read_length)


def get_read_contiguity(read_contiguity):
    logging.info("Calculating read contiguity ...")
    total_mapped = 0
    non_contiguous = 0
    for read_id, fragments in read_contiguity.items():
        last_frag = -1
        for frag_id in sorted(set(fragments)):
            total_mapped += 1
            if last_frag < 0:  # init with first fragment
                last_frag = frag_id
                continue

            if frag_id != (last_frag + 1):  # non-contiguous mapping
                non_contiguous += 1
            last_frag = frag_id
    logging.info("Out of %d fragments, %d were non-contiguous (%.1f%%).", total_mapped, non_contiguous,
                 non_contiguous / total_mapped * 100)


def write_coverage_per_read(read_coverage, read_length):
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
