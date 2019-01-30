#!/usr/bin/env python3

import argparse
import logging

import pysam


def get_read_metadata(samfile):
    reads_metadata = dict()

    logging.info("Reading SAM file %s ...", samfile)
    with pysam.AlignmentFile(samfile, "r") as sam:
        for read in sam:
            if read.is_unmapped:
                logging.debug("Skipping %s (unaligned)", read.qname)
                continue

            metadata = read_header_to_dict(read.qname)
            rdid = "{}_{}".format(metadata["Fq.Id"], metadata["Rd.Id"])

            reads_metadata.setdefault(rdid, list()).append(metadata)
    return reads_metadata


def read_header_to_dict(header):
    read_metadata = {}
    for x in str(header).split(";"):
        read_metadata[x.split(":")[0]] = x.split(":")[1]
    return read_metadata


def get_read_contiguity(metadata_collection):
    logging.info("Calculating read contiguity ...")
    total_mapped = 0
    non_contiguous = 0
    jump_sizes = dict()
    left_unmapped = 0
    right_unmapped = 0
    tails_unmapped = 0

    for read_id, read_metadata in metadata_collection.items():
        read_length = int(read_metadata[0]["Rd.Ln"])
        first_aligned = read_length + 1
        last_aligned = -1
        fragments = []
        for x in read_metadata:
            fragments.append(int(x["Fr.Id"]))
            s = int(x["Fr.St"])
            e = int(x["Fr.En"])
            first_aligned = s if s < first_aligned else first_aligned
            last_aligned = e if e > last_aligned else last_aligned

        fragments = sorted(set(fragments))
        total_mapped += len(fragments)

        if first_aligned > 50:
            left_unmapped += 1
        if read_length - last_aligned > 50:
            right_unmapped += 1
        if first_aligned > 50 and read_length - last_aligned > 50:
            tails_unmapped += 1

        last_frag = -1
        for frag_id in fragments:
            if last_frag < 0:  # init with first fragment
                last_frag = frag_id
                continue

            jump_size = frag_id - (last_frag + 1)
            if jump_size > 0:  # non-contiguous mapping
                non_contiguous += 1
                jump_sizes[jump_size] = jump_sizes.get(jump_size, 0) + 1
            last_frag = frag_id

    logging.info("Out of %d mapped fragments, %d were non-contiguous (%.1f%%).", total_mapped, non_contiguous,
                 non_contiguous / total_mapped * 100)
    total_unmapped = sum([k * v for k, v in jump_sizes.items()])
    logging.info("Number of unmapped fragments between mapped fragments: %d", total_unmapped)
    for jump in sorted(jump_sizes.keys()):
        logging.info("%s unmapped: %s (%.1f%%)", jump, jump_sizes[jump], jump_sizes[jump] / total_unmapped * 100)

    num_reads = len(metadata_collection)
    logging.info("Out of %d reads, %d had an unaligned start (%.1f%%).", num_reads, left_unmapped,
                 left_unmapped / num_reads * 100)
    logging.info("Out of %d reads, %d had an unaligned end (%.1f%%).", num_reads, right_unmapped,
                 right_unmapped / num_reads * 100)
    logging.info("Out of %d reads, %d had both the start and end unaligned (%.1f%%).", num_reads, tails_unmapped,
                 tails_unmapped / num_reads * 100)


def write_coverage_per_read(metadata_collection):
    logging.info("Writing coverage per read ...")
    print("read_id", "coverage", sep=";")
    for read_id, read_metadata in metadata_collection.items():
        # pretty inefficient, but it works
        coverage = set()
        for fragment_metadata in read_metadata:
            s = int(fragment_metadata["Fr.St"])
            e = int(fragment_metadata["Fr.En"])
            coverage |= set(range(s, e))

        covered = len(coverage) / int(read_metadata[0]["Rd.Ln"])
        print(read_id, "{:.3f}".format(covered), sep=";")


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s]: %(message)s")

    # get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="SAM/BAM file.", metavar="SAM", action="store", type=str)
    args = parser.parse_args()

    metadata = get_read_metadata(args.input_sam)

    get_read_contiguity(metadata)
    write_coverage_per_read(metadata)
