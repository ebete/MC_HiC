#!/usr/bin/env python3

import argparse
import gzip
import logging

import pysam
from Bio import SeqIO

import utils


def get_mapped_reads(samfile):
    read_ids = list()

    logging.info("Reading SAM file %s ...", samfile)
    with pysam.AlignmentFile(samfile, "r") as sam:
        for read in sam:
            if read.is_unmapped:
                logging.debug("Skipping %s (unaligned)", read.qname)
                continue

            metadata = utils.read_header_to_dict(read.qname)
            read_ids.append(utils.make_read_id(metadata))

    return set(read_ids)


def extract_mapped_read_records(fasta_file, read_ids):
    logging.info("Getting FASTA records that have mapped fragments ...")
    mappable = 0
    records = 0
    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            records += 1
            metadata = utils.read_header_to_dict(record.id)
            if utils.make_read_id(metadata) not in read_ids:
                continue
            mappable += 1
            print(record.format("fasta"), end="")
    logging.info("%d records extracted from %d total (%.1f%%).", mappable, records, mappable / records * 100)


if __name__ == '__main__':
    utils.init_logger()

    # get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="SAM/BAM file to extract mapped fragments from.", metavar="SAM",
                        action="store", type=str)
    parser.add_argument("input_fasta", help="FASTA file with all the query sequences used in generating the SAM file "
                                            "(gzipped).", metavar="FASTA", action="store", type=str)
    args = parser.parse_args()

    mapped_read_ids = get_mapped_reads(args.input_sam)
    extract_mapped_read_records(args.input_fasta, mapped_read_ids)
