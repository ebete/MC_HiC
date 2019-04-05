#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
import sys

import pysam
from Bio import SeqIO

import utils


def get_mapped_reads(samfile):
    """
    Retrieves a set of all aligned read IDs.

    :param samfile: Input SAM file.
    :return: Set containing all aligned read IDs.
    """
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


def extract_mapped_read_records(fasta_file, read_ids, mapped_out, unmapped_out):
    """
    Split a FASTA file into a FASTA with the aligned and one with the unaligned reads.

    :param fasta_file: FASTA file used for creating the SAM file.
    :param read_ids: Read IDs that were aligned
    :param mapped_out: FASTA file to write the aligned reads to.
    :param unmapped_out: FASTA file to write the unaligned reads to.
    """
    logging.info("Getting FASTA records that have mapped fragments ...")
    mappable = 0
    records = 0
    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            records += 1
            metadata = utils.read_header_to_dict(record.id)
            if utils.make_read_id(metadata) not in read_ids:
                SeqIO.write(record, unmapped_out, "fasta")
                continue
            mappable += 1
            SeqIO.write(record, mapped_out, "fasta")
    logging.info("%d records extracted from %d total (%.1f%%).", mappable, records, mappable / records * 100)


if __name__ == '__main__':
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="SAM/BAM file to extract mapped fragments from.", metavar="SAM",
                        action="store", type=str)
    parser.add_argument("input_fasta", help="FASTA file with all the query sequences used in generating the SAM file "
                                            "(gzipped).", metavar="FASTA", action="store", type=str)
    parser.add_argument("-m", "--mapped", help="Write the mappable reads to a file instead of stdout.",
                        metavar="MAPPED", action="store", type=argparse.FileType("w"), default=sys.stdout)
    parser.add_argument("-u", "--unmapped", help="Optionally output all unmappable sequences to this file.",
                        metavar="UNMAPPED", action="store", type=argparse.FileType("w"), default=os.devnull)
    args = parser.parse_args()

    try:
        mapped_read_ids = get_mapped_reads(args.input_sam)
        extract_mapped_read_records(args.input_fasta, mapped_read_ids, args.mapped, args.unmapped)
    finally:
        args.mapped.close()
        args.unmapped.close()

    logging.shutdown()
