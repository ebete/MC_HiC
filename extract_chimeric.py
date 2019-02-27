#!/usr/bin/env python3

import argparse
import logging
import sys
from tempfile import NamedTemporaryFile

import pysam

import utils


def parse_chimeric(sam_input, bam_output, min_mappings):
    with NamedTemporaryFile("wb", buffering=0) as tmp:
        logging.info("Created temporary BAM file %s", tmp.name)
        with pysam.AlignmentFile(sam_input) as sam, \
            pysam.AlignmentFile(tmp, "wb", template=sam) as fout:
            grouped_reads = dict()

            logging.info("Reading SAM file ...")
            for read in sam:
                # fragment id (Fr.Id) should be ignored
                read_metadata = utils.read_header_to_dict(read.qname)
                rdid = utils.make_read_id(read_metadata)

                grouped_reads.setdefault(rdid, list()).append(read)

            logging.info("Writing chimeric reads to temporary BAM")
            for read_id, reads in grouped_reads.items():
                if len(reads) < min_mappings:
                    continue
                for x in reads:
                    fout.write(x)

        # write to sorted BAM
        logging.info("Writing final BAM to %s ...", bam_output)
        pysam.sort("-o", bam_output, tmp.name)


if __name__ == "__main__":
    utils.init_logger()

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="SAM", action="store", nargs="?",
                        type=argparse.FileType("r"), default=sys.stdin)
    parser.add_argument("output_sam", help="Output BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("-m", "--min-mappings", metavar="N", help="Minimum number of mapping locations of a read.",
                        action="store", type=int, default=1)
    args = parser.parse_args()

    parse_chimeric(args.input_sam, args.output_sam, args.min_mappings)
    args.input_sam.close()

    logging.shutdown()
