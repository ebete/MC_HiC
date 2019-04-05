#!/usr/bin/env python3

import argparse
import logging
import os

import pysam

import utils


def print_mq_in_sam(samfile):
    """
    Write the MAPQ values of the reads to a CSV format.

    :param samfile: Input SAM file.
    """
    input_name = os.path.basename(samfile).partition(".")[0]
    print(input_name, end="")
    with pysam.AlignmentFile(samfile, "r") as f:
        logging.info("Parsing MQ from %s ...", samfile)
        for read in f.fetch():
            print(";{}".format(read.mapq), end="")
    print()


if __name__ == "__main__":
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM files.", metavar="INFILE", action="store", type=str,
                        nargs="+")
    args = parser.parse_args()

    for samfile in utils.glob_all_files(args.input_sam):
        print_mq_in_sam(samfile)

    logging.shutdown()
