#!/usr/bin/env python3

import argparse
import glob
import logging
import os

import pysam


def print_mq_in_sam(samfile):
    input_name = os.path.basename(samfile).partition(".")[0]
    print(input_name, end="")
    with pysam.AlignmentFile(samfile, "r") as f:
        logging.info("Parsing MQ from %s ...", samfile)
        for read in f.fetch():
            print(";{}".format(read.mapq), end="")
    print()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM files.", metavar="INFILE", action="store", type=str,
                        nargs="+")
    args = parser.parse_args()

    for globfile in args.input_sam:
        for samfile in glob.iglob(globfile):
            print_mq_in_sam(samfile)
