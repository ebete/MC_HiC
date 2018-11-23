#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
import random

from Bio import SeqIO

random.seed(1337)


def subsample_fasta(fasta_file, out_dir, n):
    """
    Splits a FASTA file into n FASTA files, randomly distributing the records
    in a balanced fashion.

    :param fasta_file: Path to the input FASTA file.
    :param out_dir: Path to the directory where the output FASTA files should be
        written to.
    :param n: Number of splits to make.
    """
    # read sequences into memory
    records = []
    record_bin = []
    i = 0
    logging.info("Reading records from %s ...", fasta_file)
    with gzip.open(fasta_file, "rt") as fin:
        for record in SeqIO.parse(fin, "fasta"):
            records.append(record)
            record_bin.append(i)
            i = (i + 1) % n

    # generate file handles
    logging.info("Randomly distributing records over %d bins ...", n)
    outfile_name = "sub_{}_" + os.path.basename(fasta_file)
    handles = [gzip.open(os.path.join(out_dir, outfile_name.format(i)), "wt") for i in range(n)]
    random.shuffle(record_bin)

    # subsample records
    logging.info("Writing records to %s ...", out_dir)
    for bin, record in zip(record_bin, records):
        SeqIO.write(record, handles[bin], "fasta")

    # close file handles
    logging.info("Done! Closing file handles ...")
    for x in handles:
        x.close()


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", help="Input FASTA file", metavar="FASTA", action="store", type=str)
    parser.add_argument("output_dir", help="Output directory", metavar="DIR", action="store", type=str)
    parser.add_argument("--subsamples", "-n", help="Number of subsamples to generate", metavar="n", action="store",
                        type=int, default=5)
    args = parser.parse_args()

    subsample_fasta(args.input_fasta, args.output_dir, args.subsamples)

    logging.shutdown()
