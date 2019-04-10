#!/usr/bin/env python3

import argparse
import gzip
import logging

from Bio import SeqIO
from Bio.Alphabet import IUPAC

import utils


def calculate_homopolymer_counts(fasta_file, min_repeats, min_count):
    """
    Count the length and occurrences of homopolymers within a FASTA file.

    :param fasta_file: FASTA file to scan.
    :param min_repeats: Minimum length of a homopolymer.
    :param min_count: Minimum observations of a homopolymer.
    """
    count_matrix = dict()

    logging.info("Extracting homopolymers from %s ...", fasta_file)
    with gzip.open(fasta_file, "rt") as handle:
        records = SeqIO.parse(handle, "fasta", IUPAC.unambiguous_dna)
        for record in records:
            last_base = ""
            poly_len = 0
            for base in record.seq:
                if base == last_base:
                    poly_len += 1
                    continue
                if poly_len >= min_repeats:
                    count_matrix.setdefault(last_base, dict())[poly_len] = \
                        count_matrix.get(last_base, dict()).get(poly_len, 0) + 1
                last_base = base
                poly_len = 1

    print("nucleotide", "polymer_length", "count", sep="\t")
    for nucleotide, polymer_lengths in count_matrix.items():
        for polymer_length, count in polymer_lengths.items():
            if count < min_count:
                continue
            print(nucleotide, polymer_length, count, sep="\t")


if __name__ == '__main__':
    utils.init_logger(logging.DEBUG)

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", help="FASTA file to count the homopolymers in (gzipped).", metavar="FASTA",
                        type=str)
    parser.add_argument("-r", "--min-repeats", help="Minimum length of a homopolymer.", metavar="REPEATS", type=int,
                        default=1)
    parser.add_argument("-c", "--min-count", help="Minimum occurrences of a homopolymer.", metavar="COUNTS", type=int,
                        default=1)
    args = parser.parse_args()

    calculate_homopolymer_counts(args.input_fasta, args.min_repeats, args.min_count)

    logging.shutdown()
