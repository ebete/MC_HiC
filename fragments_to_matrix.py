#!/usr/bin/env python3

import argparse
import csv
import logging

import utils


def get_fragment_counts(csv_file):
    """
    Get the occurrence rate of fragments per read from a CSV file.

    :param csv_file: CSV file to extract the mapped fragments from.
    :return: Counts of fragments per read in a dictionary.
    """
    logging.info("Extracting read fragmentation from %s ...", csv_file)
    fragments_counts = dict()

    with open(csv_file, "r", newline="") as fin:
        handle = csv.reader(fin, delimiter=";")
        next(handle)
        for row in handle:
            rdid = row[0]
            fragments_counts[rdid] = fragments_counts.get(rdid, 0) + 1

    fragments_per_read = dict()
    for fragments in fragments_counts.values():
        fragments_per_read[fragments] = fragments_per_read.get(fragments, 0) + 1

    return fragments_per_read


def create_matrix(counts, limit=5):
    """
    Prints the fragment dict as a matrix to stdout.

    :param counts: Dictionary containing counts of fragments per read
    :param limit: Get the fragments per read up until this count
    """
    print("file", "\t".join([str(i) for i in range(1, limit + 1)]), sep="\t")
    for fname, fragment_count in counts.items():
        row = fname
        for i in range(1, limit + 1):
            row += '\t' + str(fragment_count.get(i, 0))
        print(row)


if __name__ == '__main__':
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", help="Input CSV files.", metavar="CSV", action="store", type=str,
                        nargs="+")
    args = parser.parse_args()

    counts_per_file = {}
    max_count = 0
    for csvfile in utils.glob_all_files(args.input_csv):
        counts = get_fragment_counts(csvfile)
        counts_per_file[csvfile] = counts
        max_count = max(max_count, max(counts.keys()))
    create_matrix(counts_per_file, max_count)

    logging.shutdown()
