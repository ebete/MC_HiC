#!/usr/bin/env python3

import argparse
import csv
import logging

import utils


def calculate_2mer_freq(counts_file):
    count_matrix = dict()

    with open(counts_file, "r", newline="") as handle:
        records = csv.reader(handle, delimiter="\t")
        next(records)
        for row in records:
            nuc1 = str(row[0][0])
            nuc2 = str(row[0][1])
            count = int(row[1])

            left = "x{}".format(nuc2)
            right = "{}x".format(nuc1)

            count_matrix.setdefault(nuc1, dict())[left] = count
            count_matrix.setdefault(nuc2, dict())[right] = count

    lines = ""
    header = ""
    for ref, d in count_matrix.items():
        lines += ref
        for other in sorted(d.keys()):
            lines += "\t" + str(d[other])
        lines += "\n"
        header = "x\t{}\n".format("\t".join(sorted(d.keys())))
    print(header + lines)


if __name__ == '__main__':
    utils.init_logger(logging.DEBUG)

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", help="CSV file containing the 2-mer counts.", metavar="CSV",
                        action="store", type=str)
    args = parser.parse_args()

    calculate_2mer_freq(args.input_csv)
