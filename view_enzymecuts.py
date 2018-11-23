#!/usr/bin/env python3

import gzip
import re
import argparse
from Bio import SeqIO
from Bio.Restriction import DpnII

# 'GATC' with a single substitution
patt = re.compile(r"(GAT[^C])|(GA[^T]C)|(G[^A]TC)|([^G]ATC)")


def view_sites(fasta_file):
    """
    Highlights all DpnII sites and DpnII sites with a point mutation.

    :param fasta_file: Path to the FASTA file to read.
    """
    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            print(">" + record.id)
            query = record.seq
            res = str(query) \
                .replace(DpnII.site, "\033[1;97;41m{}\033[0m".format(DpnII.site))
            res = patt.sub("\033[0;35m\\g<0>\033[0m", res)
            print(res)


if __name__ == "__main__":
    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", help="Input FASTA file (gzipped).", metavar="INFILE", action="store", type=str)
    args = parser.parse_args()

    view_sites(args.input_fasta)
