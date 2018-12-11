#!/usr/bin/env python3

import argparse
import gzip
import logging
from glob import iglob

from Bio import SeqIO, Restriction


def site_frequency(fasta_file, enzyme, max_mismatch):
    """
    Counts restriction site matches and close matches.

    :type fasta_file: str
    :param fasta_file: Path to the FASTA file to read.

    :type enzyme: str
    :param enzyme: Name of the restriction enzyme

    :type max_mismatch: int
    :param max_mismatch: Maximum number of mismatches with the restriction site

    :rtype: tuple
    :return: Tuple of a dictionary and the total length
    """
    restriction_enzyme = getattr(Restriction, enzyme)
    logging.info("Counting frequency of %s (%s) in %s; allowing %d mismatches ...",
                 restriction_enzyme.site, enzyme, fasta_file, max_mismatch)

    # slow or fast counting
    count_method = get_exact_matches if max_mismatch == 0 else get_close_matches

    total_freq = {}
    total_bp = 0
    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            d = count_method(record.seq, restriction_enzyme.site, max_mismatch)
            left_merge(total_freq, d)
            total_bp += len(record.seq)

    logging.debug("Scanned %d locations", total_bp)
    for k, v in total_freq.items():
        logging.debug("Sites with %d mismatches: %d", k, v)

    return total_freq, total_bp


def get_close_matches(seq, site, mismatch_cutoff=0):
    """
    Get the number of positions where the restriction site is found in the
    sequence with a maximum number of mismatches.

    :type seq: str
    :param seq: Sequence to look in

    :type site: str
    :param site: The restriction site to look for

    :type mismatch_cutoff: int
    :param mismatch_cutoff: Maximum number of mismatches with the restriction
        site

    :rtype: dict
    :return: A dictionary containing the number of matches for each number of
        mismatches
    """
    match_pos = {}
    # sliding window approach
    for i in range(0, len(seq) - len(site)):
        # determine hamming distance
        mismatches = 0
        for j in range(len(site)):
            if seq[i + j] != site[j]:
                mismatches += 1
                if mismatches > mismatch_cutoff:
                    break
        # add site position if criterion satisfied
        if mismatches > mismatch_cutoff:
            continue
        match_pos[mismatches] = match_pos.get(mismatches, 0) + 1
    return match_pos


def get_exact_matches(seq, site, mismatch_cutoff=0):
    return {0: seq.count(site)}


def left_merge(d1: dict, d2: dict):
    for k, v in d2.items():
        d1[k] = d1.get(k, 0) + v


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s]: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", help="Input FASTA genomes (gzipped).", metavar="INFILE", action="store",
                        type=str, nargs="+")
    parser.add_argument("-e", "--enzyme", help="Restriction enzyme to calculate cutting frequency of", metavar="ENZYME",
                        action="store", type=str, default="DpnII")
    parser.add_argument("-m", "--max-mismatch", help="Maximum number of mismatches in the cutting site",
                        metavar="MISMATCHES", action="store", type=int, default=0)
    args = parser.parse_args()

    total_freq = {}
    total_len = 0
    for glob in args.input_fasta:
        for fasta in iglob(glob):
            site_freq, fasta_length = site_frequency(fasta, args.enzyme, args.max_mismatch)
            total_len += fasta_length
            left_merge(total_freq, site_freq)

    print("Total positions: {:d}".format(total_len))
    for k, v in total_freq.items():
        print("{:s} sites with {:d} mismatches: {:d} ({:.2f}%)".format(args.enzyme, k, v, v / total_len * 100))
