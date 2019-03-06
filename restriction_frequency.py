#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
import pickle
import re

from Bio import SeqIO, Restriction

import utils


def site_frequency(fasta_file, enzyme, max_mismatch):
    """
    Counts restriction site matches and close matches.

    :type fasta_file: str
    :param fasta_file: Path to the FASTA file to read.

    :type enzyme: str
    :param enzyme: Name of the restriction enzyme

    :type max_mismatch: int
    :param max_mismatch: Maximum number of mismatches with the restriction site

    :rtype: dict
    :return: Dictionary containing the restriction site positions for each
        record.
    """
    restriction_enzyme = getattr(Restriction, enzyme)
    logging.info("Finding indices of %s (%s) in %s; allowing %d mismatches ...",
                 restriction_enzyme.site, enzyme, fasta_file, max_mismatch)

    site_locations = {}
    with gzip.open(fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            site_locations[record.id] = get_close_matches(record.seq, restriction_enzyme.site, max_mismatch)

            logging.debug("[%s] Scanned %d locations", record.id, len(record.seq))
            for k, v in site_locations[record.id].items():
                logging.debug("[%s] Sites with %d mismatches: %d", record.id, k, len(v))

    return site_locations


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
    if mismatch_cutoff == 0:  # much faster than the following approach
        return {0: [m.start() for m in re.finditer("(?={:s})".format(site), str(seq), re.IGNORECASE)]}

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
        match_pos.setdefault(mismatches, []).append(i)
    return match_pos


def export_site_index(fname, enzyme_name="", overwrite=False):
    if os.path.exists(fname) and not overwrite:
        logging.warning("File %s exists; Index will not be updated.", fname)
        return

    logging.info("Creating %s index file %s ...", enzyme_name, fname)
    with open(fname, "wb") as index_file:
        pickle.dump(site_freq, index_file, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    utils.init_logger()

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", help="Input FASTA genomes (gzipped).", metavar="INFILE", action="store",
                        type=str, nargs="+")
    parser.add_argument("-e", "--enzyme", help="Restriction enzyme to calculate cutting frequency of. Can be a comma-"
                                               "separated list of multiple enzymes.", metavar="ENZYME", action="store",
                        type=str, default="DpnII")
    parser.add_argument("-o", "--output-index", help="Create restriction site index files.", action="store_true")
    parser.add_argument("-f", "--overwrite", help="Overwrite existing index files.", action="store_true")
    parser.add_argument("-m", "--max-mismatch", help="Maximum number of mismatches in the cutting site.",
                        metavar="MISMATCHES", action="store", type=int, default=0)
    args = parser.parse_args()

    enzyme_list = args.enzyme.split(',')

    for fasta in utils.glob_all_files(args.input_fasta):
        for enzyme in enzyme_list:
            site_freq = site_frequency(fasta, enzyme, args.max_mismatch)

            if not args.output_index:
                continue
            fname = "{:s}.{:s}".format(fasta, enzyme)
            export_site_index(fname, enzyme, args.overwrite)
