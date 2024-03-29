#!/usr/bin/env python3

import argparse
import logging
import pickle

import pysam

import utils


def load_site_index(fname):
    logging.info("Loading restriction site index %s ...", fname)
    with open(fname, "rb") as f:
        return pickle.load(f)


def get_closest_sites(sam_file, enzyme_index):
    """
    Find the closest restriction site to each alignment.

    :param sam_file: SAM file to take the alignments from.
    :param enzyme_index: Restriction enzyme dictionary containing the positions.
    """
    print("read_id", "start_dist", "end_dist", sep='\t')

    with pysam.AlignmentFile(sam_file, "r") as sam:
        start_idx = 1
        last_chr = ""
        enzyme_pos = []
        for read in sam.fetch():
            if read.flag & 0x4 != 0:
                logging.debug("Skipping %s (unaligned)", read.qname)
                continue

            if last_chr != read.reference_name:
                last_chr = read.reference_name
                start_idx = 1
                enzyme_pos = enzyme_index[read.reference_name][0]

            start_dst, start_idx = get_distance_to_site(read.reference_start, enzyme_pos, start_idx)
            start_idx -= 1
            end_dst, end_idx = get_distance_to_site(read.reference_end, enzyme_pos, start_idx)

            print(read.qname, start_dst, end_dst, sep='\t')


def get_distance_to_site(site_pos, enzyme_positions, start_idx=1):
    """
    Compute the shortest distance from site_pos to a restriction site.

    :param site_pos: Position to perform the lookup for.
    :param enzyme_positions: Dictionary containing the positions of restriction enzymes.
    :param start_idx: Start looking at this restriction site index (reduced search space).
    :return: The distance and index of the closest restriction site.
    """
    idx = start_idx if start_idx > 1 else 1
    while idx < len(enzyme_positions) - 1 and (site_pos - enzyme_positions[idx]) > 0:
        idx += 1

    if abs(site_pos - enzyme_positions[idx]) < site_pos - enzyme_positions[idx - 1]:
        return site_pos - enzyme_positions[idx], idx
    else:
        return site_pos - enzyme_positions[idx - 1], idx


if __name__ == '__main__':
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_index", help="Restriction site index file.", metavar="INDEX", action="store", type=str)
    parser.add_argument("input_sam", help="SAM/BAM file.", metavar="SAM", action="store", type=str)
    args = parser.parse_args()

    enzyme_sites = load_site_index(args.input_index)

    # load and sort mapped fragments
    logging.info("Finding closest sites for mapped reads in %s ...", args.input_sam)
    get_closest_sites(args.input_sam, enzyme_sites)
