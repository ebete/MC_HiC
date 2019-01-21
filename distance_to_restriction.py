#!/usr/bin/env python3

import argparse
import logging
import pickle

import pysam


def load_site_index(fname):
    logging.info("Loading restriction site index %s ...", fname)
    with open(fname, "rb") as f:
        return pickle.load(f)


def get_closest_sites(sam_file, enzyme_index):
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
    idx = start_idx
    while idx < len(enzyme_positions) and (site_pos - enzyme_positions[idx]) > 0:
        idx += 1
    return min(abs(site_pos - enzyme_positions[idx]), site_pos - enzyme_positions[idx - 1]), idx


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s]: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_index", help="Restriction site index file.", metavar="INDEX", action="store", type=str)
    parser.add_argument("input_sam", help="SAM/BAM file.", metavar="SAM", action="store", type=str)
    args = parser.parse_args()

    enzyme_sites = load_site_index(args.input_index)
    # load and sort mapped fragments

    logging.info("Finding closest sites for mapped reads in %s ...", args.input_sam)
    get_closest_sites(args.input_sam, enzyme_sites)
