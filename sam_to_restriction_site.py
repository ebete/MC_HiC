#!/usr/bin/env python3

import argparse
import logging
import pickle

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pysam
import seaborn as sns

plt.style.use("ggplot")
matplotlib.rcParams['figure.figsize'] = (16, 8)
matplotlib.rcParams['figure.dpi'] = 80
# matplotlib.rcParams["figure.autolayout"] = True
sns.set()


def load_site_index(fname):
    logging.info("Loading restriction site index %s ...", fname)
    with open(fname, "rb") as f:
        return pickle.load(f)


def get_closest_site(site_index, cut_sites, locus):
    # fast-forward to site closest to locus
    dist_to_locus = cut_sites[site_index] - locus
    while dist_to_locus < 0:
        site_index += 1
        dist_to_locus = cut_sites[site_index] - locus
    # step back if previous site is closer
    if locus - cut_sites[site_index - 1] < dist_to_locus:
        site_index -= 1
        dist_to_locus = cut_sites[site_index] - locus
    return site_index, dist_to_locus


def get_sites_per_read(input_index, input_sam, distance_cutoff=20, sam_region="."):
    logging.info("Finding closest sites for mapped reads in %s ...", input_sam)
    cut_sites = input_index[sam_region.split(":")[0]][0]  # TODO: less hack-y selection of chromosome
    sites_per_read = {}
    with pysam.AlignmentFile(input_sam, "r") as samfile:
        start_site = 1
        for read in samfile.fetch(region=sam_region):
            # check alignment start
            start_site, start_dist = get_closest_site(start_site, cut_sites, read.reference_start)
            if start_dist > distance_cutoff:
                start_site += 1
            # check alignment end
            end_site = start_site
            end_site, end_dist = get_closest_site(end_site, cut_sites, read.reference_end)
            if end_dist > distance_cutoff:
                end_site -= 1
            # assign found sites to read
            site_indices = [cut_sites[i] for i in range(start_site, end_site + 1)]
            if len(site_indices) > 0:
                sites_per_read.setdefault(read.qname, []).append(site_indices)
    return sites_per_read


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s]: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_index", help="Restriction site index file", metavar="INDEX", action="store", type=str)
    parser.add_argument("input_sam", help="Sorted SAM/BAM file to process", metavar="SAM", action="store", type=str)
    parser.add_argument("-r", "--region", help="Limit search to specific region", metavar="REGION", action="store",
                        type=str, default=".")
    parser.add_argument("-b", "--bin-size", help="Size of the bins", metavar="SIZE", action="store",
                        type=int, default=5000)
    args = parser.parse_args()

    enzyme_sites = load_site_index(args.input_index)
    sites_per_read = get_sites_per_read(enzyme_sites, args.input_sam, 20, args.region)
    del enzyme_sites

    site_interactions = []
    for k, v in sites_per_read.items():
        for i in range(len(v)):
            for j in range(len(v)):
                if i == j:
                    continue
                site_interactions.append([v[i][0], v[j][0]])
    del sites_per_read
    scatter_points = np.array(site_interactions)
    del site_interactions

    plt.scatter(scatter_points[:, 1], scatter_points[:, 1])
    plt.show()
