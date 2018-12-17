#!/usr/bin/env python3

import argparse
import logging
import pickle

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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


def get_sites_per_read(input_index, mapping_df, distance_cutoff=20):
    sites_per_read = {}
    start_site = 1
    for index, row in mapping_df.iterrows():
        cut_sites = input_index[row["chromosome"]][0]
        # check alignment start
        start_site, start_dist = get_closest_site(start_site, cut_sites, int(row["start"]))
        if start_dist > distance_cutoff:
            start_site += 1
        # check alignment end
        end_site = start_site
        end_site, end_dist = get_closest_site(end_site, cut_sites, int(row["end"]))
        if end_dist > distance_cutoff:
            end_site -= 1
        # assign found sites to read
        site_indices = [cut_sites[i] for i in range(start_site, end_site + 1)]
        if len(site_indices) > 0:
            sites_per_read.setdefault(row["read"], []).append(site_indices)
    return sites_per_read


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


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s]: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_index", help="Restriction site index file", metavar="INDEX", action="store", type=str)
    parser.add_argument("input_csv", help="Fragment mapping file", metavar="CSV", action="store", type=str)
    parser.add_argument("-b", "--bin-size", help="Size of the bins", metavar="SIZE", action="store",
                        type=int, default=500)
    args = parser.parse_args()

    enzyme_sites = load_site_index(args.input_index)
    # load and sort mapped fragments
    mapping_table = pd.read_csv(args.input_csv, sep=";", index_col=False)
    mapping_table.sort_values("start", axis=0, inplace=True, ascending=True)

    logging.info("Finding closest sites for mapped reads in %s ...", args.input_csv)
    sites_per_read = get_sites_per_read(enzyme_sites, mapping_table, 20)

    del enzyme_sites
    del mapping_table

    # create symmetric interaction matrix
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
    np.savetxt("{}_interactions.csv".format(args.input_csv), scatter_points, header="x;y", comments="", fmt="%d",
               delimiter=";", encoding="utf-8")

    # plot the matrix
    # fig, ax = plt.subplots(nrows=1, ncols=1)
    # ax.set_facecolor("k")
    # plt.scatter(scatter_points[:, 0], scatter_points[:, 1])
    # plt.hist2d(scatter_points[:, 0], scatter_points[:, 1], bins=args.bin_size, cmap=plt.get_cmap("inferno"), vmax=15)
    # plt.colorbar()
    # sns.jointplot(scatter_points[:, 0], scatter_points[:, 1], kind="kde", shade=True)
    #
    # plt.show()
