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
    last_chr = ""
    cut_sites = []
    for index, row in mapping_df.iterrows():
        if last_chr != row["chromosome"]:
            last_chr = row["chromosome"]
            logging.info("Assigning sites on %s ...", last_chr)
            start_site = 1
            cut_sites = input_index[last_chr][0]

        # check alignment start
        while True:
            if len(cut_sites) <= start_site:
                start_site -= 1
                break
            start_dist = cut_sites[start_site] - (int(row["start"]) - distance_cutoff)
            if start_dist >= 0:
                break
            start_site += 1

        # check alignment end
        end_site = start_site
        while True:
            if len(cut_sites) <= end_site:
                break
            end_dist = cut_sites[end_site] - (int(row["end"]) + distance_cutoff)
            if end_dist > 0:
                break
            end_site += 1
        end_site -= 1

        # assign found sites to read
        site_indices = [cut_sites[i] for i in range(start_site, end_site + 1)]
        if len(site_indices) > 0:
            sites_per_read \
                .setdefault(row["read"], {}) \
                .setdefault(last_chr, []) \
                .append(site_indices)
    return sites_per_read


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s]: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_index", help="Restriction site index file", metavar="INDEX", action="store", type=str)
    parser.add_argument("input_csv", help="Fragment mapping file", metavar="CSV", action="store", type=str)
    # parser.add_argument("-b", "--bin-size", help="Size of the bins", metavar="SIZE", action="store",
    #                     type=int, default=500)
    parser.add_argument("-d", "--distance-cutoff", help="Maximum distance between a fragment and a restriction site",
                        metavar="DISTANCE", action="store", type=int, default=20)
    args = parser.parse_args()

    enzyme_sites = load_site_index(args.input_index)
    # load and sort mapped fragments
    mapping_table = pd.read_csv(args.input_csv, sep=";", index_col=False)
    mapping_table.sort_values(by=["chromosome", "start"], axis=0, inplace=True, ascending=True)

    logging.info("Finding closest sites for mapped reads in %s ...", args.input_csv)
    sites_per_read = get_sites_per_read(enzyme_sites, mapping_table, args.distance_cutoff)

    del enzyme_sites
    del mapping_table

    # create symmetric interaction pairs
    logging.info("Exporting interacting pairs to %s_interactions.csv ...", args.input_csv)
    site_interactions = []
    for read, chr_mapping in sites_per_read.items():
        for c1, pos1 in chr_mapping.items():
            for c2, pos2 in chr_mapping.items():
                for i in range(len(pos1)):
                    for j in range(len(pos2)):
                        if c1 == c2 and i == j:
                            continue
                        site_interactions.append([c1, pos1[i][0], c2, pos2[j][0]])
    del sites_per_read
    scatter_points = np.array(site_interactions)
    del site_interactions
    np.savetxt("{}_interactions.csv".format(args.input_csv), scatter_points,
               header="chr.1;pos.1;chr.2;pos.2", comments="", delimiter=";", encoding="utf-8", fmt="%s")

    # plot the matrix
    # fig, ax = plt.subplots(nrows=1, ncols=1)
    # ax.set_facecolor("k")
    # plt.scatter(scatter_points[:, 0], scatter_points[:, 1])
    # plt.hist2d(scatter_points[:, 0], scatter_points[:, 1], bins=args.bin_size, cmap=plt.get_cmap("inferno"), vmax=15)
    # plt.colorbar()
    # sns.jointplot(scatter_points[:, 0], scatter_points[:, 1], kind="kde", shade=True)
    #
    # plt.show()
