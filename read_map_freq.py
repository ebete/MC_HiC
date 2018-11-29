#!/usr/bin/env python3

import argparse
import glob
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

plt.style.use("ggplot")
sns.set()


class WeightedLinkedList(object):
    """
    Linked list with weights assigned to the edges between vertices.
    """

    def __init__(self):
        """Initialise a new WeightedLinkedList object."""
        self.head = None
        self.num_nodes = 0

    def insert_head(self, fragment_start, fragment_end):
        """
        Insert a new vertex at the head of the linked list.

        :type fragment_start: int
        :param fragment_start: Start position of the fragment.

        :type fragment_end: int
        :param fragment_end: End position of the fragment.

        :rtype WeightedLinkedList
        :return: Itself.
        """
        value = WeightedLinkedListNode(fragment_start, fragment_end)
        value.set_endpoint(self.head)
        self.head = value
        self.num_nodes += 1
        return self

    def collapse_nodes(self, min_weight=1000):
        """
        Collapses all vertices that have edges with a weight less than
        min_weight.

        :type min_weight: int
        :param min_weight: Minimum weight of an edge to not collapse the vertex.

        :rtype WeightedLinkedList
        :return: Itself.
        """
        cur = self.head
        while cur:
            if cur.get_weight() is None or cur.get_weight() >= min_weight:
                cur = cur.vertex_endpoint
                continue
            cur.collapse()
            self.num_nodes -= 1
        return self

    def print_list(self):
        """
        Prints all vertices contained in the WeightedLinkedList in order.
        """
        cur = self.head
        while cur:
            print(cur)
            cur = cur.vertex_endpoint

    def __len__(self):
        return self.num_nodes

    def __str__(self):
        return "WeightedLinkedList{{items: {}}}".format(len(self))


class WeightedLinkedListNode(object):
    """
    A node in a WeightedLinkedList.
    """
    node_id = 0

    def __init__(self, fragment_start, fragment_end):
        """
        Initialises a new WeightedLinkedListNode.

        :type fragment_start: int
        :param fragment_start: Start position of the fragment.

        :type fragment_end: int
        :param fragment_end: End position of the fragment.
        """
        self.fragment_start = fragment_start
        self.fragment_end = fragment_end
        self.vertex_endpoint = None

        self.node_id = str(WeightedLinkedListNode.node_id)
        WeightedLinkedListNode.node_id += 1

    def set_endpoint(self, node):
        """
        Assign a new endpoint for this node.

        :type node: WeightedLinkedListNode
        :param node: The node to point to.

        :rtype WeightedLinkedListNode
        :return: Itself.
        """
        self.vertex_endpoint = node
        logging.debug("New endpoint set (w:{}; end:{})".format(self.get_weight(), node))
        return self

    def get_weight(self):
        """
        Gets the distance between this node and the node that is pointed to.

        :rtype int
        :return: The distance between the nodes.
        """
        return (self.fragment_start - self.vertex_endpoint.fragment_end) if (self.vertex_endpoint is not None) else None

    def collapse(self):
        """
        Merges this node with the next node in the linked list, creating a new
        fragment that spans both nodes.

        :rtype WeightedLinkedListNode
        :return: Itself.
        """
        if self.vertex_endpoint is None:
            return self
        logging.debug("Collapsing nodes {} and {}".format(self.node_id, self.vertex_endpoint.node_id))

        if self.fragment_end < self.vertex_endpoint.fragment_end:
            self.fragment_end = self.vertex_endpoint.fragment_end

        if self.fragment_start > self.vertex_endpoint.fragment_start:
            self.fragment_start = self.vertex_endpoint.fragment_start

        self.node_id += "_{}".format(self.vertex_endpoint.node_id)
        self.vertex_endpoint = self.vertex_endpoint.vertex_endpoint

        return self

    def __str__(self):
        return "WeightedLinkedListNode{{id: {}; weight: {}; value: ({}:{})}}" \
            .format(self.node_id, self.get_weight(), self.fragment_start, self.fragment_end)


def read_sam(fname, mapq_cutoff=0, sam_region="."):
    """
    Parse reads from a SAM/BAM file and Create a WeightedLinkedList per read.

    :type fname: str
    :param fname: The SAM/BAM file location.

    :type mapq_cutoff: int
    :param mapq_cutoff: Minimum mapping quality (MAPQ/MQ) before considering a mapped fragment.

    :type sam_region: str
    :param sam_region: Region in SAM format to limit the reading to.

    :rtype dict
    :return: Dictionary containing a WeightedLinkedList associated with the read.
    """
    mapped_reads = {}

    samfile = pysam.AlignmentFile(fname, "r")
    logging.info("Reading alignment file %s ...", fname)
    for read in samfile.fetch(region=sam_region):
        if read.mapq < mapq_cutoff:
            logging.debug("Skipping %s (MAPQ: %d)", read.qname, read.mapq)
            continue

        logging.debug("Read %s found at %d:%d", read.qname, read.reference_start, read.reference_end)
        mapped_reads.setdefault(read.qname, WeightedLinkedList()).insert_head(read.reference_start, read.reference_end)
    samfile.close()

    return mapped_reads


def merge_and_count_freq(reads, dist_cutoff):
    """
    Merge fragments into a single fragment that are mapped close to each other.
    The frequency of mapped fragments per read are returned as a DataFrame.

    :type reads: dict
    :param reads: Dictionary containing a WeightedLinkedList associated with the
                  read.

    :type dist_cutoff: int
    :param dist_cutoff: Minimum distance between fragments before they are no
                        longer merged.

    :rtype: pd.DataFrame
    :return: DataFrame containing the frequencies of the number of fragments per
             read mapped.
    """
    logging.info("Merging fragments fewer than %d bases removed from each other ...", dist_cutoff)
    counts = []
    for key, value in reads.items():
        counts.append(len(value.collapse_nodes(dist_cutoff)))
    bins = np.array(np.unique(counts, return_counts=True))
    return pd.DataFrame(bins[1:, :], columns=bins[0])


def plot_frequencies(frequencies, outfile):
    """
    Create a matplotlib bar plot from a pandas DataFrame.

    :type frequencies: pd.DataFrame
    :param frequencies: DataFrame containing the frequencies of the number of
                        fragments per read mapped.

    :type outfile: str
    :param outfile: Location of the PDF output file.
    """
    logging.info("Generating bar plot and writing to %s", outfile)
    with PdfPages(outfile) as pdf:
        # log scale
        frequencies.plot(kind="bar")
        plt.title("Mapped fragments per read")
        plt.xlabel("Mapped fragments")
        plt.ylabel("Frequency ($\\log_{10}$)")
        plt.semilogy()
        pdf.savefig(bbox_inches="tight")
        plt.close()
        # linear scale
        frequencies.plot(kind="bar")
        plt.title("Mapped fragments per read")
        plt.xlabel("Mapped fragments")
        plt.ylabel("Frequency")
        pdf.savefig(bbox_inches="tight")
        plt.close()
        # annotated heatmap
        df = (frequencies / frequencies.sum()).transpose()
        sns.heatmap(df * 100, annot=True, fmt=".0f", vmin=0, vmax=100, square=True, cmap="inferno")
        pdf.savefig(bbox_inches="tight")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM files.", metavar="INFILE", action="store", type=str,
                        nargs="+")
    parser.add_argument("--img-output", "-o", help="Output location of the PDF images", metavar="PDF", action="store",
                        type=str)
    parser.add_argument("--distance-cutoff", "-d", help="Minimum distance between two fragments before considering "
                                                        "them as separate",
                        metavar="CUTOFF", action="store", type=int, default=1000)
    parser.add_argument("--minimum-mapq", "-q", help="Minimum MAPQ that a fragment needs to have", metavar="MQ",
                        action="store", type=int, default=1)
    parser.add_argument("--region", "-r", help="Limit read to specific region", metavar="REGION",
                        action="store", type=str, default=".")
    args = parser.parse_args()

    bins = pd.DataFrame()
    xlab = []
    for globfile in args.input_sam:
        for samfile in glob.iglob(globfile):
            mapped_reads = read_sam(samfile, args.minimum_mapq, args.region)
            df = merge_and_count_freq(mapped_reads, args.distance_cutoff)
            bins = pd.concat([bins, df], axis=0, ignore_index=True)
            xlab.append(os.path.basename(samfile))

    bins.fillna(0, inplace=True)
    bins = bins.transpose()
    bins.columns = xlab

    plot_frequencies(bins, args.img_output)

    logging.shutdown()
