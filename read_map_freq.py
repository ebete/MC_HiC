#!/usr/bin/env python3

import argparse
import logging

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


class WeightedLinkedList(object):
    """
    Linked list with weights assigned to the edges between vertices.
    """

    def __init__(self):
        """Initialise a new WeightedLinkedList object."""
        self.head = None
        self.num_nodes = 0

    def insert_head(self, fragment_start, fragment_end, ref_name):
        """
        Insert a new vertex at the head of the linked list.

        :type fragment_start: int
        :param fragment_start: Start position of the fragment.

        :type fragment_end: int
        :param fragment_end: End position of the fragment.

        :type ref_name: str
        :param ref_name: Reference sequence the fragment lies on.

        :rtype WeightedLinkedList
        :return: Itself.
        """
        value = WeightedLinkedListNode(fragment_start, fragment_end, ref_name)
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
        for x in self:
            print(x)

    def __iter__(self):
        cur = self.head
        while cur:
            yield cur
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

    def __init__(self, fragment_start, fragment_end, ref_name):
        """
        Initialises a new WeightedLinkedListNode.

        :type fragment_start: int
        :param fragment_start: Start position of the fragment.

        :type fragment_end: int
        :param fragment_end: End position of the fragment.

        :type ref_name: str
        :param ref_name: Reference sequence the fragment lies on.
        """
        self.fragment_start = fragment_start
        self.fragment_end = fragment_end
        self.ref_name = ref_name
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
        if self.vertex_endpoint is None or self.ref_name != self.vertex_endpoint.ref_name:
            return None
        else:
            return self.fragment_start - self.vertex_endpoint.fragment_end

    def collapse(self):
        """
        Merges this node with the next node in the linked list, creating a new
        fragment that spans both nodes.

        :rtype WeightedLinkedListNode
        :return: Itself.
        """
        if self.vertex_endpoint is None:  # last in linked list
            return self
        if self.ref_name != self.vertex_endpoint.ref_name:  # next is on different chromosome
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
        return "WeightedLinkedListNode{{id: {}; weight: {}; value: ({}:{}-{})}}" \
            .format(self.node_id, self.get_weight(), self.ref_name, self.fragment_start, self.fragment_end)


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
        if read.flag & 0x4 != 0:
            logging.debug("Skipping %s (unaligned)", read.qname)
            continue
        if read.mapq < mapq_cutoff:
            logging.debug("Skipping %s (MAPQ: %d)", read.qname, read.mapq)
            continue
        if read.flag & 0x100 != 0:
            logging.debug("Skipping %s (secondary alignment)", read.qname)
            continue

        # fragment id (Fr.Id) should be ignored
        read_metadata = {}
        for x in str(read.qname).split(";"):
            read_metadata[x.split(":")[0]] = x.split(":")[1]
        rdid = "{}_{}".format(read_metadata["Fq.Id"], read_metadata["Rd.Id"])

        logging.debug("Read %s found at %s:%d-%d", read.qname, read.reference_name, read.reference_start,
                      read.reference_end)
        mapped_reads.setdefault(rdid, WeightedLinkedList()).insert_head(read.reference_start, read.reference_end,
                                                                        read.reference_name)
    samfile.close()

    return mapped_reads


def merge_and_count_freq(reads, dist_cutoff):
    """
    Merge fragments into a single fragment that are mapped close to each other.
    The frequency of mapped fragments per read are returned as a Numpy array.

    :type reads: dict
    :param reads: Dictionary containing a WeightedLinkedList associated with the
                  read.

    :type dist_cutoff: int
    :param dist_cutoff: Minimum distance between fragments before they are no
                        longer merged.

    :rtype: numpy.core.multiarray
    :return: Array containing the frequencies of the number of fragments per
             read mapped.
    """
    logging.info("Merging fragments fewer than %d bases removed from each other ...", dist_cutoff)
    counts = []
    for key, value in reads.items():
        counts.append(len(value.collapse_nodes(dist_cutoff)))
    return np.array(np.unique(counts, return_counts=True))


def export_mapped_regions(reads, outfile):
    """
    Write the mapped regions to a CSV file.

    :type reads: dict
    :param reads: Dictionary containing a WeightedLinkedList associated with the
                  read.

    :type outfile: str
    :param outfile: Where to write the CSV file to.
    :return:
    """
    with open(outfile, "wt") as handle:
        logging.info("Writing mapped regions to %s ...", outfile)
        handle.write("read;chromosome;start;end\n")
        for key, value in reads.items():
            fragment_string = ""
            for fragment in value:
                fragment_string = "{};{};{};{}\n".format(key, fragment.ref_name, fragment.fragment_start,
                                                         fragment.fragment_end) + fragment_string
            handle.write(fragment_string)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="INFILE", action="store", type=str)
    parser.add_argument("-c", "--csv-output", help="Output location of the CSV mapped regions file.", metavar="CSV",
                        action="store", type=str)
    parser.add_argument("-d", "--distance-cutoff", help="Minimum distance between two fragments before considering "
                                                        "them as separate.",
                        metavar="CUTOFF", action="store", type=int, default=1000)
    parser.add_argument("-q", "--minimum-mapq", help="Minimum MAPQ that a fragment needs to have.", metavar="MQ",
                        action="store", type=int, default=1)
    parser.add_argument("-r", "--region", help="Limit read to specific region.", metavar="REGION",
                        action="store", type=str, default=".")
    args = parser.parse_args()

    mapped_reads = read_sam(args.input_sam, args.minimum_mapq, args.region)
    arr = merge_and_count_freq(mapped_reads, args.distance_cutoff)
    if args.csv_output is not None:
        export_mapped_regions(mapped_reads, args.csv_output)

    logging.shutdown()
