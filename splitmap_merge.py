#!/usr/bin/env python3

import argparse
import logging

import utils


def compare_original_to_splitmap(original_sam, splitmap_sam):
    original_metadata_dict = utils.get_mapping_metadata(original_sam)
    splitmap_metadata_dict = utils.get_mapping_metadata(splitmap_sam)

    fmt_set = lambda x: " ".join(map(str, x))

    logging.info("Comparing original alignments to splitmapped alignments ...")
    # print("chr", "start", "end", "read_id", "src", sep="\t")
    # print("intersect", "first_run", "second_run", sep="\t")
    print("read_id", "fragment_id", "event", sep="\t")
    for read_id, orig_fragments in original_metadata_dict.items():
        first_fragments = set(orig_fragments.keys())
        if read_id not in splitmap_metadata_dict:
            for fragment_id, orig_metadata in orig_fragments.items():
                # print(orig_metadata["position"]["chr"], orig_metadata["position"]["start"],
                #       orig_metadata["position"]["end"], read_id, "original", sep="\t")
                # print(fmt_set(first_fragments), "", "", sep="\t")
                print(read_id, fragment_id, "set1", sep="\t")
            continue

        splitmap_fragments = splitmap_metadata_dict[read_id]
        second_fragments = set()
        fragment_lookup = dict()
        for frid, metadata in splitmap_fragments.items():
            srcfr = set(map(int, metadata["Src.Fr"].split(",")))
            for fragment in srcfr:
                fragment_lookup.setdefault(fragment, list()).append(frid)
            second_fragments |= srcfr

        first_only = first_fragments - second_fragments
        second_only = second_fragments - first_fragments
        intersecting = first_fragments & second_fragments

        for frid in intersecting:
            if len(fragment_lookup[frid]) > 1:
                print(read_id, frid, "multimap", sep="\t")
            elif equal_positions(orig_fragments[frid]["position"],
                                 splitmap_fragments[fragment_lookup[frid][0]]["position"]):
                print(read_id, frid, "equal", sep="\t")
            else:
                print(read_id, frid, "different", sep="\t")
        for frid in first_only:
            print(read_id, frid, "set1", sep="\t")
        for frid in second_only:
            print(read_id, frid, "set2", sep="\t")


def equal_positions(first_pos, second_pos, max_dist=0):
    """
    Returns whether the two loci are located within max_dist bases from each other.

    :type first_pos: dict
    :param first_pos: Dictionary containing the first locus.

    :type second_pos: dict
    :param second_pos: Dictionary containing the first locus.

    :type max_dist: int
    :param max_dist: Maximum distance the two loci may be removed from each other.

    :rtype bool
    :return: True if the loci are close to each other, false otherwise.
    """
    if first_pos["chr"] != second_pos["chr"]:
        return False  # different chromosome
    if second_pos["start"] - max_dist <= first_pos["start"] <= second_pos["end"] + max_dist:
        return True  # start position falls within the other fragment
    if second_pos["start"] - max_dist <= first_pos["end"] <= second_pos["end"] + max_dist:
        return True  # end position falls within the other fragment
    return False


if __name__ == "__main__":
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("original_sam", help="Original SAM/BAM file.", metavar="ORIGINAL", action="store", type=str)
    parser.add_argument("splitmap_sam", help="SAM/BAM file created by SplitMap.", metavar="SPLITMAP", action="store",
                        type=str)
    args = parser.parse_args()

    compare_original_to_splitmap(args.original_sam, args.splitmap_sam)

    logging.shutdown()
