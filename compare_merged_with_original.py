#!/usr/bin/env python3

import argparse
import logging

import utils


def compare_merged_to_original(original_sam, mergemap_sam):
    """
    Print statistics of alignments that occur in both SAM files to a CSV format.

    :param original_sam: SAM file used as input for the MergeMap SAM file.
    :param mergemap_sam: SAM file generated by MergeMap.
    """
    min_mapq_increase = 5
    min_appended_length = 10

    original_metadata_dict = utils.get_mapping_metadata(original_sam)
    mergemap_metadata_dict = utils.get_mapping_metadata(mergemap_sam)

    logging.info("Comparing original alignments to mergemapped alignments ...")
    print("original_qname",  # ID of the original fragment
          "original_length",  # unclipped length of the original fragment
          "original_mapq",  # MAPQ of the original fragment
          "mergemap_qname",  # ID of the MergeMap fragment
          "mergemap_length",  # unclipped length of the MergeMap fragment
          "mergemap_mapq",  # MAPQ of the MergeMap fragment
          "effective_appended",  # effective length added (appended fragment - clipping)
          sep="\t")

    for read_id, mergemap_fragments in mergemap_metadata_dict.items():
        original_metadata = original_metadata_dict.get(read_id, dict())
        if not original_metadata:
            logging.warning("Read %s not mapped in original dataset. Is this the correct SAM file pair?", read_id)
            continue

        for mergemap_frid, metadata in mergemap_fragments.items():
            original_fragments = [int(x) for x in metadata["Src.Fr"].split(",")]
            ori_fragment = original_fragments[0 if metadata["Src.Ori"] == "Left" else -1]
            if ori_fragment not in original_metadata:
                logging.warning("Fragment %d of read %s not mapped in original dataset. "
                                "Is this the correct SAM file pair?", ori_fragment, read_id)
                continue

            # calculate effective appended length (part of appended fragment that is actually mapped)
            appended_maplen = int(metadata["Src.Ln"]) - metadata["clipping"][
                "clip_start" if metadata["Src.Ori"] == "Right" else "clip_end"]
            mapq_increase = metadata["mapq"] - original_metadata[ori_fragment]["mapq"]
            if appended_maplen < 0 or (mapq_increase < min_mapq_increase and appended_maplen < min_appended_length):
                logging.debug("Skipping fragment %d of read %s (MAPQ increase/appended length insufficient).",
                              ori_fragment, read_id)
                continue

            print(original_metadata[ori_fragment]["qname"],
                  original_metadata[ori_fragment]["clipping"]["unclipped_length"],
                  original_metadata[ori_fragment]["mapq"],
                  metadata["qname"],
                  metadata["clipping"]["unclipped_length"],
                  metadata["mapq"],
                  appended_maplen,
                  sep="\t")


if __name__ == "__main__":
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("original_sam", help="Original SAM/BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("mergemap_sam", help="SAM/BAM file created by MergeMap.", metavar="SAM", action="store",
                        type=str)
    args = parser.parse_args()

    compare_merged_to_original(args.original_sam, args.mergemap_sam)

    logging.shutdown()
