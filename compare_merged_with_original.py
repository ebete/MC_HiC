#!/usr/bin/env python3

import argparse
import logging

import pysam


def compare_merged_to_original(original_sam, mergemap_sam):
    original_metadata_dict = get_mapping_metadata(original_sam)
    mergemap_metadata_dict = get_mapping_metadata(mergemap_sam)

    logging.info("Comparing original alignments to mergemapped alignments ...")
    print("original_qname", "original_length", "mergemap_qname", "mergemap_length", sep="\t")
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

            print(original_metadata[ori_fragment]["qname"], original_metadata[ori_fragment]["unclipped_len"],
                  metadata["qname"], metadata["unclipped_len"], sep="\t")
        break  # TODO: remove


def get_mapping_metadata(sam_input):
    alignments = dict()

    logging.info("Reading SAM file %s ...", sam_input)
    with pysam.AlignmentFile(sam_input, "r") as sam:
        for read in sam:
            if read.is_unmapped:
                continue

            metadata = read_header_to_dict(read.qname)
            metadata["unclipped_len"] = unclipped_length(read.cigar)
            metadata["qname"] = read.qname

            rdid = "{}_{}".format(metadata["Fq.Id"], metadata["Rd.Id"])

            alignments.setdefault(rdid, dict())[int(metadata["Fr.Id"])] = metadata

    return alignments


def read_header_to_dict(header):
    """
    Map the metadata stored in the FASTA header to a dictionary.

    :type header: str
    :param header: The FASTA header as generated by fastq_subsample_and_merge.py

    :rtype dict
    :return: A dictionary representing the key-value pairs in the header
    """
    read_metadata = dict()
    for x in str(header).split(";"):
        read_metadata[x.split(":")[0]] = x.split(":")[1]
    return read_metadata


def unclipped_length(cigar_tuple):
    unclipped = 0
    for operation, length in cigar_tuple:
        if operation in {4, 5}:  # skip clipped
            continue
        unclipped += length
    return unclipped


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("original_sam", help="Original SAM/BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("mergemap_sam", help="SAM/BAM file created by MergeMap.", metavar="SAM", action="store",
                        type=str)
    args = parser.parse_args()

    compare_merged_to_original(args.original_sam, args.mergemap_sam)

    logging.shutdown()
