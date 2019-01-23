#!/usr/bin/env python3

import argparse
import logging
import re

import pysam

cigar_decoder = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "=",
    8: "X"
}


def find_read_mappings(read_id, sam_file):
    alignments = {}
    with pysam.AlignmentFile(sam_file, "r") as sam:
        for read in sam.fetch():
            if read.is_unmapped:
                logging.debug("Skipping %s (unaligned)", read.qname)
                continue
            # elif read.is_reverse:
            #     logging.debug("Skipping %s (reverse-complement)", read.qname)
            #     continue

            read_metadata = read_header_to_dict(read.qname)
            rdid = "{}_{}".format(read_metadata["Fq.Id"], read_metadata["Rd.Id"])
            if rdid != read_id:
                continue

            alignments.setdefault(int(read_metadata["Fr.Id"]), list()).append(read)
    return alignments


def read_header_to_dict(header):
    read_metadata = {}
    for x in str(header).split(";"):
        read_metadata[x.split(":")[0]] = x.split(":")[1]
    return read_metadata


def mapping_to_read_coverage(read):
    md_list = ""
    for k, v in read.tags:
        if k == "MD":
            md_list = re.findall(r"(?:[A-Z])|(?:\^[A-Z]+)|(?:[0-9]+)", v)
            break
    md_string = ""
    for x in md_list:
        try:  # match
            matches = int(x)
            md_string += (">" if read.is_reverse else "<") * matches
            continue
        except ValueError:
            pass
        if x[0] == "^":  # insertion
            # md_string += x[1:].upper()
            pass
        else:  # mismatch
            md_string += x.lower()

    cov = ""
    for k, v in read.cigar:
        if k in (4, 5):  # clipped
            cov += " " * v
        elif k in (0, 7, 8):  # aligned
            # cov += cigar_decoder.get(k, "#") * v
            cov += md_string[:v]
            md_string = md_string[v:]
        elif k == 1:  # insertion
            cov += cigar_decoder.get(k, "#") * v
        # elif not md_added:
        #     cov += md_string
    return cov


def get_fragment_origin(read):
    rlen = sum(map(lambda x: 0 if x[0] in (2, 3, 6) else x[1], read.cigar))
    ori_start = 0
    ori_end = rlen

    match_type, match_length = read.cigar[0]
    if match_type in (4, 5):
        ori_start = match_length
    match_type, match_length = read.cigar[-1]
    if match_type in (4, 5):
        ori_end = rlen - match_length

    return ori_start, ori_end


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s]: %(message)s")

    # get command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("read_id", help="Read ID to look for. It is formatted like {Fq.Id}_{Rd.Id}.", metavar="ID",
                        action="store", type=str)
    parser.add_argument("input_sam", help="SAM/BAM file.", metavar="SAM", action="store", type=str)
    args = parser.parse_args()

    aln = find_read_mappings(args.read_id, args.input_sam)

    # print(">>>{:s}".format(args.read_id))

    print("ref;strand;start;end")

    seq = ""
    for fragment_id in sorted(aln.keys()):
        # print(">fragment_{:d} ({:d} aln)".format(fragment_id, len(aln[fragment_id])))

        for x in aln[fragment_id]:
            metadata = read_header_to_dict(x.qname)
            start_idx = int(metadata["Fr.St"])
            if not seq:
                seq = "_" * int(metadata["Rd.Ln"])

            if not x.is_supplementary:
                seq = "".join([seq[:start_idx], x.seq, seq[start_idx + len(x.seq):]])

            coverage = mapping_to_read_coverage(x)
            frag_start, frag_end = get_fragment_origin(x)

            # print(" "*start_idx + coverage)
            print(x.reference_name, 2 if x.is_reverse else 1, start_idx + frag_start, start_idx + frag_end, sep=";")
    print("read;3;0;{:d}".format(len(seq)))

    seq_fmt = re.sub(r"(GATC)", "\033[1;97;41m\\g<0>\033[0m", seq)
    seq_fmt = re.sub(r"(GAT[^C])|(GA[^T]C)|(G[^A]TC)|([^G]ATC)", "\033[0;97;106m\\g<0>\033[0m", seq_fmt)
    # print(seq_fmt)
