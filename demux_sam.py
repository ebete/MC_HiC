#!/usr/bin/env python3

import argparse
import glob
import logging
import os

import pysam


def demux_sam(sam_input, map_counts, outdir="./", sam_region="."):
    input_name = os.path.basename(sam_input).partition(".")[0]
    demuxed_reads = {}
    with pysam.AlignmentFile(sam_input, "r") as samfile:
        logging.info("Demuxing alignment file %s ...", sam_input)
        for read in samfile.fetch(region=sam_region):
            read_metadata = {}
            for x in str(read.qname).split(";"):
                read_metadata[x.split(":")[0]] = x.split(":")[1]

            outfile = os.path.join(outdir, "{}_{}.bam".format(input_name, read_metadata.get("Fq.Id", "Unknown")))
            demuxed_reads.setdefault(outfile, []).append(read)
            map_counts.setdefault(read_metadata.get("Fq.Id", "Unknown"), []).append(read_metadata.get("Rd.Id", "-1"))

        for outfile, reads in demuxed_reads.items():
            logging.info("Writing %d reads to %s ...", len(reads), outfile)
            with pysam.AlignmentFile(outfile, "wb", template=samfile) as handle:
                for read in reads:
                    handle.write(read)
            logging.info("Creating index of %s ...", outfile)
            pysam.index(outfile)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM files.", metavar="INFILE", action="store", type=str,
                        nargs="+")
    parser.add_argument("--output", "-o", help="Output directory of the de-multiplexed alignment files", metavar="DIR",
                        action="store", type=str, default="./")
    parser.add_argument("--region", "-r", help="Limit read exporting to specific region", metavar="REGION",
                        action="store", type=str, default=".")
    args = parser.parse_args()

    count_per_file = {}
    for globfile in args.input_sam:
        for samfile in glob.iglob(globfile):
            demux_sam(samfile, count_per_file, args.output, sam_region=args.region)

    for fq_id, rd_ids in count_per_file.items():
        counts = len(set(rd_ids))
        logging.info("Unique reads in %s: %d", fq_id, counts)