#!/usr/bin/env python3

import argparse
import gzip
import logging

from Bio import SeqIO, Restriction, SeqRecord


def get_cut_region(reference_genome, fwd_primer, rev_primer, restriction_enzyme):
    restr_enzyme = getattr(Restriction, restriction_enzyme)
    if fwd_primer["chr"] != rev_primer["chr"]:
        raise ValueError("The primers are not located on the same chromosome.")

    logging.info("Reading reference %s ...", reference_genome)
    with gzip.open(reference_genome, "rt") as handle:
        for chromosome in SeqIO.parse(handle, "fasta"):
            if chromosome.id != fwd_primer["chr"]:
                continue
            logging.info("Scanning %s for cut sites ...", chromosome.id)

            fwd_restriction_pos = str(chromosome.seq).find(restr_enzyme.site, fwd_primer["end"])
            fwd_restriction_pos = str(chromosome.seq).find(restr_enzyme.site, fwd_restriction_pos + 1)  # 2nd site
            fwd_restriction_pos += len(restr_enzyme)  # include restriction site

            rev_restriction_pos = str(chromosome.seq).rfind(restr_enzyme.site, 0, rev_primer["start"])
            rev_restriction_pos = str(chromosome.seq).rfind(restr_enzyme.site, 0, rev_restriction_pos)  # 2nd site

            soi = SeqRecord.SeqRecord(seq=chromosome.seq[rev_restriction_pos:fwd_restriction_pos],
                                      id="{:s}:{:d}-{:d}".format(
                                          chromosome.id, rev_restriction_pos, fwd_restriction_pos),
                                      description="", name="")
            print(soi.format("fasta"), end="")
            return
    logging.error("Chromosome %s does not exist in reference.", reference_genome)


def make_locus(position_string):
    chromosome = position_string.split(":")
    locus = sorted([int(pos) for pos in chromosome[1].split("-")])
    if len(locus) != 2:
        raise ValueError("Locus string is malformed.")
    return {"chr": chromosome[0], "start": locus[0], "end": locus[1]}


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")

    # Get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_genome", help="Input FASTA genome (gzipped).", metavar="FASTA", action="store", type=str)
    parser.add_argument("forward_primer", help="Start and end of the forward primer. Format like chr1:50-100.",
                        metavar="FWD", action="store", type=str)
    parser.add_argument("reverse_primer", help="Start and end of the reverse primer. Format like chr1:50-100.",
                        metavar="REV", action="store", type=str)
    parser.add_argument("-e", "--enzyme", help="Restriction enzyme used for digestion.", metavar="ENZYME",
                        action="store", type=str, default="DpnII")
    args = parser.parse_args()

    fwd_locus = make_locus(args.forward_primer)
    rev_locus = make_locus(args.reverse_primer)

    get_cut_region(args.input_genome, fwd_locus, rev_locus, args.enzyme)

    logging.shutdown()
