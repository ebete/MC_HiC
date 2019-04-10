#!/usr/bin/env python3

import argparse
import gzip
import logging

from Bio import SeqIO, Restriction, SeqRecord

import utils


def get_cut_region(reference_genome, fwd_primer, rev_primer, restriction_enzyme):
    """
    Print the viewpoint region in FASTA format based on the primer positions.
    Will extend up to 2 restriction sites away.

    :param reference_genome: Genome where to extract the viewpoint region from
    :param fwd_primer: Position of the viewpoint forward primer.
    :param rev_primer: Position of the viewpoint reverse primer.
    :param restriction_enzyme: Restriction enzyme used in digestion.
    """
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
    """
    Convert a locus string (chr:start-end) to a dictionary.

    :param position_string: Genomic locus string.
    :return: A dictionary representation of the locus.
    """
    chromosome = position_string.split(":")
    locus = sorted([int(pos) for pos in chromosome[1].split("-")])
    if len(locus) != 2:
        raise ValueError("Locus string is malformed.")
    return {"chr": chromosome[0], "start": locus[0], "end": locus[1]}


if __name__ == "__main__":
    utils.init_logger()

    # get command-line arguments
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
