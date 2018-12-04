#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
from glob import iglob

from Bio import SeqIO, Restriction
from Bio.SeqRecord import SeqRecord


def subsample_fastq(fastq_files, fasta_out, sample_size=1000, restriction_enzyme="DpnII"):
    restr_enzyme = getattr(Restriction, restriction_enzyme)

    for fastq_in in fastq_files:
        fname = os.path.basename(fastq_in).split("_")[1]

        logging.info("Sampling %d reads from %s ...", sample_size, fastq_in)
        total_fragments = 0
        with gzip.open(fastq_in, "rt") as fin:
            handle = SeqIO.parse(fin, "fastq")
            rec_idx = 0
            for record in handle:
                fgmt_idx = 0
                for digestion in restr_enzyme.catalyse(record.seq):
                    rec = SeqRecord(
                        digestion,
                        id="Fq.Id:{:s};Rd.Id:{:d};Fr.Id:{:d};Rd.Ln:{:d}".format(
                            fname, rec_idx, fgmt_idx, len(digestion)),
                        name="",
                        description=""
                    )

                    SeqIO.write(rec, fasta_out, "fasta")
                    fgmt_idx += 1
                    total_fragments += 1
                rec_idx += 1
                if rec_idx >= sample_size:
                    break
        logging.info("Digestion with %s resulted in %d fragments", restriction_enzyme, total_fragments)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("output_fa", help="Output location of the FASTA file (gzipped)", metavar="FASTA",
                        action="store", type=str)
    parser.add_argument("input_fq", help="Input FASTQ files (gzipped)", metavar="FASTQ", action="store", type=str,
                        nargs="+")
    parser.add_argument("--sample", "-n", help="Number of reads to extract from each FASTQ file", metavar="N",
                        action="store", type=int, default=1000)
    parser.add_argument("--enzyme", help="Restriction enzyme to use for digestion", metavar="ENZYME",
                        action="store", type=str, default="DpnII")
    args = parser.parse_args()

    input_files = []
    for glob in args.input_fq:
        for fastq in iglob(glob):
            input_files.append(fastq)

    with gzip.open(args.output_fa, "wt") as fout:
        subsample_fastq(input_files, fout, args.sample, args.enzyme)
