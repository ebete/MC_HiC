#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
from glob import iglob

from Bio import SeqIO, Restriction
from Bio.SeqRecord import SeqRecord


def subsample_fastq(fastq_files, fasta_out, sample_size=1000, restriction_enzyme="DpnII", min_length=-1, max_length=-1):
    restr_enzyme = getattr(Restriction, restriction_enzyme, None)
    if restr_enzyme is None:
        logging.warning("Restriction enzyme %s not found; continuing without digesting.", restriction_enzyme)

    for fastq_in in fastq_files:
        fname = os.path.basename(fastq_in).split("_")[1]

        logging.info("Sampling %s reads from %s ...", sample_size if sample_size > 0 else "all", fastq_in)
        total_fragments = 0
        with gzip.open(fastq_in, "rt") as fin:
            handle = SeqIO.parse(fin, "fastq")
            rec_idx = 0
            for record in handle:
                fgmt_idx = 0
                for digestion in (restr_enzyme.catalyse(record.seq) if restr_enzyme is not None else [record.seq]):
                    rec = SeqRecord(
                        digestion,
                        id="Fq.Id:{:s};Rd.Id:{:d};Rd.Ln:{:d};Fr.Id:{:d};Fr.Ln:{:d}".format(
                            fname, rec_idx, len(record.seq), fgmt_idx, len(digestion)),
                        name="",
                        description=""
                    )
                    if min_length >= 0 and len(digestion) < min_length:
                        logging.debug("Skipping %s; too short.", rec.id)
                        continue
                    if 0 <= max_length < len(digestion):
                        logging.debug("Skipping %s; too long.", rec.id)
                        continue

                    SeqIO.write(rec, fasta_out, "fasta")
                    fgmt_idx += 1
                    total_fragments += 1
                rec_idx += 1
                if 0 < sample_size <= rec_idx:
                    break
        logging.info("Digestion with %s resulted in %d fragments", restr_enzyme, total_fragments)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s]: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("output_fa", help="Output location of the FASTA file (gzipped)", metavar="FASTA",
                        action="store", type=str)
    parser.add_argument("input_fq", help="Input FASTQ files (gzipped)", metavar="FASTQ", action="store", type=str,
                        nargs="+")
    parser.add_argument("-n", "--sample", help="Number of reads to extract from each FASTQ file. Use a negative value "
                                               "to extract all reads.", metavar="SIZE", action="store", type=int,
                        default=1000)
    parser.add_argument("-e", "--enzyme", help="Restriction enzyme to use for digestion. Leave empty for no digestion.",
                        metavar="ENZYME", action="store", type=str, default="N/a")
    parser.add_argument("-s", "--min", help="Minimum length of a read fragment.", metavar="LENGTH", action="store",
                        type=int, default=-1)
    parser.add_argument("-l", "--max", help="Maximum length of a read fragment.", metavar="LENGTH", action="store",
                        type=int, default=-1)
    args = parser.parse_args()

    input_files = []
    for glob in args.input_fq:
        for fastq in iglob(glob):
            input_files.append(fastq)

    with gzip.open(args.output_fa, "wt") as fout:
        subsample_fastq(input_files, fout, args.sample, args.enzyme, args.min, args.max)
