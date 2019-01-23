#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
import re
from glob import iglob

from Bio import Seq, SeqIO, Restriction
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def subsample_fastq(fastq_files, fasta_out, sample_size=1000, restriction_enzyme="DpnII", fragment_length=-1,
                    fragment_interval=-1, min_length=-1, max_length=-1):
    restr_enzyme = getattr(Restriction, restriction_enzyme, None)

    if fragment_interval <= 0:
        fragment_interval = fragment_length
    if fragment_length > 0:
        logging.info("Reads will be split into %d long fragments with a %d base interval.",
                     fragment_length, fragment_interval)
        max_length = -1
    elif restr_enzyme is None:
        logging.warning("Restriction enzyme %s is not valid; continuing without fragmentation.", restriction_enzyme)
    else:
        logging.info("Using restriction enzyme %s for digestion.", restriction_enzyme)

    for fastq_in in fastq_files:
        # create a unique FASTQ identifier
        fname = os.path.basename(fastq_in).partition(".")[0]
        fname = re.sub("[\W]+", "_", fname)

        logging.info("Sampling %s reads from %s ...", sample_size if sample_size > 0 else "all", fastq_in)
        total_fragments = 0
        with gzip.open(fastq_in, "rt") as fin:
            line = 0
            rec_idx = 0
            for record in fin:
                line += 1
                if (line + 2) % 4 != 0:  # sequence lines in FASTQ
                    continue
                seq = Seq.Seq(record.strip(), alphabet=IUPAC.ambiguous_dna)
                fgmt_idx = 0
                fragments = []
                if fragment_length > 0:
                    # fixed-length fragmentation mode
                    fragments = [seq[i:i + fragment_length] for i in
                                 range(0, len(seq), fragment_interval)]
                elif restr_enzyme is not None:
                    # enzyme-based fragmentation mode
                    fragments = restr_enzyme.catalyse(seq)
                else:
                    # no fragmentation
                    fragments = [seq]

                start_pos = 0
                for digestion in fragments:
                    fr_len = len(digestion)
                    """
                    FASTA header makeup:
                    --------------------
                    Fq.Id: Source FASTQ UID
                    Rd.Id: Source read UID
                    Rd.Ln: Length of the source read
                    Fr.Id: Current fragment number of the read
                    Fr.Ln: Length of the fragment (Fr.Ln == Fr.En - Fr.St)
                    Fr.St: Zero-based start index of the fragment on the read (inclusive)
                    Fr.En: Zero-based end index of the fragment on the read (exclusive)
                    """
                    rec = SeqRecord(
                        digestion,
                        id="Fq.Id:{:s};Rd.Id:{:d};Rd.Ln:{:d};Fr.Id:{:d};Fr.Ln:{:d};Fr.St:{:d};Fr.En:{:d}".format(
                            fname, rec_idx, len(seq), fgmt_idx, fr_len, start_pos, start_pos + fr_len),
                        name="",
                        description=""
                    )
                    start_pos += fr_len

                    if min_length > 0 and len(digestion) < min_length:
                        logging.debug("Skipping %s; too short.", rec.id)
                        continue
                    if 0 < max_length < len(digestion):
                        logging.debug("Skipping %s; too long.", rec.id)
                        continue

                    SeqIO.write(rec, fasta_out, "fasta")
                    fgmt_idx += 1
                    total_fragments += 1
                rec_idx += 1
                if 0 < sample_size <= rec_idx:
                    break
        logging.info("Fragmentation resulted in %d fragments", total_fragments)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s]: %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("output_fa", help="Output location of the FASTA file (gzipped)", metavar="FASTA",
                        action="store", type=str)
    parser.add_argument("input_fq", help="Input FASTQ files (gzipped)", metavar="FASTQ", action="store", type=str,
                        nargs="+")
    parser.add_argument("-n", "--sample", help="Limit the number of reads to extract from each FASTQ file.",
                        metavar="SIZE", action="store", type=int, default=-1)
    parser.add_argument("-e", "--enzyme", help="Restriction enzyme to use for digestion. Leave empty for no digestion.",
                        metavar="ENZYME", action="store", type=str, default="N/a")
    parser.add_argument("-f", "--fragment-length", help="Length of the fragments. This will override --enzyme and "
                                                        "--max-length.",
                        metavar="LENGTH", action="store", type=int, default=-1)
    parser.add_argument("-i", "--fragment-interval", help="Interval between the fragments. This value is only used "
                                                          "when --fragment-length is also specified. Defaults to "
                                                          "--fragment-length",
                        metavar="LENGTH", action="store", type=int, default=-1)
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
        subsample_fastq(input_files, fout, args.sample, args.enzyme, args.fragment_length, args.fragment_interval,
                        args.min, args.max)
