#!/usr/bin/env python3

import gzip
import os
import sys

from Bio import SeqIO


def subsample_fastq(indir, fastq_out):
    for fastq_in in os.listdir(indir):
        fastq_in = os.path.join(indir, fastq_in)
        fname = os.path.basename(fastq_in).split("_")[1]
        with gzip.open(fastq_in, "rt") as fin:
            handle = SeqIO.parse(fin, "fastq")
            rec_idx = 0
            for record in handle:
                record.id = "Fq.Id:{:s};Rd.Id:{:d};Rd.Ln:{:d}".format(fname, rec_idx, len(record.seq))
                record.description = ""
                SeqIO.write(record, fastq_out, "fastq")

                rec_idx += 1
                if rec_idx >= 1000:
                    break


if __name__ == '__main__':
    with gzip.open(sys.argv[2], "wt") as fout:
        subsample_fastq(sys.argv[1], fout)
