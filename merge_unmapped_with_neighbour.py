#!/usr/bin/env python3

import argparse
import gzip
import logging

import pysam
from Bio import Seq, SeqIO, SeqRecord
from Bio.Alphabet import IUPAC


def get_mapped_fragments(sam_input):
    alignments = dict()

    logging.info("Reading SAM file %s ...", sam_input)
    with pysam.AlignmentFile(sam_input, "r") as sam:
        for read in sam:
            if read.is_unmapped:
                continue

            metadata = read_header_to_dict(read.qname)
            rdid = "{}_{}".format(metadata["Fq.Id"], metadata["Rd.Id"])

            alignments.setdefault(rdid, list()).append(int(metadata["Fr.Id"]))

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


def combine_unmapped_and_mapped(fasta_file, mapped_fragments):
    logging.info("Parsing fragments from %s and comparing to mapped ...", fasta_file)
    total_reads = 0
    unmapped_reads = 0
    with gzip.open(fasta_file, "rt") as handle:
        read_records = list()
        last_id = ""
        for record in SeqIO.parse(handle, "fasta", IUPAC.unambiguous_dna):
            record_metadata = read_header_to_dict(record.id)
            fasta_id = "{}_{}".format(record_metadata["Fq.Id"], record_metadata["Rd.Id"])

            if fasta_id == last_id:
                read_records.append(record)
                continue

            total_reads += 1
            if last_id in mapped_fragments:
                mapped_fragments[last_id] = sorted(set(mapped_fragments[last_id]))
                fragment_id = 0
                for new_record in do_merge(read_records, mapped_fragments[last_id]):
                    new_record.id = "Fq.Id:{:s};Rd.Id:{:s};Rd.Ln:{:s};Fr.Id:{:d};Fr.Ln:{:d};{:s}".format(
                        record_metadata["Fq.Id"], record_metadata["Rd.Id"], record_metadata["Rd.Ln"], fragment_id,
                        len(new_record.seq), new_record.id)
                    print(new_record.format("fasta"), end="")
                    fragment_id += 1
            else:
                unmapped_reads += 1

            last_id = fasta_id
            read_records = [record]
    logging.info("Out of %d reads, %d had no mapped fragments (%.1f%%).", total_reads, unmapped_reads,
                 unmapped_reads / total_reads * 100)


def do_merge(fasta_records, mapped_fragments):
    new_seqs = list()

    mapped_count = len(mapped_fragments)
    if mapped_count < 2:  # no in-between unmapped fragments exist
        return new_seqs

    last_id = mapped_fragments[0]
    for fragment, idx in zip(mapped_fragments, range(mapped_count)):
        jump_size = fragment - last_id - 1  # number of fragments skipped
        last_id = fragment

        if jump_size <= 0 or idx < 1:  # skip contiguous and first mapped fragment
            continue

        jump_region = list(range(mapped_fragments[idx - 1], fragment + 1))
        for i in range(len(jump_region) - 2):
            new_seqs.append(SeqRecord.SeqRecord(  # Extend from left
                Seq.Seq("".join([str(fasta_records[x].seq) for x in jump_region[:i + 2]]),
                        alphabet=IUPAC.unambiguous_dna),
                id="Src.Fr:{:s};Src.Ori:Left".format(",".join([str(x) for x in jump_region[:i + 2]])),
                description="", name=""
            ))
            new_seqs.append(SeqRecord.SeqRecord(  # Extend from right
                Seq.Seq("".join([str(fasta_records[x].seq) for x in jump_region[i + 1:]]),
                        alphabet=IUPAC.unambiguous_dna),
                id="Src.Fr:{:s};Src.Ori:Right".format(",".join([str(x) for x in jump_region[i + 1:]])),
                description="", name=""
            ))
    return new_seqs


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("input_fasta", help="Input FASTA file.", metavar="FASTA", action="store", type=str)
    args = parser.parse_args()

    mapped_fragments = get_mapped_fragments(args.input_sam)
    combine_unmapped_and_mapped(args.input_fasta, mapped_fragments)

    logging.shutdown()
