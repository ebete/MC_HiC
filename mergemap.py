#!/usr/bin/env python3

import argparse
import gzip
import logging

import pysam
from Bio import SeqIO, SeqRecord
from Bio.Alphabet import IUPAC

import utils


def get_mapped_fragments(sam_input):
    alignments = dict()

    logging.info("Reading SAM file %s ...", sam_input)
    with pysam.AlignmentFile(sam_input, "r") as sam:
        for read in sam:
            if read.is_unmapped:
                continue

            metadata = utils.read_header_to_dict(read.qname)
            rdid = utils.make_read_id(metadata)

            alignments.setdefault(rdid, list()).append(int(metadata["Fr.Id"]))

    return alignments


def combine_unmapped_and_mapped(fasta_file, mapped_fragments):
    logging.info("Parsing fragments from %s and comparing to mapped ...", fasta_file)
    total_reads = 0
    unmapped_reads = -1  # initialisation will set it to zero
    with gzip.open(fasta_file, "rt") as handle:
        read_records = list()
        last_id = ""
        last_metadata = dict()
        for record in SeqIO.parse(handle, "fasta", IUPAC.unambiguous_dna):
            record_metadata = utils.read_header_to_dict(record.id)
            fasta_id = utils.make_read_id(record_metadata)

            if fasta_id == last_id:
                read_records.append(record)
                continue

            total_reads += 1
            if last_id not in mapped_fragments:
                unmapped_reads += 1
            mapped_fragments[last_id] = sorted(set(mapped_fragments.get(last_id, list())))

            fragment_id = 0
            for new_record in merge_unmapped(read_records, mapped_fragments[last_id]):
                new_record.id = "Fq.Id:{:s};Rd.Id:{:s};Rd.Ln:{:s};Fr.Id:{:d};Fr.Ln:{:d};{:s}".format(
                    last_metadata["Fq.Id"], last_metadata["Rd.Id"], last_metadata["Rd.Ln"], fragment_id,
                    len(new_record.seq), new_record.id)
                print(new_record.format("fasta"), end="")
                fragment_id += 1

            last_id = fasta_id
            last_metadata = record_metadata
            read_records = [record]
    logging.info("Out of %d reads, %d had no mapped fragments (%.1f%%).", total_reads, unmapped_reads,
                 unmapped_reads / total_reads * 100)


def merge_unmapped(fasta_records, mapped_fragments):
    """
    Generate new fragments from two unmapped ones. All adjoined unmapped
    fragments will be concatenated to create a single, larger fragment.

    :type fasta_records: list
    :param fasta_records: A list containing all the SeqRecords generated from a
        single read.

    :type mapped_fragments: list
    :param mapped_fragments: A list containing the Fr.Ids of the mapped
        fragments.

    :rtype list
    :return: A list containing the new SeqRecords.
    """
    new_seqs = list()

    if len(fasta_records) < 2:
        # no fragments to merge
        return new_seqs

    for fragment_idx in range(1, len(fasta_records)):
        if fragment_idx in mapped_fragments or fragment_idx - 1 in mapped_fragments:
            # not both unmapped
            continue

        prev_fragment = fasta_records[fragment_idx - 1]
        cur_fragment = fasta_records[fragment_idx]

        seq = prev_fragment.seq + cur_fragment.seq
        new_seqs.append(SeqRecord.SeqRecord(
            seq,
            id="Src.Op:MergeMap;Src.Fr:{:d},{:d}".format(fragment_idx - 1, fragment_idx),
            description="", name=""
        ))

    return new_seqs


if __name__ == "__main__":
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("input_fasta", help="Input FASTA file.", metavar="FASTA", action="store", type=str)
    args = parser.parse_args()

    mapped_fragments = get_mapped_fragments(args.input_sam)
    combine_unmapped_and_mapped(args.input_fasta, mapped_fragments)

    logging.shutdown()
