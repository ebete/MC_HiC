#!/usr/bin/env python3

import argparse
import gzip
import logging

from Bio import SeqIO, SeqRecord
from Bio.Alphabet import IUPAC

import utils


def combine_unmapped_and_mapped(fasta_file, mapped_fragments, max_fragment_slice_length):
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
            if last_id in mapped_fragments:
                fragment_id = 0
                for new_record in do_split(read_records, mapped_fragments[last_id],
                                           add_length_cutoff=max_fragment_slice_length):
                    new_record.id = "Fq.Id:{:s};Rd.Id:{:s};Rd.Ln:{:s};Fr.Id:{:d};Fr.Ln:{:d};{:s}".format(
                        last_metadata["Fq.Id"], last_metadata["Rd.Id"], last_metadata["Rd.Ln"], fragment_id,
                        len(new_record.seq), new_record.id)
                    print(new_record.format("fasta"), end="")
                    fragment_id += 1
            else:
                unmapped_reads += 1

            last_id = fasta_id
            last_metadata = record_metadata
            read_records = [record]
    logging.info("Out of %d reads, %d had no mapped fragments (%.1f%%).", total_reads, unmapped_reads,
                 unmapped_reads / total_reads * 100)


def do_split(fasta_records, mapped_fragments, add_length_cutoff=-1):
    """
    Create new fragments by splitting the clipped part from fragments and
    merging them with the next clipped part. Also create fragments by taking
    the end parts of unmapped fragments and merging them together.

    :type fasta_records: list
    :param fasta_records: A list containing all the SeqRecords generated from a
        single read. The list is ordered based on fragment ID (i.e. index 3 is
        also fragment 3).

    :type mapped_fragments: dict
    :param mapped_fragments: A dictionary containing the Fr.Ids of the mapped
        fragments and the associated metadata.

    :type add_length_cutoff: int
    :param add_length_cutoff: Limit the maximum size of the part taken from a
        fragment.

    :rtype list
    :return: A list containing the new SeqRecords.
    """
    new_seqs = list()

    mapped_frid = set(mapped_fragments.keys())
    total_fragments = len(fasta_records)
    for idx in range(2, total_fragments - 1):  # skip first two and last fragment
        prev_part = fasta_records[idx - 1].seq[-add_length_cutoff:] if add_length_cutoff > 0 else fasta_records[
            idx - 1].seq
        cur_part = fasta_records[idx].seq[:add_length_cutoff] if add_length_cutoff > 0 else fasta_records[idx].seq
        operation = ""
        new_seq = None
        unmapped_added = -1

        if idx - 1 in mapped_frid and idx in mapped_frid:
            # both fragments are mapped
            operation = "BothMapped"
            new_seq = prev_part + cur_part
            unmapped_added = 0
        elif idx - 1 in mapped_frid:
            # only previous fragment mapped
            operation = "LeftMapped"
            new_seq = fasta_records[idx - 1].seq + cur_part
            unmapped_added = len(cur_part)
        elif idx in mapped_frid:
            # only current fragment mapped
            operation = "RightMapped"
            new_seq = prev_part + fasta_records[idx].seq
            unmapped_added = len(prev_part)
        else:
            # both unmapped
            operation = "BothUnmapped"
            new_seq = fasta_records[idx - 1].seq + fasta_records[idx].seq
            unmapped_added = len(new_seq)

        new_seqs.append(SeqRecord.SeqRecord(
            seq=new_seq,
            id="Src.Op:{:s};Src.Fr:{:d},{:d};Src.Ln:{:d}".format(operation, idx - 1, idx, unmapped_added),
            description="", name=""
        ))

    return new_seqs


if __name__ == "__main__":
    utils.init_logger()

    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam", help="Input SAM/BAM file.", metavar="SAM", action="store", type=str)
    parser.add_argument("input_fasta", help="Input FASTA file.", metavar="FASTA", action="store", type=str)
    parser.add_argument("-m", "--max-slice-len", help="Maximum size of the part taken from a fragment.", metavar="N",
                        action="store", type=int, default=-1)
    args = parser.parse_args()

    mapped_fragments = utils.get_mapping_metadata(args.input_sam)
    combine_unmapped_and_mapped(args.input_fasta, mapped_fragments, args.max_slice_len)

    logging.shutdown()
