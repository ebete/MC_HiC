#!/usr/bin/env python3

import argparse
import gzip
import logging

import pysam
from Bio import Seq, SeqIO, SeqRecord
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
            if last_id in mapped_fragments:
                mapped_fragments[last_id] = sorted(set(mapped_fragments[last_id]))
                fragment_id = 0
                for new_record in do_merge(read_records, mapped_fragments[last_id]):
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


def do_merge(fasta_records, mapped_fragments, add_length_cutoff=50):
    """
    When an unmapped region exists betweeen two mapped fragments, generate new
    fragments by extending both mapped fragments into the unmapped region. The
    extension will stop when it reaches the other mapped fragment or when the
    length threshold is exceeded. Only the longest fragments are returned from
    both extensions.

    :type fasta_records: list
    :param fasta_records: A list containing all the SeqRecords generated from a
        single read.

    :type mapped_fragments: list
    :param mapped_fragments: A list containing the Fr.Ids of the mapped
        fragments.

    :type add_length_cutoff: int
    :param add_length_cutoff: When the fragment is extended by more than this
        value, stop extending.

    :rtype list
    :return: A list containing the new SeqRecords.
    """
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

        # extend from left
        left_extend = None
        left_orig_len = len(fasta_records[jump_region[0]])
        for i in range(2, len(jump_region)):
            seq = Seq.Seq("".join([str(fasta_records[x].seq) for x in jump_region[:i]]), alphabet=IUPAC.unambiguous_dna)
            left_extend = SeqRecord.SeqRecord(
                seq[:left_orig_len + add_length_cutoff],
                id="Src.Op:MergeMap;Src.Fr:{:s};Src.Ori:Left".format(",".join([str(x) for x in jump_region[:i]])),
                description="", name=""
            )
            appended_length = len(left_extend) - left_orig_len
            left_extend.id = "{:s};Src.Ln:{:d}".format(left_extend.id, appended_length)

            if appended_length >= add_length_cutoff:
                break
        new_seqs.append(left_extend)

        # extend from right
        right_extend = None
        right_orig_len = len(fasta_records[jump_region[-1]])
        for i in range(len(jump_region) - 2, 0, -1):
            seq = Seq.Seq("".join([str(fasta_records[x].seq) for x in jump_region[i:]]), alphabet=IUPAC.unambiguous_dna)
            right_extend = SeqRecord.SeqRecord(
                seq[-(right_orig_len + add_length_cutoff):],
                id="Src.Op:MergeMap;Src.Fr:{:s};Src.Ori:Right".format(",".join([str(x) for x in jump_region[i:]])),
                description="", name=""
            )
            appended_length = len(right_extend) - right_orig_len
            right_extend.id = "{:s};Src.Ln:{:d}".format(right_extend.id, appended_length)

            if appended_length >= add_length_cutoff:
                break
        new_seqs.append(right_extend)

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
