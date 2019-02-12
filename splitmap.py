#!/usr/bin/env python3

import argparse
import gzip
import logging

import pysam
from Bio import SeqIO, SeqRecord
from Bio.Alphabet import IUPAC

import utils


def get_mapped_fragments(sam_input):
    """
    Generate a dictionary containing the metadata of the mapped read fragments
    from a SAM file. Only the primary alignment is considered.

    The dictionary is structured as follows:

    {
      read_id: {
        fragment_id: {
          metadata_key: metadata_value
        }
      }
    }

    :type sam_input str
    :param sam_input: The SAM file to parse all the mapped fragments from.

    :rtype dict
    :return: A dictionary containing mapped reads metadata.
    """
    alignments = dict()

    logging.info("Reading SAM file %s ...", sam_input)
    with pysam.AlignmentFile(sam_input, "r") as sam:
        for read in sam:
            if read.is_unmapped or read.is_secondary:
                continue

            metadata = utils.read_header_to_dict(read.qname)
            metadata["cigar"] = read.cigar
            frid = int(metadata["Fr.Id"])

            alignments.setdefault(utils.make_read_id(metadata), dict())[frid] = metadata

    return alignments


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
    for idx in range(2, total_fragments - 1):  # skip first two and last
        if idx - 1 in mapped_frid and idx in mapped_frid:  # run clipped merge
            prev_fragment = mapped_fragments[idx - 1]
            cur_fragment = mapped_fragments[idx]

            prev_clipped = prev_fragment["cigar"][-1][1] if prev_fragment["cigar"][-1][0] in {4, 5} else 0
            cur_clipped = cur_fragment["cigar"][0][1] if cur_fragment["cigar"][0][0] in {4, 5} else 0
            if prev_clipped < 20 or cur_clipped < 20:  # clipped sequences too short
                continue

            # restrict maximum cut size
            prev_clipped = min(add_length_cutoff, prev_clipped) if add_length_cutoff > 0 else prev_clipped
            cur_clipped = min(add_length_cutoff, cur_clipped) if add_length_cutoff > 0 else cur_clipped

            # create new fragment by merging clipped parts
            new_seqs.append(SeqRecord.SeqRecord(
                seq=fasta_records[idx - 1].seq[-prev_clipped:] + fasta_records[idx].seq[:cur_clipped],
                id="Src.Op:MergeClipped;Src.Fr:{:d},{:d}".format(idx - 1, idx),
                description="", name=""
            ))
            continue

        # a fragment was not mapped: make small fragment from ends
        seq = None
        if add_length_cutoff > 0:
            seq = fasta_records[idx - 1].seq[-add_length_cutoff:] + fasta_records[idx].seq[:add_length_cutoff]
        else:
            seq = fasta_records[idx - 1].seq + fasta_records[idx].seq

        new_seqs.append(SeqRecord.SeqRecord(
            seq=seq,
            id="Src.Op:MergeEnds;Src.Fr:{:d},{:d}".format(idx - 1, idx),
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

    mapped_fragments = get_mapped_fragments(args.input_sam)
    combine_unmapped_and_mapped(args.input_fasta, mapped_fragments, args.max_slice_len)

    logging.shutdown()