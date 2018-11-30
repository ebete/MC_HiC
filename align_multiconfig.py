#!/usr/bin/env python3

from __future__ import print_function

import argparse
import csv
import logging
import os
import subprocess
import sys

# add current conda environment to the search path
_exec_env = os.environ.copy()
_exec_env["PATH"] = "{:s}:{:s}".format(os.path.dirname(sys.executable), _exec_env["PATH"])


def load_configs(config_file):
    """
    Parse alignment options from a TSV file.

    :type config_file: str
    :param config_file: Path to the TSV file containing the alignment configuration

    :rtype dict
    :return: A dictionary with the alignment settings
    """
    logging.info("Parsing configuration file %s ...", config_file)
    cfg = {}
    with open(config_file, "r", newline="") as fin:
        handle = csv.reader(fin, delimiter=";")
        colnames = [str(x).strip().lower() for x in next(handle, [])]
        for row in handle:
            try:
                if not row or len(row) != len(colnames) or row[0][0] == "#":
                    continue
            except IndexError:  # empty first column
                continue
            # create key, value pairs of the csv headers and the parameter values
            cfg[row[0]] = {k: str(v).strip() for k, v in zip(colnames[1:], row[1:])}
            # Find and add all the paths to the FASTA files
            cfg[row[0]]["input_fasta"] = get_fasta_splits(cfg[row[0]]["input_fasta"])

            logging.debug("Found config %s: %s", row[0], cfg[row[0]])
    return cfg


def run_aligners(config, output_dir, force_run):
    """
    Run the aligners using the given configuration.

    :type config: dict
    :param config: Dictionary containing the run parameters

    :type output_dir: str
    :param output_dir: Path to write the resulting BAM files to
    """
    for cfg_name, cfg_params in config.items():
        logging.info("Running configuration %s ...", cfg_name)
        fasta_input = cfg_params["input_fasta"][0]

        bam_out = os.path.join(output_dir, "{}_{}.bam".format(cfg_name, cfg_params["aligner"]))
        if not force_run and os.path.exists(bam_out):
            logging.warning("File %s exists; %s will be skipped. Use --force-run to override this check.", bam_out,
                            cfg_name)
            continue

        if "bwa" in cfg_params["aligner"]:
            run_bwa(fasta_input, bam_out, cfg_params)
        elif "bowtie" in cfg_params["aligner"]:
            run_bowtie2(fasta_input, bam_out, cfg_params)
        elif "last" in cfg_params["aligner"]:
            run_last(fasta_input, bam_out, cfg_params)
        else:
            logging.warning("Unknown aligner %s specified. Skipping %s ...", cfg_params["aligner"], cfg_name)


def run_bwa(in_fasta, out_bam, params):
    logging.info("Executing BWA ...")
    # @formatter:off
    align_out = subprocess.Popen(("bwa", "mem",
        "-x", "ont2d",  # Oxford Nanopore 2D reads
        "-t", params["threads"],  # number of processing threads
        "-k", params["seed_length"],  # minimum seed length
        "-A", params["sw_match"],  # match score
        "-B", params["sw_mismatch"],  # mismatch score
        "-O", params["query_gap_open"],  # gap open penalty
        "-E", params["query_gap_extend"],  # gap extension penalty
        "-w", params["max_gap_length"],  # maximum gap length
        "-T", params["score_threshold"],  # output alignment score threshold
        params["reference"],  # reference genome
        in_fasta  # FASTA with reads
    ), stdout=subprocess.PIPE, env=_exec_env)
    # @formatter:on
    print(" ".join(align_out.args))

    export_to_bam(align_out.stdout, out_bam)
    align_out.wait()

    if align_out.returncode != 0:
        raise subprocess.SubprocessError("BWA exited with a non-zero status ({})".format(align_out.returncode))
    logging.info("BWA alignment file written to %s", out_bam)


def run_bowtie2(in_fasta, out_bam, params):
    logging.info("Executing Bowtie2 ...")
    # @formatter:off
    # functions: (C)onstant (L)inear (S)qrt lo(G)_e
    align_out = subprocess.Popen(("bowtie2",
        "--threads", params["threads"],  # number of processing threads
        "-x", params["reference"],  # reference genome
        "-U", in_fasta, "-q" if "fastq" in in_fasta.split(".") else "-f",  # FASTA/Q with reads
        "--local",  # use Smith-Waterman alignment
        "-D", params["consec_seed_fails"],  # consecutive seed extensions that may fail
        "-R", params["max_reseeds"],  # max re-seeds allowed
        "-N", params["seed_mismatch"],  # max mismatches in seeds
        "-L", params["seed_length"],  # length of seeds
        "-i", params["seed_interval_fun"],  # function used to calculate intervals between seeds
        "--ma", params["sw_match"],  # match bonus during alignment
        "--mp", "{0:},{0:}".format(params["sw_mismatch"]),  # mismatch penalty during alignment
        "--score-min", params["min_score_fun"],  # function used to calculate miniumum alignment score
        "--rdg", "{},{}".format(params["query_gap_open"], params["query_gap_extend"]),  # read gap open/extension penalties
        "--rfg", "{},{}".format(params["ref_gap_open"], params["ref_gap_extend"]),  # reference gap open/extension penalties
        "--ignore-quals",  # ignore Phred scores for mismatch scoring
        #"-a",  # report all valid alignments (very slow)
        "-t"  # write run times to stderr
    ), stdout=subprocess.PIPE, env=_exec_env)
    # @formatter:on
    print(" ".join(align_out.args))

    export_to_bam(align_out.stdout, out_bam)
    align_out.wait()

    if align_out.returncode != 0:
        raise subprocess.SubprocessError("Bowtie2 exited with a non-zero status ({})".format(align_out.returncode))
    logging.info("Bowtie2 alignment file written to %s", out_bam)


def run_last(in_fasta, out_bam, params):
    logging.info("Executing LAST ...")
    # LAST aligner
    with open(out_bam + ".maf", "w") as tmp_maf:
        # @formatter:off
        last_out = subprocess.Popen(("lastal",
            "-P", params["threads"],  # number of processing threads
            "-a", params["query_gap_open"],  # gap open penalty
            "-b", params["query_gap_extend"],  # gap extension penalty
            "-A", params["ref_gap_open"],  # reference gap open penalty
            "-B", params["ref_gap_extend"],  # reference gap extension penalty
            "-l", params["seed_length"],  # minimum length of seeds
            "-T", "0",  # local alignment
            "-s", "2",  # look on both strands
            "-Q", "1" if "fastq" in in_fasta.split(".") else "0",  # input data type
            "-e", params["score_threshold"],  # minimum alignment score
            params["reference"],  # reference genome
            in_fasta  # FASTA with reads
        ), stdout=tmp_maf, env=_exec_env)
        # @formatter:on
        print(" ".join(last_out.args))
        last_out.wait()
    # estimate split alignments
    with open(out_bam + ".maf", "r") as tmp_maf:
        # @formatter:off
        split_out = subprocess.Popen(("last-split",
            "-d", "2"  # strandedness is unknown
        ), stdin=tmp_maf, stdout=subprocess.PIPE, env=_exec_env)
        # @formatter:on
        print(" ".join(split_out.args))
        # convert MAF+ to SAM
        maf_out = subprocess.Popen(("maf-convert", "-n", "sam", "-"), stdin=split_out.stdout, stdout=subprocess.PIPE,
                                   env=_exec_env)
        print(" ".join(maf_out.args))
        # add @SQ header to SAM (reference sequence)
        samtools_out = subprocess.Popen(("samtools", "view", "-bt", params["reference"] + ".fai", "-"),
                                        stdin=maf_out.stdout,
                                        stdout=subprocess.PIPE, env=_exec_env)
        print(" ".join(samtools_out.args))

        export_to_bam(samtools_out.stdout, out_bam)

        if last_out.returncode != 0:
            raise subprocess.SubprocessError("LAST exited with a non-zero status ({})".format(last_out.returncode))
    logging.info("LAST alignment file written to %s", out_bam)


def export_to_bam(mapper_output, out_bam):
    """
    Takes the output stream of an aligner and writes it to a sorted+indexed BAM file.

    :type mapper_output: subprocess.IO
    :param mapper_output: Stdout stream of an aligner

    :type out_bam: str
    :param out_bam: Path of the output BAM file
    """
    with open(out_bam, "wb") as bamfile:
        bam_unsorted = subprocess.Popen(("samtools", "view", "-b", "-"), stdin=mapper_output, stdout=subprocess.PIPE,
                                        env=_exec_env)
        print(" ".join(bam_unsorted.args))
        bam_sorted = subprocess.Popen(("samtools", "sort", "-"), stdin=bam_unsorted.stdout, stdout=bamfile,
                                      env=_exec_env)
        print(" ".join(bam_sorted.args))
        bam_sorted.wait()
        bam_index = subprocess.Popen(("samtools", "index", out_bam), stdout=subprocess.DEVNULL, env=_exec_env)
        print(" ".join(bam_index.args))
        bam_index.wait()


def get_fasta_splits(split_path):
    """
    Finds the splitted FASTA files that match the path. If no splits are found,
    it will return the given path again.

    :type split_path: str
    :param split_path: Path to the splitted FASTA files

    :rtype list
    :return: All found paths
    """
    fasta_splits = []
    n = 0
    while True:
        fasta_path = os.path.join(os.path.dirname(split_path),
                                  "sub_{:d}_{:s}".format(n, os.path.basename(split_path)))
        if not os.path.isfile(fasta_path):
            break
        fasta_splits.append(fasta_path)
        n += 1
    return fasta_splits if len(fasta_splits) > 0 else [split_path]


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s] %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_config", help="Input configuration file", metavar="CONFIG", action="store", type=str)
    parser.add_argument("--output", "-o", help="BAM output directory", metavar="DIR", action="store", type=str,
                        default="./")
    parser.add_argument("--force-run", "-f", help="Force-overwrite existing BAM result files", action="store_true")
    args = parser.parse_args()

    cfg = load_configs(args.input_config)
    run_aligners(cfg, args.output, args.force_run)

    logging.shutdown()
