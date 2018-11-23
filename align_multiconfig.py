#!/usr/bin/env python3

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
        handle = csv.reader(fin, delimiter="\t")
        colnames = next(handle, None)
        for row in handle:
            try:
                if not row or len(row) < 4 or row[0][0] == "#":
                    continue
            except IndexError:
                continue
            cfg[row[0]] = {
                "fasta": get_fasta_splits(row[1]),
                "ref": row[2],
                # create key, value pairs of the csv headers and the parameter values
                "args": {k: v for k, v in zip(colnames[3:], row[3:])}
            }
            logging.debug("Found config %s: %s", row[0], cfg[row[0]])
    return cfg


def run_aligners(config, output_dir):
    """
    Run the aligners using the given configuration.

    :type config: dict
    :param config: Dictionary containing the run parameters

    :type output_dir: str
    :param output_dir: Path to write the resulting BAM files to
    """
    for cfg_name, cfg_params in config.items():
        logging.info("Running configuration %s ...", cfg_name)
        fasta_input = cfg_params["fasta"][0]

        # BWA
        bam_out = os.path.join(output_dir, "{}_{}.bam".format(cfg_name, "bwa"))
        run_bwa(fasta_input, cfg_params["ref"], bam_out, cfg_params["args"])
        # Bowtie2
        bam_out = os.path.join(output_dir, "{}_{}.bam".format(cfg_name, "bowtie2"))
        run_bowtie2(fasta_input, cfg_params["ref"], bam_out, cfg_params["args"])
        # LAST
        bam_out = os.path.join(output_dir, "{}_{}.bam".format(cfg_name, "last"))
        run_last(fasta_input, cfg_params["ref"], bam_out, cfg_params["args"])


def run_bwa(in_fasta, in_reference, out_bam, params):
    logging.info("Executing BWA ...")
    # @formatter:off
    align_out = subprocess.Popen(("bwa", "mem",
        "-t", params["threads"],  # number of processing threads
        "-T", params["score_threshold"],  # output alignment score threshold
        "-C", in_reference,  # reference genome
        in_fasta  # FASTA with reads
    ), stdout=subprocess.PIPE, env=_exec_env)
    # @formatter:on
    logging.info(" ".join(align_out.args))

    export_to_bam(align_out.stdout, out_bam)
    align_out.wait()

    if align_out.returncode != 0:
        raise subprocess.SubprocessError("BWA exited with a non-zero status ({})".format(align_out.returncode))
    logging.info("BWA alignment file written to %s", out_bam)


def run_bowtie2(in_fasta, in_reference, out_bam, params):
    logging.info("Executing Bowtie2 ...")
    # @formatter:off
    # functions: (C)onstant (L)inear (S)qrt lo(G)_e
    align_out = subprocess.Popen(("bowtie2",
        "--threads", params["threads"],  # number of processing threads
        "-x", in_reference,  # reference genome
        "-U", in_fasta, "-f",  # FASTA with reads
        "--local",  # use Smith-Waterman alignment
        "-D", params["consec_seed_fails"],  # consecutive seed extensions that may fail
        "-R", params["max_reseeds"],  # max re-seeds allowed
        "-N", params["seed_mismatch"],  # max mismatches in seeds
        "-L", params["seed_length"],  # length of seeds
        "-i", params["seed_interval_fun"],  # function used to calculate intervals between seeds
        "--ma", params["sw_match"],  # match bonus during alignment
        "--score-min", params["min_score_fun"],  # function used to calculate miniumum alignment score
        "--rdg", params["query_gap"],  # read gap open/extension penalties
        "--rfg", params["ref_gap"],  # reference gap open/extension penalties
        "-a",  # report all valid alignments (very slow)
        "-t",  # write run times to stderr
    ), stdout=subprocess.PIPE, env=_exec_env)
    # @formatter:on
    logging.info(" ".join(align_out.args))

    export_to_bam(align_out.stdout, out_bam)
    align_out.wait()

    if align_out.returncode != 0:
        raise subprocess.SubprocessError("Bowtie2 exited with a non-zero status ({})".format(align_out.returncode))
    logging.info("Bowtie2 alignment file written to %s", out_bam)


def run_last(in_fasta, in_reference, out_bam, params):
    logging.info("Executing LAST ...")
    # LAST aligner
    # @formatter:off
    last_out = subprocess.Popen(("lastal",
        "-P", params["threads"],  # number of processing threads
        in_reference,  # reference genome
        in_fasta,  # FASTA with reads
    ), stdout=subprocess.PIPE, env=_exec_env)
    # @formatter:on
    logging.info(" ".join(last_out.args))
    # estimate split alignments
    split_out = subprocess.Popen(("last-split"), stdin=last_out.stdout, stdout=subprocess.PIPE, env=_exec_env)
    logging.info(" ".join(split_out.args))
    # convert MAF+ to SAM
    maf_out = subprocess.Popen(("maf-convert", "-n", "sam", "-"), stdin=split_out.stdout, stdout=subprocess.PIPE,
                               env=_exec_env)
    logging.info(" ".join(maf_out.args))
    # add @SQ header to SAM (reference sequence)
    samtools_out = subprocess.Popen(("samtools", "view", "-bt", in_reference + ".fai", "-"), stdin=maf_out.stdout,
                                    stdout=subprocess.PIPE, env=_exec_env)
    logging.info(" ".join(samtools_out.args))

    export_to_bam(samtools_out.stdout, out_bam)
    last_out.wait()

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
        logging.info(" ".join(bam_unsorted.args))
        bam_sorted = subprocess.Popen(("samtools", "sort", "-"), stdin=bam_unsorted.stdout, stdout=bamfile,
                                      env=_exec_env)
        logging.info(" ".join(bam_sorted.args))
        bam_sorted.wait()
        bam_index = subprocess.Popen(("samtools", "index", out_bam), stdout=subprocess.DEVNULL, env=_exec_env)
        logging.info(" ".join(bam_index.args))
        bam_index.wait()


def get_fasta_splits(split_path):
    """
    Finds the splitted FASTA files that match the path.

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
    return fasta_splits


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s] %(message)s")

    # Get command argument
    parser = argparse.ArgumentParser()
    parser.add_argument("input_config", help="Input configuration file", metavar="CONFIG", action="store", type=str)
    parser.add_argument("--output", "-o", help="BAM output directory", metavar="DIR", action="store", type=str,
                        default="./")
    args = parser.parse_args()

    cfg = load_configs(args.input_config)
    run_aligners(cfg, args.output)

    logging.shutdown()
