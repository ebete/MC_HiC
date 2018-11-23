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
        next(handle)  # skip header
        for row in handle:
            try:
                if not row or len(row) < 4 or row[0][0] == "#":
                    continue
            except IndexError:
                continue
            cfg[row[0]] = {
                "fasta": get_fasta_splits(row[1]),
                "ref": row[2],
                "args": tuple(row[3:])
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
        fasta_input = cfg_params["fasta"][0]

        logging.debug("%s %s", cfg_name, cfg_params)

        # BWA
        bam_out = os.path.join(output_dir, "{}_{}.bam".format(cfg_name, "bwa"))
        run_bwa(fasta_input, cfg_params["ref"], bam_out, cfg_params["args"])
        # Bowtie2
        bam_out = os.path.join(output_dir, "{}_{}.bam".format(cfg_name, "bowtie2"))
        run_bowtie2(fasta_input, cfg_params["ref"], bam_out, cfg_params["args"])


def run_bwa(in_fasta, in_reference, out_bam, params):
    logging.info("Executing BWA ...")
    align_out = subprocess.Popen(
        ("bwa", "mem", "-t", params[0], "-T", params[1], "-C", in_reference, in_fasta),
        stdout=subprocess.PIPE, env=_exec_env)
    logging.info(" ".join(align_out.args))

    export_to_bam(align_out.stdout, out_bam)
    align_out.wait()

    if align_out.returncode != 0:
        raise subprocess.SubprocessError("BWA exited with a non-zero status ({})".format(align_out.returncode))


def run_bowtie2(in_fasta, in_reference, out_bam, params):
    logging.info("Executing Bowtie2 ...")
    align_out = subprocess.Popen(
        ("bowtie2", "--threads", params[0], "-x", in_reference, "-U", in_fasta, "-f", "--local",
         "--very-sensitive-local", "-a", "-t"),
        stdout=subprocess.PIPE, env=_exec_env)
    logging.info(" ".join(align_out.args))

    export_to_bam(align_out.stdout, out_bam)
    align_out.wait()

    if align_out.returncode != 0:
        raise subprocess.SubprocessError("Bowtie2 exited with a non-zero status ({})".format(align_out.returncode))


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
