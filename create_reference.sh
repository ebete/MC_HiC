#!/bin/bash -e

THREADS=4

# Bowtie2
bowtie2-build -f --threads "${THREADS}" "$1" "$1"

# BWA
bwa index -a bwtsw "$1"

# LAST
lastdb -S 2 -v "$1" "$1"
