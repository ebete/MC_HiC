#!/bin/sh

# BED header
echo 'track name="covered_regions" description="Regions with a decent coverage"'

tmpfile="$(mktemp)"

>&2 echo "Calculating depth ..."
samtools depth -d 0 "$1" \
	| awk -F'\t' '$3>=3' \
	> "$tmpfile"

>&2 echo "Writing regions in BED format ..."
python3 /home/thom/PycharmProjects/McHiC/collapse_coverage.py -d 10 -l 50 "$tmpfile" \
	| awk -F'[\:\-\\t]' 'BEGIN{OFS="\t";} NR>1{print $1,$2,$3;}'

rm "$tmpfile"
