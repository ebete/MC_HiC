#!/bin/bash

for arg in "$@"; do
	OUT_COUNTS="${arg%%.*}_2mer.csv"
	OUT_MTX="${arg%%.*}_mtx.csv"
	TMP_JF="$(mktemp)"

	jellyfish count -t 4 -m 2 -s 100 -o "${TMP_JF}" <(gzip -cd "${arg}")
	jellyfish dump "${TMP_JF}"  \
		| bioawk -c fastx 'BEGIN{OFS="\t"; print "first","second","count";} {print substr($seq, 1, 1),substr($seq, 2, 1),$name;}' \
		> "${OUT_COUNTS}"

	python3 "/home/thom/PycharmProjects/McHiC/calculate_2mer_freq.py" "${OUT_COUNTS}" \
		> ${OUT_MTX}
	echo "Created 2-mer counts ${OUT_COUNTS} and matrix ${OUT_MTX}"

	rm "${TMP_JF}"
done
