#!/bin/sh

for arg in "$@"; do
	zcat "${arg}" \
		| bioawk -c fastx '{print ">" $name "_start " $comment "\n" substr($seq, 0, 50) "\n" ">" $name "_end " $comment "\n" substr($seq, length($seq)-49, length($seq));}' \
		| gzip \
		> "primers_${arg%%.*}.fa.gz"
done
