#!/usr/bin/env python3

from os import path
import sys

from mchc_tools import get_vp_info

# Initialization
run_nid = int(sys.argv[1])
vp_info = get_vp_info(run_nid)

# Map split fragments to genome
fasta_fname = './split_files/spf_' + vp_info['run_id'] + '.fasta.gz'
bam_fname = './bam_files/bam_' + vp_info['run_id'] + '.bam'

# '~/bulk/env/bwa/bwa mem -k 15 -T 15 -t 12 ' +
cmd_str = \
    'bwa bwasw -b 5 -q 2 -r 1 -z 10 -T 15 -t 6 ' + \
    '/data0/repository/{0:s}/{0:s}.fa {1:s} '.format(vp_info['genome'], fasta_fname) + \
    '| samtools view -q 1 -hbS - | samtools sort ' + \
    '> ' + bam_fname
print('>>> You can run:\n$ {:s}'.format(cmd_str))

# Make an index file
# bwa index ./chrAll.fa

# Delft: sinter bigmem --ntasks=1 --cpus-per-task=64  bash
# UMC: qlogin -l h_rt=08:00:00 -l h_vmem=50G -pe threaded 1
# qsub -N <job_name> -l h_rt=08:00:00 -l h_vmem=50G -pe threaded 1 ~/bulk/bin/run_script.sh "../../Useful_Sample_Codes/BWA/bwa/bwa bwasw a.fasta > b.sam";

# Reason for BWA-SW: Split-mapped fragments have flag=0x100 (=2048) as set
# [main] Version: 0.7.16a-r1181
# [main] CMD: /Users/aallahyar/mw/Useful_Sample_Codes/BWA/bwa/bwa bwasw /Users/aallahyar/Technical/Dataset/Genome_Assembly/Mus_Musculus/mm9/chrAll.fa_BWA ./seq.fasta
# [main] Real time: 4.295 sec; CPU: 4.012 sec
# read1	16	chr7	110976097	15	280S166M3I181M
# read1	16	chr7	110975046	119	10M1I42M2D228M349S
#
# [main] CMD: /Users/aallahyar/mw/Useful_Sample_Codes/BWA/bwa/bwa mem /Users/aallahyar/Technical/Dataset/Genome_Assembly/Mus_Musculus/mm9/chrAll.fa_BWA ./seq.fasta
# [main] Real time: 4.663 sec; CPU: 4.184 sec
# read1	16	chr7	110976097	33	280S166M3I181M
# read1	2064	chr7	110975046	33	10M1I42M2D228M349H
