#!/usr/bin/env python3
# compress raw: tar -czf ESCs-V122-PM_raw.tar.gz  ./deLaatUDNHiC5Pro3
# check integrity: tar -tzf my_tar.tar.gz >/dev/null
# combine fasq: cat *.fastq | gzip > ../../ESCs-V121-PM.fastq.gz
# run: ri=19; qsub -N splt${ri} -l h_rt=40:00:00 -l h_vmem=20G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 ./01_combine_raw_files.py ${ri}"

import argparse
import gzip
# Initialization
from glob import glob

# Get command argument
parser = argparse.ArgumentParser()
parser.add_argument("input_fastq", help="Input FASTQ file (gzipped).", metavar="INFILE", action="store", type=str)
parser.add_argument("output_fasta", help="Output FASTA file (gzipped).", metavar="OUTFILE", action="store", type=str)
args = parser.parse_args()

# Read gz files
print('Writing sequencing data to: {:s}'.format(args.output_fasta))

with gzip.open(args.output_fasta, 'wt') as fastq_fid:
    for fastq_ind, fastq_id in enumerate(vp_info['seq_file_indices'].split(';')):
        raw_ptr = './raw_files/raw_' + fastq_id + '_*.fastq.gz'
        raw_flst = glob(raw_ptr)
        if len(raw_flst) != 1:
            Exception('Can not find a unique source file: {:s}.'.format(raw_ptr))

        print('Reading: {:s}'.format(raw_flst[0]))
        rd_ind = 0
        with gzip.open(raw_flst[0], 'rt') as raw_fid:
            while True:
                rd_ind = rd_ind + 1
                rd_oid = raw_fid.readline().rstrip('\n')
                rd_seq = raw_fid.readline().rstrip('\n')
                rd_plus = raw_fid.readline().rstrip('\n')
                rd_pred = raw_fid.readline().rstrip('\n')
                if rd_oid == '':
                    break
                if rd_oid[0] != '@' or rd_plus != '+':
                    raise Exception('The file is corrupted.\n' +
                                    'Read #{:d}:\n'.format(rd_ind) +
                                    '\tID: [{:s}],\n\tplus: [{:s}]'.format(rd_oid, rd_plus))
                if rd_ind % 10000 == 0:
                    print('Processed {:d} reads.'.format(rd_ind))

                rd_sid = 'Fl.Id:{:d};Rd.Id:{:d};Rd.Ln:{:d}'.format(fastq_ind + 1, rd_ind, len(rd_seq))
                fastq_fid.write('>' + rd_sid + '\n')
                fastq_fid.write(rd_seq + '\n')
