#!/usr/bin/env python3
# ri=19; qsub -N splt${ri} -l h_rt=40:00:00 -l h_vmem=20G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 ./02_restriction_site_split.py ${ri}"

import re
import gzip
from os import path
import sys

from mchc_tools import get_vp_info
from seq_tools import get_re_info

# Initialization
run_nid = int(sys.argv[1])
MAX_FRG_SIZE = 2000

# Load vp info file
vp_info = get_vp_info(run_nid)
re_seq = [get_re_info(re_name=re_name, property='seq') for re_name in vp_info['res_cutters'].split(';')]
regex_ptr = '|'.join(re_seq)

# Read the the fasta file
fasta_name = './fasta_files/fa_' + vp_info['run_id'] + '.fasta.gz'
splt_name = './split_files/spf_' + vp_info['run_id'] + '.fasta.gz'
rd_ind = 1
frg_ind = 1
n_reduced = 0
print('Reading reads from: {:s}'.format(fasta_name))
print('Writing fragments to: {:s}'.format(splt_name))
assert not path.isfile(splt_name)
with gzip.open(fasta_name, 'rt') as fasta_fid, \
        gzip.open(splt_name, 'wt') as splt_fid:
    while True:
        rd_sid = fasta_fid.readline().rstrip('\n')
        rd_seq = fasta_fid.readline().rstrip('\n')
        if rd_sid == '':
            break
        if rd_ind % 10000 == 0:
            print('Processed {:d} reads and produced {:d} fragments.'.format(rd_ind, frg_ind))

        frg_be = 0
        for re_item in re.finditer(regex_ptr, rd_seq):
            frg_en = re_item.end()
            if frg_en - frg_be > MAX_FRG_SIZE:
                n_reduced += 1
                frg_en = MAX_FRG_SIZE
            splt_fid.write('{:s};Fr.Id:{:d};Fr.SBp:{:d};Fr.EBp:{:d}\n'.format(rd_sid, frg_ind, frg_be, frg_en))
            splt_fid.write(rd_seq[frg_be:frg_en] + '\n')
            frg_ind = frg_ind + 1
            frg_be = re_item.end() - len(re_item.group())
        frg_en = len(rd_seq)

        if frg_en - frg_be > MAX_FRG_SIZE:
            n_reduced += 1
            frg_en = MAX_FRG_SIZE
        splt_fid.write('{:s};Fr.Id:{:d};Fr.SBp:{:d};Fr.EBp:{:d}\n'.format(rd_sid, frg_ind, frg_be, frg_en))
        splt_fid.write(rd_seq[frg_be:frg_en] + '\n')
        frg_ind = frg_ind + 1
        rd_ind = rd_ind + 1

if n_reduced != 0:
    print('[w] Warning: {:d} fragments are reduced to {:d}bp.'.format(n_reduced, MAX_FRG_SIZE))
