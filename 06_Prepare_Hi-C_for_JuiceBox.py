#!/usr/bin/env python

# run: qsub -N HiCProc -l h_rt=50:00:00 -l h_vmem=50G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 06_Prepare_Hi-C_for_JuiceBox.py";

import gzip
import numpy as np
import os
from seq_tools import get_chr_info
from _utilities import accum_array
from mchc_tools import get_vp_info, get_mchc_data

# Initialization
run_list = [0]
vp_info_lst = [get_vp_info(vpi) for vpi in run_list]
# run_id = ','.join([vp_info_lst[i]['run_id'] for i in range(len(vp_info_lst))])
run_id = '{:s}-{:s}_R{:d}-{:d}'.format(vp_info_lst[0]['run_id'], vp_info_lst[-1]['run_id'], run_list[0], run_list[-1])
chr_lst = get_chr_info(vp_info_lst[0]['genome'], 'chr_name')

# Load data
frg_pd = get_mchc_data(vp_info_lst, min_MQ=20)
print('Sorting circles')
frg_pd = frg_pd.iloc[np.lexsort([frg_pd['SeqStart'], frg_pd['Chr']])]

# Output source files
hc_fname = './hic_files/hic_' + run_id + '.tsv.gz'
with gzip.open(hc_fname, 'wb') as hc_fid:
    cir_uid = np.unique(frg_pd['ReadID'].values, return_inverse=True)[1]
    cir_grp = accum_array(cir_uid, frg_pd.values)
    n_cir = len(cir_grp)
    print('Looping over [{:d}] circles'.format(n_cir))
    for cir_idx, cir in enumerate(cir_grp):
        if cir_idx % 100000 == 0:
            print('{:d}/{:d} reads are processed.'.format(cir_idx, n_cir))
        grp_size = len(cir)
        if grp_size <= 0:
            continue
        out_str = ''
        for f1 in range(grp_size):
            for f2 in range(f1+1, grp_size):
                out_str = out_str + '\t'.join([
                    'Cir{:07d}'.format(cir[f1, 0]),     # readname
                    '{:d}'.format(cir[f1, 4] == 1),     # strand1
                    '{:d}'.format(cir[f1, 1]),          # chr1
                    '{:d}'.format(cir[f1, 2]),          # pos1
                    '1',                                # frag1
                    '{:d}'.format(cir[f2, 4] == 1),     # strand2
                    '{:d}'.format(cir[f2, 1]),          # chr2
                    '{:d}'.format(cir[f2, 2]),          # pos2
                    '2',                                # frag2
                    '{:d}'.format(cir[f1, 5]),          # MQ1
                    '{:d}'.format(cir[f2, 5]),          # MQ2
                ]) + '\n'
        if len(out_str) > 0:
            hc_fid.write(out_str)
print()

# Sorting gz file
sgz_fname = './hic_files/hic_' + run_id + '_sorted.tsv.gz'
cmd_str = 'zcat ' + hc_fname + ' | sort -V -S 80% -k3,3 -k7,7 -T ./ --parallel=12 | gzip > ' + sgz_fname  # Using 40% of memory for buffer

print('Running: ' + cmd_str)
os.system(cmd_str)

# Convert to .hic
hic_fname = './hic_files/hic_' + run_id + '.hic'
cmd_str = 'java -Xmx10g -jar ./hic_files/juicer_tools.1.8.9_jcuda.0.8.jar pre ' + \
          sgz_fname + ' ' + hic_fname + ' ' + vp_info_lst[0]['genome']
print('Running: ' + cmd_str)
os.system(cmd_str)
