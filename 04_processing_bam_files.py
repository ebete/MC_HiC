#!/usr/bin/env python2

# ri=10; qsub -N Prc${ri} -l h_rt=50:00:00 -l h_vmem=20G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 04_processing_bam_files.py ${ri}";

# Initialize
import h5py
import numpy as np
import gzip
import pysam
import pandas as pd
from pandas import read_csv
from os import path
import sys

from mchc_tools import get_vp_info
from seq_tools import get_chr_info, hasOL


def saveRead(gz_fid, frg_lst):
    # Check if each fragment is a true hop
    n_frg = frg_lst.shape[0]
    if n_frg > 1:
        for fi in range(n_frg):
            if frg_lst[fi, 11] == 0:
                continue
            nei_idx = np.where(hasOL(frg_lst[fi, 1:4], frg_lst[:, 1:4], offset=5000))[0]
            if len(nei_idx) > 1:
                top_MQi = np.argmax(frg_lst[nei_idx, 5])
                frg_lst[nei_idx, 11] = 0
                frg_lst[nei_idx[top_MQi], 11] = 1  # Setting the top fragment as true hop

    # Save fragment set
    for frg in frg_lst:
        gz_fid.write(frg_template.format(*frg))


# Main code
run_nid = int(sys.argv[1])
vp_info = get_vp_info(run_nid)
chr_lst = get_chr_info(vp_info['genome'], 'chr_name')
chr_map = dict(zip(chr_lst, np.arange(len(chr_lst)) + 1))

# Read file line by line
bam_fname = './bam_files/bam_{:s}.bam'.format(vp_info['run_id'])
gz_fname = './frag_files/frg_{:s}.tsv.gz'.format(vp_info['run_id'])
# assert not path.isfile(gz_fname)  TODO: uncomment
print('Processing: {:s}\n'.format(bam_fname) +
      'Saving to: {:s}\n'.format(gz_fname))
ReadID_old = -1
FrgID_old = -1
n_frg_info = 12
n_processed = 0
frg_template = '\t'.join(['{:d}'] * n_frg_info) + '\n'
frg_set = np.empty([0, n_frg_info], dtype=np.int64)
with pysam.AlignmentFile(bam_fname, 'rb') as bam_fid, gzip.open(gz_fname, 'wt') as gz_fid:
    gz_fid.write(
        '\t'.join(['ReadID', 'Chr', 'RefStart', 'RefEnd', 'Strand', 'MQ',
                   'FileID', 'FrgID', 'SeqStart', 'SeqEnd', 'ReadLength', 'TrueHop']) + '\n'
    )
    for que_idx, que_line in enumerate(bam_fid):
        if que_idx % 100000 == 0:
            print('Processed {:d} fragments in {:d} reads.'.format(que_idx + 1, n_processed))
        if (np.bitwise_and(que_line.flag, 0x800) == 0x800) or (que_line.reference_name not in chr_lst):
            continue
        FileID, ReadID, ReadLength, FrgID, SeqStart, SeqEnd = \
            [int(x.split(':')[1]) for x in que_line.query_name.split(';')]
        RefChrNid = chr_map[que_line.reference_name]
        RefStart = que_line.reference_start
        RefEnd = que_line.reference_end
        RefStrand = 1 - (que_line.is_reverse * 2)
        frg_info = np.array([ReadID, RefChrNid, RefStart, RefEnd, RefStrand, que_line.mapping_quality,
                             FileID, FrgID, SeqStart, SeqEnd, ReadLength, 1]).reshape([1, -1])

        # Check order of fragments
        if FrgID < FrgID_old:
            raise Exception('Order of fragments are lost.')
        FrgID_old = FrgID

        # Check if read has ended
        if ReadID == ReadID_old:
            frg_set = np.vstack([frg_set, frg_info])
            continue

        # Save the read
        saveRead(gz_fid, frg_set)
        frg_set = frg_info.copy()
        ReadID_old = ReadID
        n_processed += 1

    if frg_set.shape[0] != 0:  # Saving the last read after file has finished
        saveRead(gz_fid, frg_set)

# saving results in HDF5
print('Reading: {:s}'.format(gz_fname))
frg_pd = read_csv(gz_fname, delimiter='\t', compression='gzip')
frg_pd = frg_pd.iloc[np.lexsort([frg_pd['SeqStart'], frg_pd['ReadID']])]
h5_fname = './frag_files/frg_{:s}.hdf5'.format(vp_info['run_id'])
print('Writing: {:s}'.format(h5_fname))
h5_fid = h5py.File(h5_fname, 'w')
h5_fid.create_dataset('frg_np', data=frg_pd.values, compression='gzip', compression_opts=5,
                      chunks=(10, frg_pd.shape[1]))
h5_fid.create_dataset('frg_np_header_lst', data=list(frg_pd.columns.values))
h5_fid.create_dataset('chr_lst', data=chr_lst)
h5_fid.close()
