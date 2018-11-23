#!/usr/bin/env python3

import numpy as np
from pandas import read_csv
import h5py
import pandas
import os


def get_vp_info(run_id):
    vpi_lst = read_csv('vp_info.tsv', delimiter='\t')
    if isinstance(run_id, str):
        run_id = int(np.where(vpi_lst['run_id'] == run_id)[0])
    vp_info = vpi_lst.iloc[run_id].to_dict()
    vp_info['row_index'] = run_id
    if not vp_info:
        raise Exception('VP information could not be found.')
    return vp_info


def get_mchc_data(vp_info_lst, target_field='frg_np', data_path='./frag_files/',
                  min_MQ=20, only_hops=False, reindex_reads=True, max_rows=np.inf):
    MAX_N_CIR = 100000000
    out_pd = pandas.DataFrame()
    if not isinstance(vp_info_lst, list):
        vp_info_lst = [vp_info_lst]

    header_lst = []
    for vi, vp_info in enumerate(vp_info_lst):
        tsv_name = data_path + 'frg_' + vp_info['run_id'] + '.hdf5'
        print('Loading [{:d}: {:s}] data from: {:s}'.format(vp_info['row_index'], vp_info['run_id'], tsv_name))

        if not os.path.isfile(tsv_name):
            continue

        h5_fid = h5py.File(tsv_name, 'r')
        if np.isinf(max_rows):
            data_np = h5_fid[target_field].value
        else:
            print('Selecting only top [{:d}] rows in dataset'.format(max_rows))
            data_np = h5_fid[target_field][:max_rows]

        header_lst = list(h5_fid[target_field + '_header_lst'].value)
        h5_fid.close()
        part_pd = pandas.DataFrame(data_np, columns=header_lst)

        # Filtering fragments
        if min_MQ:
            part_pd = part_pd.loc[part_pd['MQ'] >= min_MQ]
        if only_hops:
            part_pd = part_pd.loc[part_pd['TrueHop'] > 0]

        # Adjust Read IDs
        assert np.max(part_pd['ReadID']) < MAX_N_CIR
        part_pd['ReadID'] = part_pd['ReadID'] + (vi + 1) * MAX_N_CIR

        # Append the part
        out_pd = out_pd.append(part_pd, ignore_index=True)
    out_pd = out_pd[header_lst]

    if reindex_reads:
        header_lst.append('ReadID_original')
        out_pd[header_lst[-1]] = out_pd['ReadID'].copy()
        out_pd['ReadID'] = np.unique(out_pd['ReadID'], return_inverse=True)[1] + 1
        print('Got [{:,d}] reads and [{:d}] fragments after re-indexing.'.format(
            np.max(out_pd['ReadID']), out_pd.shape[0]))

    return out_pd[header_lst]
