#!/usr/bin/env python

# Initialization
import numpy as np
import platform
if platform.system() == 'Linux':
    import matplotlib
    matplotlib.use('Agg')
from matplotlib import pyplot as plt, cm
from seq_tools import get_chr_info, hasOL
from _utilities import accum_array, showprogress
from mchc_tools import get_vp_info, get_mchc_data
from matplotlib.colors import LinearSegmentedColormap

run_list = [0] # range(6)
vp_info_lst = [get_vp_info(vpi) for vpi in run_list]
chr_lst = get_chr_info(vp_info_lst[0]['genome'], 'chr_name')[:-1]
chr_size = get_chr_info(vp_info_lst[0]['genome'], 'chr_size')[:-1]
run_id = ','.join([vp_info_lst[i]['run_id'] for i in range(len(vp_info_lst))])

# Bin list
bin_w = 1000000
MAX_FREQ = 50
hic_index = []
for chr_ind, chr_name in enumerate(chr_lst):
    bin_lst = np.arange(0, chr_size[chr_ind], bin_w).reshape([-1, 1])
    bin_arr = np.hstack([bin_lst[:-1], bin_lst[1:]])
    n_bin = bin_arr.shape[0]
    hic_index.append(np.hstack([np.tile(chr_ind+1, [n_bin, 1]), bin_arr]))
hic_index = np.vstack(hic_index)
n_bin = hic_index.shape[0]

# Load data
frg_pd = get_mchc_data(vp_info_lst, min_MQ=20)
frg_pd = frg_pd[['ReadID', 'Chr', 'RefStart', 'RefEnd', 'Strand']]

# Loop over circles
hic_mat = np.zeros([n_bin, n_bin], dtype=np.int16)
cir_uid = np.unique(frg_pd['ReadID'].values, return_inverse=True)[1]
cir_grp = accum_array(cir_uid, frg_pd.values)
n_cir = len(cir_grp)
print('Mapping [{:d}] circles on {:d} bins:'.format(n_cir, n_bin))
for cir_idx, cir in enumerate(cir_grp):
    if len(cir) == 0:
        continue

    index_lst = []
    for fi in range(cir.shape[0]):
        index_lst.append(np.where(hasOL(cir[fi, 1:4], hic_index))[0]) # sometimes, a fragment can overlap with multiple bins
    index_lst = np.hstack(index_lst)
    if len(index_lst) < 2:
        continue
    showprogress(cir_idx, n_cir, n_step=25)

    for bi in index_lst:
        hic_mat[bi, index_lst] = hic_mat[bi, index_lst] + 1
        hic_mat[index_lst, bi] = hic_mat[index_lst, bi] + 1
np.fill_diagonal(hic_mat, 0)
# hic_mat[hic_mat == 0] = np.nan

# Plotting
plt.figure(figsize=(31,30))
clr_map = [cm.hot(ci) for ci in np.linspace(1.0, 0.0, 10)]
cmap = LinearSegmentedColormap.from_list('test', clr_map, N=10)
plt.imshow(hic_mat, cmap=cmap, origin='upper', interpolation='none')
cbar_h = plt.colorbar()
cbar_h.ax.tick_params(labelsize=14)
plt.clim(0, MAX_FREQ)

# Add chromosome lines
boundary_lst = list(np.unique(hic_index[:,0], return_index=True)[1])
x_tick = [int(np.mean(boundary_lst[i:i+2])) for i in range(len(boundary_lst))]
x_tick_lbl = [chr_lst[i] for i in hic_index[boundary_lst, 0]-1]
plt.xticks(x_tick, x_tick_lbl, fontsize=12, rotation=25)
plt.yticks(x_tick, x_tick_lbl, fontsize=12, rotation=0)
for xi in boundary_lst + [n_bin]:
    plt.plot([0, n_bin], [xi, xi], color='#acacac', linewidth=0.5, alpha=0.5)
    plt.plot([xi, xi], [0, n_bin], color='#acacac', linewidth=0.5, alpha=0.5)
plt.title('{:s}\n'.format(run_id) +
          '#read={:d}, bin-w={:0.1f}kb'.format(n_cir, bin_w/1e3))
plt.xlim([0, n_bin])
plt.ylim([0, n_bin])
plt.savefig('./plots/05_Hi-C_Maps_{:s}_bw{:0.0f}kb_mf{:0.0f}.pdf'.format(run_id, bin_w / 1e3, MAX_FREQ), bbox_inches='tight', transparent=True)
