import numpy as np


def get_chr_info(genome_str, property='chr_name'):
    chr_details = dict({
        'hg19': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
                'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,
                135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                155270560, 59373566, 16571
            ]
        }),
        'mm9': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172,
                129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430,
                166650296, 15902555, 16299
            ]
        }),
        'mm10': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110,
                130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566,
                171031299, 91744698, 16299,
            ]
        })
    })
    return chr_details[genome_str][property]


def get_re_info(re_name='DpnII', property='seq', genome_str=None):
    re_details = dict({
        'DpnII': dict({'seq': 'GATC'}),
        'Csp6I': dict({'seq': 'GTAC'}),
        'NlaIII': dict({'seq': 'CATG'}),
        'HindIII': dict({'seq': 'AAGCTT'})
    })

    if property == 'pos':
        re_fname = './renz_files/{:s}_{:s}.npz'.format(genome_str, re_name)
        return np.load(re_fname)['arr_0'][0]
    else:
        return re_details[re_name][property]


def hasOL(que_item, ref_lst, include_ref_left=False, include_ref_right=False, offset=0):
    if isinstance(que_item, list):
        que_item = np.array(que_item)
    que_dim = que_item.shape[0]
    [n_ref, ref_dim] = np.shape(ref_lst)
    result = np.ones(n_ref, dtype=bool)
    if que_dim != ref_dim or que_item.ndim != 1:
        raise ValueError('Query or reference are inconsistent')
    crd_ind = 0

    if que_dim == 4:  # Orientation
        result = que_item[3] == ref_lst[:, 3]
    if que_dim >= 3:  # Chromosome
        result = np.logical_and(result, que_item[0] == ref_lst[:, 0])
        crd_ind = 1
    if include_ref_left:
        OvlL = ref_lst[:, crd_ind] <= que_item[crd_ind+1] + offset
    else:
        OvlL = ref_lst[:, crd_ind] <  que_item[crd_ind+1] + offset
    if include_ref_right:
        OvlR = ref_lst[:, crd_ind+1] >= que_item[crd_ind] - offset
    else:
        OvlR = ref_lst[:, crd_ind+1] >  que_item[crd_ind] - offset
    result = np.logical_and(result, np.logical_and(OvlL, OvlR))
    return result

def findReferenceRestSites(refFile, ReSeq_lst, lineLen=50):
    # courtesy of R.Straver: https://github.com/rstraver
    import re

    ReEnz_ptr='|'.join(ReSeq_lst)
    restSitesDict = dict()

    with open(refFile,'r') as reference:
        offset = -lineLen
        matched_lst = []
        buffer_seq = ''
        for line in reference:
            if line[0] == '>':
                matched_lst.extend([[x.start()+offset, x.end()+offset] for x in (re.finditer(ReEnz_ptr, buffer_seq)) if x.start()])
                buffer_seq = 'N'*lineLen*2
                offset = -lineLen
                matched_lst = []
                restSitesDict[line[1:].rsplit()[0]] = matched_lst
            else:
                buffer_seq = buffer_seq[lineLen:] + line.rsplit()[0].upper()
                matched_lst.extend([[x.start()+offset, x.end()+offset] for x in (re.finditer(ReEnz_ptr, buffer_seq)) if x.start() < lineLen])
                offset += lineLen
        matched_lst.extend([[x.start()+offset, x.end()+offset] for x in (re.finditer(ReEnz_ptr, buffer_seq)) if x.start()])

    return restSitesDict


def load_track2bin(tsv_fname, genome, agr_func=np.mean, usecols=(0,1,2), bin_w=25000, delimiter='\t', compression='gzip'):
    from pandas import read_csv
    from _utilities import accum_array

    # Initialization
    chr_lst = get_chr_info(genome=genome, property='chr_name')
    chr_size = get_chr_info(genome=genome, property='chr_size')
    chr_map = dict(zip(chr_lst, range(1, len(chr_lst) + 1)))

    # Load track
    trk_pd = read_csv(tsv_fname, usecols=usecols, delimiter=delimiter, compression=compression)
    trk_pd.iloc[:, 0].replace(chr_map, inplace=True, regex=False)
    if str(trk_pd.iat[0,0].dtype) != 'int64':
        raise Exception('Unknown chromosome is detected.')
    trk_np = trk_pd.values
    del trk_pd

    if trk_np.shape[1] < 4:  # chr, coordinate, score
        trk_np = np.hstack([trk_np[:, 0:2], trk_np[:, 2:]])
    if np.any(trk_np[:,1] > trk_np[:,2]):
        raise Warning('Some rows have unsorted coordinates.')

    # Loop over chromosomes
    chr_idx = np.unique(trk_np[:, 0], return_inverse=True)[1]
    trk_grp = accum_array(chr_idx, trk_np)
    agr_lst = []
    for trk_cis in trk_grp:
        chr_nid = int(trk_cis[0, 0])
        bin_lst = np.arange(0, chr_size[chr_nid - 1], bin_w).reshape([-1, 1])
        n_bin = len(bin_lst) - 1
        bin_arr = np.hstack([bin_lst[:-1], bin_lst[1:]])
        bin_score = np.zeros([n_bin, 1])

        # Sorting
        s_idx = np.argsort(trk_cis[:, 1])
        trk_arr = trk_cis[s_idx, 1:3]
        trk_score = trk_cis[s_idx, 3]
        del trk_cis
        n_trk = trk_arr.shape[0]

        # Overlap check
        bi = tb = 0
        while bi < n_bin:
            while tb < n_trk and trk_arr[tb, 1] < bin_arr[bi, 0]:  # both frg indices are behind
                tb = tb + 1
            te = tb
            while te < n_trk and trk_arr[te, 0] < bin_arr[bi, 1]:  # frg end index is behind
                te = te + 1
            if tb < te:
                bin_score[bi] = agr_func(trk_score[tb:te])
            bi = bi + 1
        agr_lst.append(np.hstack([np.tile(chr_nid, [n_bin,1]), bin_arr, bin_score]))

    return np.vstack(agr_lst)
