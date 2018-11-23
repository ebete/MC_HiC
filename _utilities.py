import numpy as np

def showprogress(iter, n_iter, n_step=10, output_format='{:1.0f}%, '):
    iter = iter + 1
    if ((iter % (n_iter / float(n_step))) - ((iter - 1) % (n_iter / float(n_step))) < 0) or (n_iter / float(n_step) <= 1):
        print(output_format.format(iter * 100 / n_iter)),
        if iter == n_iter:
            print

def accum_array(group_idx, arr, func=None, default_value=None, min_n_group=None, rebuild_index=False):
    """groups a by indices, and then applies func to each group in turn.
    e.g. func example: [func=lambda g: g] or [func=np.sum] or None for speed
    based on https://github.com/ml31415/numpy-groupies
    """

    if rebuild_index:
        group_idx = np.unique(group_idx.copy(), return_inverse=True)[1]
    if not min_n_group:
        min_n_group = np.max(group_idx) + 1

    order_group_idx = np.argsort(group_idx, kind='mergesort')
    counts = np.bincount(group_idx, minlength=min_n_group)

    if isinstance(arr, np.ndarray):
        groups = np.split(arr[order_group_idx], np.cumsum(counts)[:-1], axis=0)
    else:  # If arr is a Pandas DataFrame
        groups = np.split(arr.loc[order_group_idx,:], np.cumsum(counts)[:-1], axis=0)

    if func:
        ret = [default_value] * min_n_group
        for i, grp in enumerate(groups):
            if len(grp) > 0:
                ret[i] = func(grp)
        return ret
    else:
        return groups


def mapEx(list_idx, arr, func=None, default_value=None):
    ret = [default_value] * len(list_idx)
    for i, set_idx in enumerate(list_idx):
        if len(set_idx) > 0:
            ret[i] = func(arr[set_idx])
    return ret