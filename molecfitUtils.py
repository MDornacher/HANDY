import numpy as np


# TODO this should be the place for all my funky new functions, to keep the HANDY core clean
def forward_fill_ifsame(x):
    # Get mask of non-zeros and then use it to forward-filled indices
    mask = x!=0
    idx = np.where(mask, np.arange(len(x)), 0)
    np.maximum.accumulate(idx, axis=0, out=idx)

    # Fill only if the previous and next one is the same
    x1 = x.copy()
    # Get non-zero elements
    xm = x1[mask]
    # Off the selected elements, we need to assign zeros to the previous places
    # that don't have their correspnding next ones different
    xm[:-1][xm[1:] != xm[:-1]] = 0
    # Assign the valid ones to x1. Invalid ones become zero.
    x1[mask] = x

    # Use idx for indexing to do the forward filling
    out = x1[idx]
    # For the invalid ones, keep the previous masked elements
    out[mask] = x[mask]
    return out


def regions2mask(wave, regions):
    waveTemplate = wave
    cid = np.zeros(len(wave))
    # Loop over regions in reversed order to get ascending corders simply by adding the updated mask
    for region in reversed(regions):
        for lower, upper in region:
            waveTemplate = np.ma.masked_where((waveTemplate >= lower) & (waveTemplate <= upper), waveTemplate)
        cid += waveTemplate.mask.astype(int)
    cmask = waveTemplate.mask.astype(int)
    return cmask, cid


def loadDIBs():
    with open("dibs") as f:
        data = f.read().splitlines()
    data = [float(i) * 1000 for i in data[2:]]
    return data