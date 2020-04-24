import os
import shutil

import numpy as np
from astropy.io import fits


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
    x1[mask] = xm

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
    data = [float(i) * 1000 for i in data[1:]]
    return np.array(data)


def appendToFITSdata(fileNameOut, columnName, dataFormat, dataArray):
    # Wow is this thing ugly...
    # All this workaround is necessary because the files opened with astropy fits are still open when closed(?!)
    # https://github.com/astropy/astropy/issues/7404
    with open(fileNameOut, 'rb') as f_in:
        fitsFile = fits.open(f_in, memmap=False)
        hduIndex = 1  # for molecfit

        origTable = fitsFile[hduIndex].data
        origCols = origTable.columns

        newCols = fits.ColDefs([fits.Column(name=columnName, format=dataFormat, array=dataArray)])
        hdu = fits.BinTableHDU.from_columns(origCols + newCols)
        fitsFile[hduIndex] = hdu

        with open('tmp.fits', 'wb') as f_out:
            fitsFile.writeto(f_out, overwrite=True)

    os.remove(fileNameOut)
    shutil.move('tmp.fits', fileNameOut)


def updateFITSdata(fileNameOut, columnName, dataFormat, dataArray):
    # Same story as with appendToFITSdata()
    with open(fileNameOut, 'rb') as f_in:
        fitsFile = fits.open(f_in, memmap=False)
        hduIndex = 1  # for molecfit

        fitsFile[hduIndex].data[columnName] = dataArray

        with open('tmp.fits', 'wb') as f_out:
            fitsFile.writeto(f_out, overwrite=True)

    os.remove(fileNameOut)
    shutil.move('tmp.fits', fileNameOut)
    fitsFile.close()


def updateFITSheader(fileNameOut, cpolys):
    with fits.open(fileNameOut, mode='update') as fitsFile:
        hduIndex = 1  # for molecfit
        header = fitsFile[hduIndex].header

        # Remove header old header entries
        if "NPOLY" in header:
            header = header[header.index("NPOLY")]
        fitsFile[hduIndex].header.append(("NPOLY", len(cpolys), "Number of subsequent polynomials used for norming"))

        for p_count, coefficients in enumerate(cpolys, 1):
            entry_name = f"P{p_count:03}"
            fitsFile[hduIndex].header[entry_name] = len(coefficients)-1
            fitsFile[hduIndex].header[entry_name].comment = "Polynomial ID with order"

            for c_count, coefficient in enumerate(coefficients):
                fitsFile[hduIndex].header[f"{entry_name}_{c_count}"] = coefficient
