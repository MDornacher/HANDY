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


def appendToFITSdata(fileNameOut, columnName, dataFormat, dataArray, hduIndex):
    # Wow is this thing ugly...
    # All this workaround is necessary because the files opened with astropy fits are still open when closed(?!)
    # https://github.com/astropy/astropy/issues/7404
    with open(fileNameOut, 'rb') as f_in:
        fitsFile = fits.open(f_in, memmap=False)

        origTable = fitsFile[hduIndex].data
        origCols = origTable.columns

        newCols = fits.ColDefs([fits.Column(name=columnName, format=dataFormat, array=dataArray)])
        hdu = fits.BinTableHDU.from_columns(origCols + newCols)
        fitsFile[hduIndex] = hdu

        with open('tmp.fits', 'wb') as f_out:
            fitsFile.writeto(f_out, overwrite=True)

    os.remove(fileNameOut)
    shutil.move('tmp.fits', fileNameOut)


def updateFITSdata(fileNameOut, columnName, dataFormat, dataArray, hduIndex):
    # Same story as with appendToFITSdata()
    with open(fileNameOut, 'rb') as f_in:
        fitsFile = fits.open(f_in, memmap=False)

        fitsFile[hduIndex].data[columnName] = dataArray

        with open('tmp.fits', 'wb') as f_out:
            fitsFile.writeto(f_out, overwrite=True)

    os.remove(fileNameOut)
    shutil.move('tmp.fits', fileNameOut)
    fitsFile.close()


def newFITSdata(fileNameOut, columnName, dataFormat, dataArray):
    # Wow is this thing ugly...
    # All this workaround is necessary because the files opened with astropy fits are still open when closed(?!)
    # https://github.com/astropy/astropy/issues/7404
    with open(fileNameOut, 'rb') as f_in:
        fitsFile = fits.open(f_in, memmap=False)
        hduIndex = len(fitsFile)

        newCols = fits.ColDefs([fits.Column(name=columnName, format=dataFormat, array=dataArray)])
        hdu = fits.BinTableHDU.from_columns(newCols)
        fitsFile.append(hdu)

        with open('tmp.fits', 'wb') as f_out:
            fitsFile.writeto(f_out, overwrite=True)

    os.remove(fileNameOut)
    shutil.move('tmp.fits', fileNameOut)
    return hduIndex


def updateFITSheader(fileNameOut, cpolys, hduIndex):
    with fits.open(fileNameOut, mode='update') as fitsFile:
        header = fitsFile[hduIndex].header

        # Remove header old header entries
        if "NPOLY" in header:
            header = header[:header.index("NPOLY")]
        header.append(("NPOLY", len(cpolys), "Number of subsequent polynomials"))

        for p_count, coefficients in enumerate(cpolys, 1):
            entry_name = f"P{p_count:03}"
            header.append((entry_name, len(coefficients)-1, "Order of chebyshev polynomial"))

            for c_count, coefficient in enumerate(coefficients):
                header.append((f"{entry_name}_{c_count}", coefficient, f"{c_count} Chebyshev coefficient"))


def wrapFITSdata(fileNameOut, dataArrays, dataColumnNames, hduIndex):
    for columnName, dataArray in dataArrays.items():
        dataFormat = "D"  # TODO: ༼ つ ◕_◕ ༽つ GIVE FORMAT FROM ARRAY PLS
        if columnName in dataColumnNames:
            # Again, very ugly but prevents updating of wavelength and flux
            # which is only there because of non-molecfit FITS files
            if columnName not in ["lambda", "flux"]:
                updateFITSdata(fileNameOut, columnName, dataFormat, dataArray, hduIndex)
        elif hduIndex:
            appendToFITSdata(fileNameOut, columnName, dataFormat, dataArray, hduIndex)
        else:
            # This should only be necessary once
            hduIndex = newFITSdata(fileNameOut, columnName, dataFormat, dataArray)


def findHDUIndex(fileName):
    hduIndex = None
    tmpIndex = 0
    dataColumnNames = []

    # Search FITS data for Molecfit or HANDY keywords
    dataKeys = ["lambda", "flux", "cflux", "norm", "cont", "cmask", "pmask"]

    fitsFile = fits.open(fileName, memmap=False)
    for hdu in fitsFile:
        if hdu.data is None or not hasattr(hdu.data, 'names'):
            tmpIndex += 1
            continue
        if any(name in hdu.data.names for name in dataKeys):
            hduIndex = tmpIndex
            dataColumnNames = hdu.data.names
            break
        tmpIndex += 1
    return hduIndex, dataColumnNames


def refitContinuum(normAppLogic, cmask, pmask):
    cpolys = []
    for id_count, degree in enumerate(normAppLogic.continuumRegionsLogic.orders, 1):
        mask = (pmask != id_count) & (cmask != 1)
        coefficients = np.polynomial.chebyshev.chebfit(np.ma.masked_array(normAppLogic.spectrum.wave, mask=mask),
                                                       np.ma.masked_array(normAppLogic.spectrum.flux, mask=mask),
                                                       degree)
        cpolys.append(coefficients)
    return cpolys

def prepareOutputDataArrays(normAppLogic, cmask, pmask):
    dataArrays = {"lambda": normAppLogic.spectrum.wave,
                  "flux": normAppLogic.spectrum.flux,
                  "norm": normAppLogic.normedSpectrum.flux,
                  "cont": normAppLogic.continuum.flux,
                  "cmask": cmask,
                  "pmask": pmask,
                  # "spoints": np.zeros(len(self.normedSpectrum.flux))  # TODO: Special points for things like balmer jump
                  }  # TODO: Is there a better solution for this?
    return dataArrays