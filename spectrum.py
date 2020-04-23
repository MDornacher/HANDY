#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
DESCRIPTION

definition of Spectrum class
"""
import os

import pandas as pd
from astropy.io import fits

def readSpectrum(filename,colWave=0,colFlux=1,skipRows=0):
    #print(filename)
    df = pd.read_csv(filename,header=None,delim_whitespace=True,comment='#',skiprows=skipRows)
    return Spectrum(wave=df[colWave].values,flux=df[colFlux].values,name=filename)

def saveSpectrum(filename,spectrum):
    print(spectrum.wave)
    if spectrum.wave is None:
        print("Spectrum must have wavelength table.")
    else:
        saveSpec=pd.DataFrame({'wave': spectrum.wave, 'flux': spectrum.flux})
        roundDict = {"wave": 4,\
                     "flux" : 6,\
                     }
        saveSpec = saveSpec.round(roundDict)
        with open(filename, 'w') as f:
            #f.write('# wave , flux\n')
            saveSpec.to_csv(f,columns=['wave','flux'],index=None,sep=' ',header=True)

def appendToFITSdata(fileName, columnName, dataFormat, dataArray):
    # TODO: the whole fits open/close thing is a bit wonky, probably better to use two with open
    fitsFile = fits.open(fileName)
    hduIndex = 1  # for molecfit
    origTable = fitsFile[hduIndex].data
    origCols = origTable.columns
    newCols = fits.ColDefs([fits.Column(name=columnName, format=dataFormat, array=dataArray)])
    hdu = fits.BinTableHDU.from_columns(origCols + newCols)

    fitsFile[hduIndex] = hdu
    outFileName = os.path.join(os.path.dirname(fileName), f'handy_{os.path.basename(fileName)}')
    with open(outFileName, 'wb') as f:
        fitsFile.writeto(f, overwrite=True)
    fitsFile.close()

def updateFITSdata(fileName, columnName, dataFormat, dataArray):
    # TODO: the whole fits open/close thing is a bit wonky, probably better to use two with open
    fitsFile = fits.open(fileName)
    hduIndex = 1  # for molecfit
    updatedCol = fits.ColDefs(fits.Column(name=columnName, format=dataFormat, array=dataArray))
    fitsFile[hduIndex].data[columnName] = updatedCol  # TODO: this probably needs more testing

    outFileName = os.path.join(os.path.dirname(fileName), f'handy_{os.path.basename(fileName)}')
    with open(outFileName, 'wb') as f:
        fitsFile.writeto(f, overwrite=True)
    fitsFile.close()

def updateFITSheader(fileName, cpolys):
    pass

class Spectrum:

    def __init__(self,name=None,\
                 wave=None,\
                 flux=None,\
                 ):
        self.name=name
        self.wave=wave
        self.flux=flux

    def __repr__(self):
        return "SPECTRUM: " + str(self.name)



################################################################################
### TESTS
################################################################################
def main():
    a=readSpectrum("exampleData/803432iuw.txt",skipRows=1)
    print(a)
    b=[]
    for i in range(10):
        b.append(Spectrum(name='lalal'))
    for i in range(10):
        print(id(b[i]))


if __name__ == '__main__':
	main()
