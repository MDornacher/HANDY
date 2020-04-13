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

def appendToFITS(fileName, column_name, data_format, data_array):
    fits_file = fits.open(fileName)
    hduIndex = 1  # for molecfit
    orig_table = fits_file[hduIndex].data
    orig_cols = orig_table.columns
    new_cols = fits.ColDefs([fits.Column(name=column_name, format=data_format, array=data_array)])
    hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)

    fits_file[hduIndex] = hdu
    outFileName = os.path.join(os.path.dirname(fileName), f'handy_{os.path.basename(fileName)}')
    with open(outFileName, 'wb') as f:
        fits_file.writeto(f, overwrite=True)
    fits_file.close()

def updateFITS(fileName, column_name, data_format, data_array):
    fits_file = fits.open(fileName)
    hduIndex = 1  # for molecfit
    updated_col = fits.ColDefs(fits.Column(name=column_name, format=data_format, array=data_array))
    fits_file[hduIndex].data[column_name] = updated_col  # TODO: this probably needs more testing

    outFileName = os.path.join(os.path.dirname(fileName), f'handy_{os.path.basename(fileName)}')
    with open(outFileName, 'wb') as f:
        fits_file.writeto(f, overwrite=True)
    fits_file.close()

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
