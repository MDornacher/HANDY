#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import copy
from scipy.interpolate import Akima1DInterpolator
import matplotlib.pyplot as plt
from PyAstronomy import pyasl
from astropy.io import fits

import spectrum as sp
import regionLogic
import radialVelocity
import specInterface
import gridDefinitionsRead

import tkinter

"""
DESCRIPTION

"""

class normAppLogic:

    def __init__(self,):
        self.folderHANDY = os.path.dirname(os.path.abspath(__file__))
        gridDefinitionsFile = os.path.join(self.folderHANDY, "gridsDefinitions.yaml")

        self.continuumRegionsLogic = regionLogic.RegionLogic()
        self.radialVelocityEstimator = radialVelocity.RadialVelocity()
        self.specSynthesizer = specInterface.SynthesizeSpectrum()
        self.gridDefinitions = gridDefinitionsRead.gridDefinition(gridDefinitionsFile)

        self.spectrum = sp.Spectrum()
        self.theoreticalSpectrum = sp.Spectrum(wave=[],\
                                               flux=[])
        self.continuum = sp.Spectrum(wave=[],\
                                     flux=[])
        self.normedSpectrum = sp.Spectrum(wave=[],\
                                          flux=[])
        self.radialVelocity = 0.0
        self.oryginalWavelength = None

    def _ask_multiple_choice_question(self,prompt,options):
        # source: https://stackoverflow.com/questions/42581016/how-do-i-display-a-dialog-that-asks-the-user-multi-choice-question-using-tkinter
        pass
        #win = tkinter.Toplevel()
        #if prompt:
        #    tkinter.Label(win, text=prompt).pack()
        #v = tkinter.IntVar()
        #for i, option in enumerate(options):
        #    tkinter.Radiobutton(win, text=option, variable=v, value=i).pack(anchor="w")
        #tkinter.Button(text="Submit", command=win.destroy).pack()
        #win.mainloop()
        #if v.get() == 0:
        #    return None
        #return options[v.get()]


    def readSpectrum(self,fileName,colWave=0,colFlux=1,skipRows=0):
        """
        Reading spectrum in text or FITS format
        and update regions and points
        """
        if not ".fits" in fileName:
            self.spectrum = sp.readSpectrum(fileName,\
                                                  colWave=colWave,\
                                                  colFlux=colFlux,\
                                                  skipRows=skipRows)

        else:
            self.spectrum = sp.Spectrum()
            """ Check more at
            http://archive.eso.org/cms/eso-data/help/1dspectra.html
            https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/readFitsSpec.html
            """
            if "_tpl" in fileName:  # TODO: this tpl is actually a typo and should be tac
                fits_file = fits.open(fileName)

                # HARD CODED
                hduIndex = 1  # for molecfit
                waveKey = "lambda"
                fluxKey = "cflux"

                fitsData = fits_file[hduIndex].data
                self.spectrum.wave = fitsData[waveKey]
                self.spectrum.flux = fitsData[fluxKey]

                self.spectrum.wave = self.spectrum.wave * 1000

                # DYNAMIC GUI SELECTION
                # hduIndex = self._ask_multiple_choice_question("Select HDU", range(len(fits_file)))
                # fitsData = fits_file[hduIndex].data

                # waveKey = self._ask_multiple_choice_question("Select wave key", fitsData.names)
                # self.spectrum.wave = fitsData[waveKey]

                # fluxKey = self._ask_multiple_choice_question("Select flux key", fitsData.names)
                # self.spectrum.flux = fitsData[fluxKey]

                fits_file.close()
            else:
                self.spectrum.wave, self.spectrum.flux = pyasl.read1dFitsSpec(fileName)
                # self.spectrum.wave = self.spectrum.wave.byteswap().newbyteorder()
                self.spectrum.flux = self.spectrum.flux.byteswap().newbyteorder()  # TODO PyAstronomy bug
            self.spectrum.name = fileName
        self.radialVelocity = 0.0
        self.oryginalWavelength = copy.deepcopy(self.spectrum.wave)


    def saveSpectrum(self,fileName):
        sp.saveSpectrum(fileName,self.spectrum)
        print("INFO : %s saved!"%fileName)

    def readTheoreticalSpectrum(self,fileName,colWave=0,colFlux=1,skipRows=0):
        self.theoreticalSpectrum = sp.readSpectrum(fileName,\
                                                   colWave=colWave,\
                                                   colFlux=colFlux,\
                                                   skipRows=skipRows)

    def computeTheoreticalSpectrum(self,teff,logg,vmic,me,vsini,vmac,resolution):
        parameters = teff,logg,vmic,me,vsini,vmac,resolution
        # There was a bug - sometimes data from TKinter comes as string, so:
        parameters = [p if p is not str else float(p.replace(",","."))for p in parameters]
        try:
            self.theoreticalSpectrum = self.specSynthesizer.synthesizeSpectrum(parameters,minWave = 3500, maxWave = 7000)
        except:
            print("Spectrum out of grid or some interpolation bug...")

    def saveNormedSpectrum(self,fileName,correctForRadialVelocity):
        saveSpectrum = copy.deepcopy(self.normedSpectrum)
        # print(correctForRadialVelocity)
        if not correctForRadialVelocity:
            print("Modifying to oryginal wavelength.")
            saveSpectrum.wave = self.oryginalWavelength
        else:
            print("Saving corrected for radial velocity.")
        sp.saveSpectrum(fileName,saveSpectrum)
        print("INFO : %s saved!"%fileName)

    def saveTheoreticalSpectrum(self,fileName):
        sp.saveSpectrum(fileName,self.theoreticalSpectrum)
        print("INFO : %s saved!"%fileName)

    def saveToFITS(self,fileName):
        # Save Normed Spectrum, Continuum and Continuum Mask to Molecfit FITS file
        hduIndex = 1  # for molecfit
        with fits.open(fileName) as fits_file:
            data_column_names = fits_file[hduIndex].data.names

        data_arrays = {'norm': self.normedSpectrum.flux,
                       'cont': self.continuum.flux,
                       'cmask': np.zeros(self.spectrum.wave.size),
                       'corder': np.zeros(self.spectrum.wave.size),
                       }  # TODO: Do we need more? Is there a better solution for this?

        for column_name, data_array in data_arrays.items():
            data_format = 'D'  # TODO: ༼ つ ◕_◕ ༽つ GIVE FORMAT FROM ARRAY PLS
            # TODO: each loop overwrites output of previous save
            # TODO: the input fits file should be overwritten instead of creating a new one => handy_handy_handy_XXX.fits
            if column_name in data_column_names:
                sp.updateFITS(fileName, column_name, data_format, data_array)
            else:
                sp.appendToFITS(fileName, column_name, data_format, data_array)

    def plotSpectrum(self):
        if self.spectrum.wave is not None:
            plt.plot(self.spectrum.wave,self.spectrum.flux)
            plt.show()
        else:
            print("WARNING: normAppLogic.plotSpectrum first read spectrum")

    def getValuesForPlotSpecialPoints(self):
        absPoints = self.continuumRegionsLogic.getAbsolutePoints(self.spectrum)
        x,y = zip(*absPoints)
        return x,y


    def getContinuumRangesForPlot(self):
        contRegionsWaveAndFlux = [[[[],[]]]]
        if self.spectrum.wave is not None:
            contRegionsWaveAndFlux = self.continuumRegionsLogic.waveToSpectrumParts(self.spectrum)
        else:
            print("WARNING: normAppLogic.getContinuumRangesForPlot:\n"\
                 +"first load spectrum for norming")
        return contRegionsWaveAndFlux


    def normSpectrum(self):
        sep = 1
        w,f = self.fitFunctionToRegions(separation = sep)
        # Insert special points before interpolation
        absolutePoints = self.continuumRegionsLogic.getAbsolutePoints(self.spectrum)
        for x,y in absolutePoints:
            idx = np.searchsorted(w,x,side='left')
            w.insert(idx,x)
            f.insert(idx,y)
        # -----
        if len(w)>1:
            interp = Akima1DInterpolator(w, f)
            self.continuum.wave = self.spectrum.wave
            self.continuum.flux = interp(self.spectrum.wave,extrapolate=False)
            self.normedSpectrum.flux = self.spectrum.flux / self.continuum.flux
            np.nan_to_num(self.normedSpectrum.flux,copy=False)
        else:
            self.continuum.wave = []
            self.continuum.flux = []
            self.normedSpectrum.flux = []
            self.normedSpectrum.wave = []

    def ifOnNormedSpectrum(self,workOnNormedSpectrum):
        if workOnNormedSpectrum:
            self.continuumRegionsLogic.clearAll()
            self.continuum.wave = []
            self.continuum.flux = []
            self.normedSpectrum = copy.deepcopy(self.spectrum)
        else:
            self.normedSpectrum.wave = []
            self.normedSpectrum.flux = []

    def fitFunctionToRegions(self,separation=1):
        self.normedSpectrum.wave = copy.deepcopy(self.spectrum.wave)
        contRegionsWaveAndFlux = self.continuumRegionsLogic.waveToSpectrumParts(self.spectrum)
        wOut = []
        fOut = []
        for ord,region in zip(self.continuumRegionsLogic.orders,contRegionsWaveAndFlux):
            wReg = []
            fReg = []
            for w,f in region:
                wReg.extend(w)
                fReg.extend(f)
            wRegOut = np.linspace(wReg[0],wReg[-1],int(max((wReg[-1]-wReg[0])/separation,1)))
            try:
                fRegOut = self.fitFunction(wReg,fReg,wRegOut,1,ord)
            except Exception as e:
                fRegOut = []
                wRegOut = []
                print(e)
                print("WARNING: Unable to fit, try making region shorter")
            wOut.extend(wRegOut)
            fOut.extend(fRegOut)

        return wOut,fOut


    def fitFunction(self,xIn,yIn,xOut,fitType,degree):
        """
        fitType = 1 : Chebyshev
        fitType = 2 : Legendre
        """
        yOut=[]
        if fitType == 1:
            fit = np.polynomial.chebyshev.chebfit(xIn,yIn, degree)
            yOut = np.polynomial.chebyshev.chebval(xOut,fit)
        else: #elif fitType == 2:
            fit = np.polynomial.legendre.Legendre.fit(xIn,yIn, degree)
            yOut = np.polynomial.legendre.legval(xOut,fit)
        return yOut

    def applyRadialVelocity(self,radVel):
        self.spectrum.wave, self.spectrum.flux = self.radialVelocityEstimator.applyRadialVelocity(\
                                                         self.spectrum.wave,\
                                                         self.spectrum.flux,\
                                                         radVel)
        self.radialVelocity+=radVel

    def updateOrderOfActiveRegion(self,order):
        self.continuumRegionsLogic.setOrderOfActiveRegion(order)

################################################################################
### TESTS
################################################################################
def testLoadSaveSpectrum():
    nal=normAppLogic()

    nal.readSpectrum(os.path.join("exampleData", "803432iuw.txt"),skipRows=1)
    nal.saveSpectrum(os.path.join("exampleData", "saveTest.txt"))

    print(nal.spectrum)
    #nal.plotSpectrum()

def testGetContinuum():
    nal=normAppLogic()
    nal.readSpectrum(os.path.join("exampleData", "803432iuw.txt"),skipRows=1)

    nal.continuumRegionsLogic.addRegion([4850,4890])
    nal.continuumRegionsLogic.addRegion([5000,5100])
    nal.continuumRegionsLogic.addRegion([5500,5600])
    nal.continuumRegionsLogic.printRegions()
    indexRegions = nal.continuumRegionsLogic.waveToIndexRegions(nal.spectrum.wave)
    print(indexRegions)

    waveCont = []
    fluxCont = []
    for reg in indexRegions:
        for ran in reg:
            print(ran)
            waveCont.extend(nal.spectrum.wave[ran[0]:ran[1]])
            fluxCont.extend(nal.spectrum.flux[ran[0]:ran[1]])
    plt.plot(waveCont,fluxCont)
    plt.show()
def main():
    #testLoadSaveSpectrum()
    testGetContinuum()

if __name__ == '__main__':
	main()
