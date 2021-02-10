import sys
from PyQt5.QtWidgets import  *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from astropy.visualization import ZScaleInterval, ImageNormalize
from matplotlib.patches import Rectangle
from astropy.nddata import CCDData
from fitImageTableWidget import fitImageTableWidget, tableModel, fileOpener, openFitData, zimshow, znorm
from fitInfo import currentFileInfo
import os
from skimage.feature import peak_local_max
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import sigma_clip, gaussian_fwhm_to_sigma
from astropy.modeling.models import Gaussian1D, Chebyshev2D
from mpl_toolkits.mplot3d import Axes3D
from reIdentificatoinWidget import reIdentifier, fittedWavelength
from scipy import stats
from scipy.optimize import curve_fit
from scipy import signal

import math
#Todo 읽어볼것
# http://articles.adsabs.harvard.edu/pdf/1986PASP...98..609H

#받아야 하는거
#pixel to wavelength method <- regFactor

#preprocessed image file
#FWHM





B = fits.open("./Spectroscopy_Example/20181023/combine/comp_15_15.0.fit")
matchList = pd.DataFrame(dict(Pixel=[403, 374, 362, 265,
                                     245, 238, 223, 213,
                                     178, 169, 144, 131, 53],
                              Wavelength=[8780.6, 8495.4, 8377.6, 7438.9,
                                          7245.2, 7173.9, 7032.4, 6929.5,
                                          6599.0, 6507.0, 6266.5, 6143.1, 5400.6]))
flux = B[0].data
reId = reIdentifier(matchList, flux)
reId.doFit()
regFactor = reId.regFactor
identificationMethod = 'linear'
FWHM = 2

C = fits.open("./Spectroscopy_Example/20181023/reduced/r_HD18247_60.0_1.fit")
# Find
flux = C[0].data[40:110,300:800]


#Todo fitting 관련된거 다 scipy optimize curvefit으로 바꾸자
def gaussian(x, amplitude, mean, stddev ):
    return (amplitude *stddev*np.sqrt(2*np.pi) / np.exp(0))* 1/ (stddev * np.sqrt(2*np.pi )) * np.exp( - (x-mean)**2 / (2*stddev**2) )
def linearSky(x, slope, intercept):
    return x*slope + intercept

def skyImage(x, slope, intercept, amplitude, mean, stddev):
    return linearSky(x, slope, intercept) + gaussian(x, amplitude, mean, stddev)



C = fits.open("./Spectroscopy_Example/20181023/reduced/r_HD216604_60.0_9.fit")


flux = C[0].data[40:110,300:800]
xFlux = np.arange(len(flux[0]))
xxFlux =  np.arange(0, len(flux[0]), 0.01)
yFlux = np.arange(len(flux))
xyFlux = np.arange(0, len(flux), 0.01)
numYSlicer = 8

skyFit = linearSky


#wavelength = fittedWavelength(xFlux, yFlux, regFactor, fitMethod)
skysubtractedFlux = []
apertureCoeffs = []
apertures = []
aperturePointsFig = plt.figure()
imgAx = aperturePointsFig.add_subplot(411)
guessAx = aperturePointsFig.add_subplot(412)
fitAx = aperturePointsFig.add_subplot(413)
subImgAx = aperturePointsFig.add_subplot(414)
imgAx.imshow(flux)
for i in range(len(flux[0])):

    guessAx.clear()
    fitAx.clear()
    fluxNow = flux[:,i]
    #find initial guess for mean value in image gaussian mean
    peaks,_ = signal.find_peaks(fluxNow)
    prominences = signal.peak_prominences(fluxNow, peaks)[0]
    peakPixGuess = peaks[np.argmax(prominences)]
    imgAx.plot(i, peakPixGuess, 'r+', ms = 5, label = 'initial guess on gauss peak')
    #find initial guess for linear sky
    #extract gaussian and linear fit
    #peakPix 주변 FWHM 3배 만큼의 픽셀을 뺀 후에 linearfit
    yObj = np.arange(peakPixGuess - FWHM * 2, peakPixGuess + FWHM * 2)
    ySky = np.delete(yFlux, yObj)
    fluxSky = np.delete(fluxNow, yObj)

    guessAx.plot(yFlux, fluxNow , 'r-')

    interceptGuess = fluxSky[0]
    slopeGuess = (fluxSky[-1] - fluxSky[0]) / len(fluxSky)

    try:
        poptSky, pcov = curve_fit(skyFit, ySky, fluxSky, [slopeGuess, interceptGuess])
    except:
        poptSky = [0,0]
    slopeGuess = poptSky[0]
    interceptGuess = poptSky[1]
    guessAx.plot(xyFlux, skyFit(xyFlux, *poptSky), 'y--')

    #gaussGuess plot after extraction of sky
    gaussGuessFlux = fluxNow -skyFit(yFlux, *poptSky)

    try:
        poptGauss, pcov = curve_fit(gaussian, yFlux, gaussGuessFlux, [gaussGuessFlux[peakPixGuess], peakPixGuess, FWHM*gaussian_fwhm_to_sigma])
    except:
        poptGauss= [0,0,0]

    guessAx.plot(xyFlux, gaussian(xyFlux, *poptGauss), 'y--')
    guessAx.plot(xyFlux, skyImage(xyFlux, *poptSky, *poptGauss), 'b--')

    amplitudeGuess = poptGauss[0]
    meanGuess = poptGauss[1]
    stddevGuess = poptGauss[2]

    try:
        poptFin, pcov = curve_fit(skyImage, yFlux, fluxNow, [slopeGuess,interceptGuess, amplitudeGuess,meanGuess,stddevGuess ])
    except:
        poptFin= [0,0,0,0,0]

    fitAx.plot(yFlux, fluxNow, 'r-')
    fitAx.plot(xyFlux, skyImage(xyFlux, *poptFin), 'b--')

    '''
    #plt.plot(fluxNow)
    #plt.axvline(peak_local_max(flux[:,i])[0][0])
    x = np.arange(len(fluxNow))
    fwhm = 2

    fluxLen = len(fluxNow)
    slopes = []
    intercepts = []
    for fluxSlice, xSlice in zip(np.array_split(fluxNow, numYSlicer), np.array_split(yFlux, numYSlicer)):
        result = stats.linregress(xSlice, fluxSlice)
        plt.plot(xSlice, xSlice*result.slope+result.intercept)
        slopes.append(result.slope)
        intercepts.append(result.intercept)
    slopes = np.array(slopes)
    slopeMask = np.argmax(np.abs(slopes - np.median(slopes)))
    slopes = np.delete(slopes, slopeMask)
    intercepts = np.delete(intercepts, slopeMask)

    slope = np.mean(slopes)
    intercept = np.mean(intercepts)
    skysubtractedFluxNow = fluxNow - (slope*yFlux + intercept)

    plt.plot(skysubtractedFluxNow)
    peakPix = peak_local_max(skysubtractedFluxNow)[0][0]


    try:
        popt, pcov = curve_fit(gaussian, x, skysubtractedFluxNow, [skysubtractedFluxNow[peakPix], peakPix, 2*gaussian_fwhm_to_sigma])
    except:
        popt= [0,0,0]
    '''


    #plt.plot(fluxNow)
    #print(popt[3])
    imgAx.plot( i, poptFin[3], 'k+', ms = 3, label = 'final Aperture')
    skysubtractedFluxNow = fluxNow - (poptFin[0] * yFlux + poptFin[1])

    apertureCoeffs.append(poptFin)
    apertures.append(poptFin[3])
    skysubtractedFlux.append(skysubtractedFluxNow)


skysubtractedFlux = np.array(skysubtractedFlux)
apertures = np.array(apertures)
apertureCoeffs = np.array(apertureCoeffs)
skysubtractedFlux = skysubtractedFlux.T
subImgAx.imshow(skysubtractedFlux)


fitMethod = 'chebyshev'
fitDeg = 3
sigmaATFitMask = 2
itersATFitMask = 2
itersATFit = 3
# APTrace with chebyshev(or polynomial)
# Extract itersATFit iteratively as mask the sigma( Maybe not useful after iter>3)

apertureTraceFig = plt.figure()


firstAx = apertureTraceFig.add_subplot(311)


if (fitMethod=='chebyshev'):
    fitter = np.polynomial.chebyshev.chebfit
    valFunc = np.polynomial.chebyshev.chebval
elif (fitMethod == 'polynomial'):
    fitter = np.polynomial.polynomial.polyfit
    valFunc = np.polynomial.polynomial.polyval

#first
fitted = fitter(xFlux, apertures, fitDeg)
xxFlux = np.arange(0, len(flux[0]), 0.001)
fitVal = valFunc(xxFlux, fitted)
firstAx.plot(xFlux, apertures, 'k+')
firstAx.plot(xxFlux, fitVal, 'g--')
firstAx.imshow(skysubtractedFlux)

secondAx = apertureTraceFig.add_subplot(312)
#Extract abnormals by sigmaClip and refit
for iATFit in range(itersATFit):
    clip_mask = sigma_clip(apertures - valFunc(xFlux, fitted), sigma=sigmaATFitMask, maxiters=itersATFitMask).mask
    fitted = fitter(xFlux[~clip_mask], apertures[~clip_mask], fitDeg)
    xxFlux = np.arange(0, len(flux[0]), 0.001)
    xfitVal = valFunc(xxFlux, fitted)
    fitVal = valFunc(xFlux, fitted)

secondAx.plot(xFlux[~clip_mask], apertures[~clip_mask], 'k+')
secondAx.plot(xFlux[clip_mask], apertures[clip_mask], 'rx')
secondAx.plot(xxFlux, xfitVal, 'g--')
secondAx.imshow(skysubtractedFlux)

residualAx = apertureTraceFig.add_subplot(313)
residualAx.plot(xFlux[~clip_mask], apertures[~clip_mask] - fitVal[~clip_mask], 'k+')
residualAx.plot(xFlux[clip_mask], apertures[clip_mask] - fitVal[clip_mask], 'rx')
residualAx.axhline(0, color= 'g', linestyle = '--')
residualAx.set_ylim(-0.5,0.5)

wavelengthFig = plt.figure()
ax = wavelengthFig.add_subplot(111)
# gauss Sum
wavelength = []
aduSum = []
for x, y in enumerate(fitVal):
    wavelength.append(fittedWavelength(x, y, regFactor, identificationMethod))
    aduSum.append(1 / np.sqrt(2 * np.pi) * apertureCoeffs[x][2])
wavelength = np.array(wavelength)
aduSum = np.array(aduSum)
ax.plot(wavelength,aduSum)

#real sum
wavelength = []
aduSum = []
for x, y in enumerate(fitVal):
    wavelength.append(fittedWavelength(x, y, regFactor, identificationMethod))
    aduSum.append(np.sum(skysubtractedFlux[:,x]))
wavelength = np.array(wavelength)
aduSum = np.array(aduSum)
ax.plot(wavelength,aduSum)


#작은 x에 대해 값이 작아진다. 왜지? preprocessing 문제?



#M = fits.open("./Spectroscopy_Example/20181023/reduced/r_moon_0.1_1.fit")
#flux = M[0].data[40:110,300:800]