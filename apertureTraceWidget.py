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

#%%





#점원일 경우 sky를 지정하거나, gaussian+1차함수로 sky/data를 찾고



#점원이 아닐 경우는 sky를 지정해서 찾자 -> sky 지정만 자동화해줄수 있으면 될거같은데,
#자 생각해보자
# sky를 ypixel에 비례하는 값(상수 포함) 이라고 가정하자.
#일반적으로 sky는 양쪽 끝 중 하나에서 부터는 반드시 나올것이다.
#이미지를 10부분으로 나눠서 왼쪽 끝/ 오른쪽 끝에서부터 각 부분에 대해 linregress
#slope가 급격히 증가하거나 감소하는 부분을 이미지 시작과 끝으로 여긴다(어차피 피팅할거니까 실제 이미지랑 거리는
# (어차피 크게 멀지도 않겠지만) 그렇게 상관 없을것
#그리고 이미지 부분 말고 나머지에 linregression해서 sky를 찾고
#이를 이미지에 빼서 object를 찾고 이걸 가우시안(점상)으로 뽑아 더하거나 그냥 더해서 - 면상의 이미지같은 경우는 각 픽셀이 가우시안으로 퍼지는걸 고려해서 해야되나?
#extract
#이때 각 픽셀당 주어지는 값을 gaussian으로 그냥 더하거나 optimal extraction 해서 extraction 한다.

#과정
#aptrace -> peak 값과 sky fitting을 사용해 aperture의 중심과 gaussian을 찾고, aperture의 중심을 수식으로 표현한다(fitting한다)
#지금은 aperture의 중심에 따라 일관된 wavelength값을 y slice마다 사용하는데,
#근데 생각을 해보면, 만약 y 값에 따라 wavelength가 일정하지 않다면, 같은 x값에서 일정한 wavelength를 사용하는것도 틀리지 않을까 싶다
#그리고 aperture 중심에 대해 sky가 빠진 구간을 gauss fitting 해서 한 yslice의 wavelength에 대한 adu 총합을 구한다. -> extraction

#>> 표준화
#이걸 exp로 나누고(adu니까) 표준성의 spectrum과 비교해서 wavelength에 따른 sensitivity 그래프를 그린 다음에
#sensitivity 를 관측대상 스펙트럼에 곱해서
#결과를 도출한다.
#이상.

#1. Aptrace. 받는건 regFactor, Image(Cutted well)
# 내놓는건 aperture fit and/or aperture pixels

#점상에 대해 먼저
#Todo fitting 관련된거 다 scipy optimize curvefit으로 바꾸자
def gaussian(x, amplitude, mean, stddev ):
    return (amplitude *stddev*np.sqrt(2*np.pi) / np.exp(0))* 1/ (stddev * np.sqrt(2*np.pi )) * np.exp( - (x-mean)**2 / (2*stddev**2) )
def linearSky(x, slope, intercept):
    return x*slope + intercept

def skyImage(x, slope, intercept, amplitude, mean, stddev):
    return linearSky(x, slope, intercept) + gaussian(x, amplitude, mean, stddev)



C = fits.open("./Spectroscopy_Example/20181023/reduced/r_HD18247_60.0_1.fit")
# Find AP points
flux = C[0].data[40:110,300:800]
xFlux = np.arange(len(flux[0]))
yFlux = np.arange(len(flux))
#wavelength = fittedWavelength(xFlux, yFlux, regFactor, fitMethod)
skysubtractedFlux = []
apertureCoeffs = []
apertures = []
for i in range(len(flux[0])):
    fluxNow = flux[:,i]
    #plt.plot(fluxNow)
    #plt.axvline(peak_local_max(flux[:,i])[0][0])
    x = np.arange(len(fluxNow))
    fwhm = 2
    peakPix = peak_local_max(flux[:,i])[0][0]
    popt, pcov = curve_fit(skyImage, x, fluxNow, [0,0, fluxNow[peakPix], peakPix, 2*gaussian_fwhm_to_sigma])
    xx  = np.arange(0, len(fluxNow), 0.01)
    #plt.plot(xx, skyImage(xx, *popt))
    #plt.plot(xx, linearSky(xx, *popt[:2]) )
    #plt.plot(xx, gaussian(xx, *popt[2:]) )
    skysubtractedFlux.append(fluxNow - linearSky(x, *popt[:2]))
    apertureCoeffs.append(popt)
    apertures.append(popt[3])
    #plt.plot(fluxNow)
    #print(popt[3])
    plt.plot( i, popt[3], 'k+', ms = 3)
skysubtractedFlux = np.array(skysubtractedFlux)
apertures = np.array(apertures)
apertureCoeffs = np.array(apertureCoeffs)
skysubtractedFlux = skysubtractedFlux.T
plt.imshow(skysubtractedFlux)

fitMethod = 'chebyshev'
fitDeg = 3
sigmaATFitMask = 3
itersATFitMask = 2
itersATFit = 3
# APTrace with chebyshev(or polynomial)
# Extract itersATFit iteratively as mask the sigma( Maybe not useful after iter>3)


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
plt.plot(xFlux, apertures, 'k+')
plt.plot(xxFlux, fitVal)

#Extract abnormals by sigmaClip and refit
for iATFit in range(itersATFit):
    clip_mask = sigma_clip(apertures - valFunc(xFlux, fitted), sigma=sigmaATFitMask, maxiters=itersATFitMask).mask
    fitted = fitter(xFlux[~clip_mask], apertures[~clip_mask], fitDeg)
    xxFlux = np.arange(0, len(flux[0]), 0.001)
    xfitVal = valFunc(xxFlux, fitted)
    fitVal = valFunc(xFlux, fitted)
plt.plot(xFlux[~clip_mask], apertures[~clip_mask], 'k+')
plt.plot(xFlux[clip_mask], apertures[clip_mask], 'rx')
plt.plot(xxFlux, xfitVal)
# gauss Sum
wavelength = []
aduSum = []
for x, y in enumerate(fitVal):
    wavelength.append(fittedWavelength(x, y, regFactor, identificationMethod))
    aduSum.append(1 / np.sqrt(2 * np.pi) * apertureCoeffs[x][2])
wavelength = np.array(wavelength)
aduSum = np.array(aduSum)
plt.plot(wavelength,aduSum)

#real sum
wavelength = []
aduSum = []
for x, y in enumerate(fitVal):
    wavelength.append(fittedWavelength(x, y, regFactor, identificationMethod))
    aduSum.append(np.sum(skysubtractedFlux[:,x]))
wavelength = np.array(wavelength)
aduSum = np.array(aduSum)
plt.plot(wavelength,aduSum)


#M = fits.open("./Spectroscopy_Example/20181023/reduced/r_moon_0.1_1.fit")
#flux = M[0].data[40:110,300:800]