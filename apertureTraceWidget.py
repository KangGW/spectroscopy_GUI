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
from reIdentificationWidget import reIdentifier, fittedWavelength
from scipy import stats
from scipy.optimize import curve_fit
from scipy import signal
from fitInfo import editInfo
from editWidget import editWidget
from matplotlib.figure import Figure

import math
#Todo 읽어볼것
# http://articles.adsabs.harvard.edu/pdf/1986PASP...98..609H

#받아야 하는거
#pixel to wavelength method <- regFactor

#preprocessed image file
#FWHM



# 결과적으로 피팅이 완료된 상태의 이미지만 넣는다.

# image Figures

#이미지 Axes < -기본 이미지
#skysubtracted image Axes < sky substracted image
#sky image< sky images

# Plot Figures < - Before aperture fitting

# Plot 이미지 Axes <-  여기에 피팅 중심점(Aperture)이랑 가우스 1시그마(얇게)
# Todo 여기서 sky /object 결정할수 있는 interactive한 요소를 넣는다.
#yCut fitting plot
#residual

# Plot Figures  < - After aperture fitting

#fitted이미지 Axes <-  여기에 Aperture랑  Tracing 된 line, 사용되지 않는 point들
# Todo 여기서 사용할/ 사용하지 않을 Aperture 결정할 수 있는 interactive 요소를 넣는다.
#yCut fitting plot ?
#residual plot

#Plot Widget < - Before aperture fitting
#X Column and fitting 정보, 클릭해서 볼 수 있게.

#Table Widget <- After aperture fitting
#X Column and fitting 정보, 클릭해서 볼 수 있게/ 선택되거나 선택되지 않은 aperture들을 표시.




def gaussian(x, amplitude, mean, stddev ):
    return (amplitude *stddev*np.sqrt(2*np.pi) / np.exp(0))* 1/ (stddev * np.sqrt(2*np.pi )) * np.exp( - (x-mean)**2 / (2*stddev**2) )

def linearSky(x, slope, intercept):
    return x*slope + intercept

def skyImage(x, slope, intercept, amplitude, mean, stddev):
    return linearSky(x, slope, intercept) + gaussian(x, amplitude, mean, stddev)


class apertureTracer(QWidget):
    progressChangeSignal = pyqtSignal(float)

    def __init__(self, flux, skyFit, objFit, fig=None, waveFig = None, skyFitInfo = None):
        super().__init__()

        self.flux = flux
        self.skyFit = skyFit
        self.objFit = objFit
        self.FWHM = 4 #Todo add method for getting FHWM, gain, readoutnoise from image hdr
        self.cutExist = False
        self.norm = znorm(flux)
        self.apertureFitMethod = 'chebyshev'
        self.pickDistance = 3
        self.skyFitInfo = skyFitInfo
        print(self.norm)
        if fig is None: self.fig = plt.figure(figsize = (12,5))
        else: self.fig = fig
        if waveFig is None: self.waveFig = plt.figure()
        else: self.waveFig = waveFig

    def findAperturePoints(self):
        flux = self.flux
        fig = self.fig
        FWHM = self.FWHM
        skyFit = self.skyFit
        objFit = self.objFit

        xFlux = np.arange(len(flux[0]))
        xxFlux = np.arange(0, len(flux[0]), 0.01)
        yFlux = np.arange(len(flux))
        xyFlux = np.arange(0, len(flux), 0.01)
        skysubtractedFlux = []
        skyFlux = []
        apertureCoeffs = []
        apertures = []

        self.imgAx = fig.add_subplot(321)
        self.subImgAx = fig.add_subplot(323, sharex = self.imgAx )
        self.skyImgAx = fig.add_subplot(325, sharex = self.imgAx )
        self.apertureAx = fig.add_subplot(322, sharey=self.imgAx)
        self.fitAx = fig.add_subplot(324)
        self.residualAx = fig.add_subplot(326)


        zimshow(self.imgAx, self.flux, normalize=self.norm)
        zimshow(self.apertureAx, self.flux, normalize=self.norm)
        if self.skyFitInfo[4]=='Auto':

            for xPix in range(len(flux[0])):
                fluxNow = flux[:, xPix]
                # find initial guess for mean value in image gaussian mean
                # Peak을 prominence가 가장 큰걸 기준으로 찾아서 sky중 가장 큰 값이 peak이 아닐 경우를 대비.

                peaks, _ = signal.find_peaks(fluxNow)
                prominences = signal.peak_prominences(fluxNow, peaks)[0]
                peakPixGuess = peaks[np.argmax(prominences)]

                #imgAx.plot(xPix, peakPixGuess, 'r+', ms=5, label='initial guess on gauss peak')
                # find initial guess for linear sky
                # extract gaussian and linear fit
                # peakPix 주변 FWHM 2배 만큼의 픽셀을 뺀 후에 linearfit


                objMin = peakPixGuess - FWHM * 2 if peakPixGuess - FWHM * 2>0 else 0
                objMax = peakPixGuess + FWHM * 2 if peakPixGuess + FWHM * 2<len(fluxNow) else len(fluxNow)

                yObj = np.arange(objMin, objMax)
                ySky = np.delete(yFlux, yObj)
                fluxSky = np.delete(fluxNow, yObj)
                #guessAx.plot(yFlux, fluxNow, 'r-')
                interceptGuess = fluxSky[0]
                slopeGuess = (fluxSky[-1] - fluxSky[0]) / len(fluxSky)
                try:
                    poptSky, pcov = curve_fit(skyFit, ySky, fluxSky, [slopeGuess, interceptGuess])
                except:
                    poptSky = [0, 0]
                slopeGuess = poptSky[0]
                interceptGuess = poptSky[1]
                # 더 빨리? ->
                #SlopeGuess = np.median(fluxNow[-1:-11]) - np.median(fluxNow[0:10])
                #interceptGuess = np.median(fluxNow[0:10])

                #guessAx.plot(xyFlux, skyFit(xyFlux, *poptSky), 'y--')

                # gaussGuess plot after extraction of sky
                gaussGuessFlux = fluxNow - skyFit(yFlux, *poptSky)

                try:
                    poptGauss, pcov = curve_fit(gaussian, yFlux, gaussGuessFlux,
                                                [gaussGuessFlux[peakPixGuess], peakPixGuess, FWHM * gaussian_fwhm_to_sigma],
                                                bounds = ((0, peakPixGuess- FWHM,0),(gaussGuessFlux[peakPixGuess]*2,peakPixGuess+ FWHM ,FWHM)))
                except:
                    poptGauss = np.array([0, peakPixGuess, FWHM * gaussian_fwhm_to_sigma])

                #guessAx.plot(xyFlux, gaussian(xyFlux, *poptGauss), 'y--')
                #guessAx.plot(xyFlux, skyImage(xyFlux, *poptSky, *poptGauss), 'b--')

                amplitudeGuess = poptGauss[0]
                meanGuess = poptGauss[1]
                stddevGuess = poptGauss[2]
                # 더 빨리? ->
                #amplitudeGuess = gaussGuessFlux[peakPixGuess]
                #meanGuess = peakPixGuess
                #stddevGuess = FWHM * gaussian_fwhm_to_sigma
                try:
                    poptFin, pcov = curve_fit(skyImage, yFlux, fluxNow,
                                              [slopeGuess, interceptGuess, amplitudeGuess, meanGuess, stddevGuess],
                                              bounds = ((-np.inf, -np.inf, 0, peakPixGuess-FWHM, 0.001), (np.inf, np.inf, 2 * gaussGuessFlux[peakPixGuess], peakPixGuess+ FWHM, FWHM)))

                except:
                    poptSky, pcov = curve_fit(skyFit, yFlux, fluxNow, [slopeGuess, interceptGuess])
                    poptFin = np.array([poptSky[0], poptSky[1], 0, peakPixGuess, FWHM * gaussian_fwhm_to_sigma])

                #fitAx.plot(yFlux, fluxNow, 'r-')
                #fitAx.plot(xyFlux, skyImage(xyFlux, *poptFin), 'b--')
                # plt.plot(fluxNow)
                # print(popt[3])
                #imgAx.plot(xPix, poptFin[3], 'k+', ms=3, label='final Aperture')
                self.apertureAx.plot([xPix, xPix], [poptFin[3]-poptFin[4]*3, poptFin[3]+poptFin[4]*3], 'k', lw=1)
                self.apertureAx.plot(xPix,  poptFin[3], 'b,', ms=3, label='final Aperture points')
                skyFluxNow = (poptFin[0] * yFlux + poptFin[1])
                skysubtractedFluxNow = fluxNow - skyFluxNow

                apertureCoeffs.append(poptFin)
                apertures.append(poptFin[3])
                skysubtractedFlux.append(skysubtractedFluxNow)
                skyFlux.append(skyFluxNow)
                self.progressChangeSignal.emit(xPix/len(flux[0])*100)
        else:

            if self.skyFitInfo[4]=='polynomial':
                fitter = np.polynomial.polynomial.polyfit
                valFunc = np.polynomial.polynomial.polyval

            elif self.skyFitInfo[4] == 'chebyshev':
                fitter = np.polynomial.chebyshev.chebfit
                valFunc = np.polynomial.chebyshev.chebval
            for xPix in range(len(flux[0])):
                xCut1, xCut2, xCut3, xCut4 = self.skyFitInfo[:4].astype(int)
                orderFit = int(self.skyFitInfo[5])
                fluxNow = flux[:, xPix]
                xSky = np.append(np.arange(xCut1, xCut2), np.arange(xCut3, xCut4))
                fluxSky = fluxNow[xSky]
                coeff, fitfull = fitter(xSky, fluxSky, deg=orderFit, full=True)
                skyFluxNow = valFunc(np.arange(len(fluxNow)), coeff)
                skysubtractedFluxNow = fluxNow - skyFluxNow
                fluxObj = skysubtractedFluxNow[np.arange(xCut2, xCut3)]
                maxObj = np.max(fluxObj)
                peakPix = np.where(skysubtractedFluxNow == maxObj)[0]
                try:
                    poptGauss, pcov = curve_fit(gaussian, yFlux, skysubtractedFluxNow,
                                            [skysubtractedFluxNow[peakPix], peakPix, FWHM * gaussian_fwhm_to_sigma],
                                            bounds=((-np.inf, peakPix - FWHM, 0),
                                                    (np.inf, peakPix + FWHM, FWHM)))
                except:
                    poptGauss = np.array([0, peakPix, FWHM * gaussian_fwhm_to_sigma])

                aperture = poptGauss[1]
                apertureCoeff = np.append(coeff, poptGauss)
                apertureCoeffs.append(apertureCoeff)
                apertures.append(aperture)
                skysubtractedFlux.append(skysubtractedFluxNow)
                skyFlux.append(skyFluxNow)
                self.progressChangeSignal.emit(xPix / len(flux[0]) * 100)


        self.progressChangeSignal.emit(100)

        skysubtractedFlux = np.array(skysubtractedFlux)
        skyFlux = np.array(skyFlux)
        self.apertures = np.array(apertures).astype(float)
        self.apertureCoeffs = np.array(apertureCoeffs)
        self.skysubtractedFlux = skysubtractedFlux.T
        self.skyFlux = skyFlux.T
        zimshow(self.subImgAx, self.skysubtractedFlux, normalize=self.norm)
        zimshow(self.skyImgAx, self.skyFlux, normalize=self.norm)

        self.imgAx.set_xlim(0, len(flux[0]))
        self.imgAx.set_ylim(0, len(flux))
        self.apertureAx.set_xlim(0, len(flux[0]))
        self.apertureAx.set_ylim(0, len(flux))

    def apertureFittingCut(self, xCut):
        apertureAx = self.apertureAx
        flux = self.flux
        fig = self.fig
        FWHM = self.FWHM
        skyFit = self.skyFit
        objFit = self.objFit

        xFlux = np.arange(len(flux[0]))
        xxFlux = np.arange(0, len(flux[0]), 0.01)
        yFlux = np.arange(len(flux))
        xyFlux = np.arange(0, len(flux), 0.01)
        fluxNow = flux[:, xCut]
        popts = self.apertureCoeffs[xCut]

        if self.cutExist == True :  self.xCutter.remove()

        self.fitAx.clear()
        self.residualAx.clear()
        self.xCutter = apertureAx.axvline(xCut, color='r', picker=True, pickradius = self.pickDistance)
        self.fitAx.plot(yFlux, fluxNow, 'r-')
        if self.skyFitInfo[4]=='Auto':
            self.fitAx.plot(xyFlux, skyImage(xyFlux, *popts), 'b--')
            self.fitAx.plot(xyFlux, skyFit(xyFlux, *popts[:2]), 'y--')
            self.fitAx.plot(xyFlux, objFit(xyFlux, *popts[2:]), 'y--')

            self.residualAx.plot(yFlux, fluxNow-skyImage(yFlux, *popts), 'b--')
            self.residualAx.axhline(0, color='r', linestyle='-')

        else:
            if self.skyFitInfo[4] == 'polynomial':
                valFunc = np.polynomial.polynomial.polyval
            elif self.skyFitInfo[4] == 'chebyshev':
                valFunc = np.polynomial.chebyshev.chebval
                gaussCoeff = popts[-3:]
            coeff = popts[:-3]
            self.fitAx.plot(xyFlux, valFunc(xyFlux, coeff) + gaussian(xyFlux, *gaussCoeff), 'b--')
            self.fitAx.plot(xyFlux, valFunc(xyFlux, coeff), 'y--')
            self.fitAx.plot(xyFlux, gaussian(xyFlux, *gaussCoeff), 'y--')
            self.residualAx.plot(yFlux, fluxNow - valFunc(yFlux, coeff) - gaussian(yFlux, *gaussCoeff), 'b--')
            self.residualAx.axhline(0, color='r', linestyle='-')

        self.cutExist = True


    def apertureTrace(self, fitMethod = None):
        if fitMethod is not None: self.apertureFitMethod = fitMethod
        apertureAx = self.apertureAx
        fitAx = self.fitAx
        residualAx = self.residualAx
        flux = self.flux
        apertureAx.clear()
        fitAx.clear()
        residualAx.clear()

        xFlux = np.arange(len(flux[0]))
        xxFlux = np.arange(0, len(flux[0]), 0.01)
        yFlux = np.arange(len(flux))
        xyFlux = np.arange(0, len(flux), 0.01)
        skysubtractedFlux = self.skysubtractedFlux
        apertures = self.apertures
        apertureCoeffs = self.apertureCoeffs
        fitMethod = self.apertureFitMethod
        norm = znorm(skysubtractedFlux)
        fitDeg = 3
        sigmaATFitMask = 2
        itersATFitMask = 2
        itersATFit = 3
        # APTrace with chebyshev(or polynomial)
        # Extract itersATFit iteratively as mask the sigma( Maybe not useful after iter>3

        if (fitMethod == 'chebyshev'):
            fitter = np.polynomial.chebyshev.chebfit
            valFunc = np.polynomial.chebyshev.chebval
        elif (fitMethod == 'polynomial'):
            fitter = np.polynomial.polynomial.polyfit
            valFunc = np.polynomial.polynomial.polyval

        # first
        fitted = fitter(xFlux, apertures, fitDeg)
        xxFlux = np.arange(0, len(flux[0]), 0.001)
        fitVal = valFunc(xxFlux, fitted)

        # Extract abnormals by sigmaClip and refit
        for iATFit in range(itersATFit):
            clip_mask = sigma_clip(apertures - valFunc(xFlux, fitted), sigma=sigmaATFitMask,
                                   maxiters=itersATFitMask).mask
            fitted = fitter(xFlux[~clip_mask], apertures[~clip_mask], fitDeg)
            xxFlux = np.arange(0, len(flux[0]), 0.001)
            xfitVal = valFunc(xxFlux, fitted)
            fitVal = valFunc(xFlux, fitted)

        self.fitVal = fitVal

        zimshow(apertureAx, skysubtractedFlux, normalize=norm)

        apertureAx.plot(xFlux[~clip_mask], apertures[~clip_mask], 'k,')
        apertureAx.plot(xFlux[clip_mask], apertures[clip_mask], 'rx')

        zimshow(fitAx, skysubtractedFlux, normalize=norm)
        fitAx.plot(xFlux[~clip_mask], apertures[~clip_mask], 'k+', ms=5)
        fitAx.plot(xFlux[clip_mask], apertures[clip_mask], 'rx')
        fitAx.plot(xxFlux, xfitVal, 'b--')


        residualAx.plot(xFlux[~clip_mask], apertures[~clip_mask] - fitVal[~clip_mask], 'k+', ms=5)
        residualAx.plot(xFlux[clip_mask], apertures[clip_mask] - fitVal[clip_mask], 'rx')
        residualAx.axhline(0, color='b', linestyle='--')
        residualAx.set_ylim(-0.5, 0.5)

        apertureAx.set_xlim(0, len(flux[0]))
        fitAx.set_xlim(0, len(flux[0]))
        residualAx.set_xlim(0, len(flux[0]))

    def apertureExtract(self, regFactor, identificationMethod, sumMethod):
        fitVal = self.fitVal
        apertureCoeffs = self.apertureCoeffs
        skysubtractedFlux = self.skysubtractedFlux

        wavelengthFig = self.waveFig
        ax = wavelengthFig.add_subplot(111)

        # gauss Sum

        if sumMethod == 'gauss':
            wavelength = []
            aduSum = []
            for x, y in enumerate(fitVal):
                wavelength.append(fittedWavelength(x, y, regFactor, identificationMethod))
                aduSum.append(1 / np.sqrt(2 * np.pi) * apertureCoeffs[x][2]*apertureCoeffs[x][4])
            wavelength = np.array(wavelength)
            aduSum = np.array(aduSum)
            ax.plot(wavelength, aduSum)

        elif sumMethod == 'real':
            wavelength = []
            aduSum = []
            for x, y in enumerate(fitVal):
                wavelength.append(fittedWavelength(x, y, regFactor, identificationMethod))
                aduSum.append(np.sum(skysubtractedFlux[:, x])) # <3 sigma along the point
            wavelength = np.array(wavelength)
            aduSum = np.array(aduSum)
            ax.plot(wavelength, aduSum)

        #Todo unit caluation
        spec = np.vstack((wavelength,aduSum))

        return spec

#우선 이미지를 받아와서 다시 크롭한다

#Todo need to implement method that reject some point
class apertureSkyPickerWidget(QWidget):

    '''
    Sky pick widget
    select sky with linear fitting
    and
    with select sky by bar and fit with selected(chebychev or linear, #Todo add any more method)
    '''
    '''
    Algorithm
    First, find sky-estimation with linearsky fitting and set sky region 1 and 2 (red and blue, each region is with 2 bars, It moves, interactively in graph region
    Second, choose sky with Change of sky region, interactively. and fit with fit choose button
    Click Finish button to go to next step
    
    [    fitting canvas      ] [fit choose 
    [     image canvas       ]   buttons  ]
    [           finish button             ]
    
    '''

    skyfitSignal = pyqtSignal(object)
    """
    This signal is emitted when the widget finishes skyfit.
    It emits 6-length np.array with 4-length sky region, small to big, int and skyfit method str, and fitOrder int,
    """

    def __init__(self, currentFileLocation = ''):
        super().__init__()
        self.filename = currentFileLocation
        self.isPressed = False
        self.isXcutterPicked=False
        self.skyPlot = None
        self.objPlot = None
        _,self.data = openFitData(self.filename)
        self.pickDistance = 3
        self.findInitsky()
        self.initUI()

    def initUI(self):
        self.layout = QGridLayout()

        self.fittingfig = plt.Figure()
        self.fig = plt.Figure()
        self.fittingcanvas = FigureCanvas(self.fittingfig)
        self.canvas = FigureCanvas(self.fig)
        self.fittingAx = self.fittingfig.add_subplot(111)
        self.ax = self.fig.add_subplot(111)
        self.fittingAx.plot(np.arange(len(self.skyImageFlux)),self.skyImageFlux)

        xCut1, xCut2, xCut3, xCut4 = self.skyFitInfo[:4].astype(int)

        self.xCutter1 = self.fittingAx.axvline(xCut1, color='r', picker=True, pickradius=self.pickDistance)
        self.xCutter2 = self.fittingAx.axvline(xCut2, color='r', picker=True, pickradius=self.pickDistance)
        self.xCutter3 = self.fittingAx.axvline(xCut3, color='b', picker=True, pickradius=self.pickDistance)
        self.xCutter4 = self.fittingAx.axvline(xCut4, color='b', picker=True, pickradius=self.pickDistance)

        self.yCutter1 = self.ax.axhline(xCut1, color='r', lw = 0.7, picker=True, pickradius=self.pickDistance)
        self.yCutter2 = self.ax.axhline(xCut2, color='r', lw = 0.7, picker=True, pickradius=self.pickDistance)
        self.yCutter3 = self.ax.axhline(xCut3, color='b', lw = 0.7, picker=True, pickradius=self.pickDistance)
        self.yCutter4 = self.ax.axhline(xCut4, color='b', lw = 0.7, picker=True, pickradius=self.pickDistance)

        self.fittingcanvas.mpl_connect("pick_event", self.onPickAtfittingCanvas)
        self.fittingcanvas.mpl_connect("motion_notify_event", self.onMoveAtfittingCanvas)
        self.fittingcanvas.mpl_connect("button_release_event", self.onReleaseAtfittingCanvas)
        self.polyBtn = QPushButton('&polynomialFitting')
        self.chevBtn = QPushButton('&chebyshevFitting')
        self.finBtn = QPushButton('&Finish')
        self.polyBtn.clicked.connect(self.onPolyBtnClicked)
        self.chevBtn.clicked.connect(self.onChevBtnClicked)
        self.finBtn.clicked.connect(self.onFinish)
        self.argText = QLineEdit(self)
        self.argText.setPlaceholderText("Order (default = 1)")

        self.layout.addWidget(self.fittingcanvas, 0,0)
        self.fitBtn = QWidget(self)
        self.fitBtnLayout = QHBoxLayout()
        self.fitBtnLayout.addWidget(self.polyBtn)
        self.fitBtnLayout.addWidget(self.chevBtn)
        self.fitBtnLayout.addWidget(self.argText)
        self.fitBtnLayout.addWidget(self.finBtn)
        self.fitBtn.setLayout(self.fitBtnLayout)
        self.layout.addWidget(self.fitBtn, 1, 0)
        self.layout.addWidget(self.canvas, 2,0)
        zimshow(self.ax, self.data)
        self.fittingcanvas.draw()
        self.canvas.draw()

        self.setLayout(self.layout)

    def onPolyBtnClicked(self):
        self.skyFitInfo[4]='polynomial'
        self.fitSky()

    def onChevBtnClicked(self):
        self.skyFitInfo[4]='chebyshev'
        self.fitSky()

    def fitSky(self):
        if self.skyPlot is not None:
            self.skyPlot.pop(0).remove()
        if self.objPlot is not None:
            self.objPlot.pop(0).remove()
        if self.argText.text() == '':
            orderFit = 1
        else:
            orderFit = int(self.argText.text())
        self.skyFitInfo[5] = orderFit

        xCut1, xCut2, xCut3, xCut4 = self.skyFitInfo[:4].astype(int)


        self.xSky = np.append(np.arange(xCut1,xCut2), np.arange(xCut3,xCut4))

        self.fluxSky = self.skyImageFlux[self.xSky]

        if self.skyFitInfo[4]=='polynomial':
            fitter = np.polynomial.polynomial.polyfit
            valFunc = np.polynomial.polynomial.polyval

        elif self.skyFitInfo[4]=='chebyshev':
            fitter = np.polynomial.chebyshev.chebfit
            valFunc = np.polynomial.chebyshev.chebval

        else:
            print('choose proper skyfit method')

        coeff, fitfull = fitter(self.xSky, self.fluxSky, deg=orderFit, full=True)

        sky = valFunc(np.arange(len(self.skyImageFlux)),coeff)
        skysub = self.skyImageFlux- sky

        self.skyPlot = self.fittingAx.plot(sky,'g-')
        self.objPlot = self.fittingAx.plot(skysub,'g-')
        self.fittingAx.set_ylim(0- np.max(self.skyImageFlux)*0.2,np.max(self.skyImageFlux)*1.5)
        self.fittingcanvas.draw()

    def findInitsky(self):
        flux = self.data
        FWHM = 4
        skyFit = linearSky
        yFlux = np.arange(len(flux))
        fluxNow = np.median(flux[:, len(flux[0])//2-len(flux[0])//5:len(flux[0])//2+len(flux[0])//5], axis=1)
        self.skyImageFlux = fluxNow
        peaks, _ = signal.find_peaks(fluxNow)
        prominences = signal.peak_prominences(fluxNow, peaks)[0]
        peakPixGuess = peaks[np.argmax(prominences)]
        objMin = peakPixGuess - FWHM * 2 if peakPixGuess - FWHM * 2 > 0 else 0
        objMax = peakPixGuess + FWHM * 2 if peakPixGuess + FWHM * 2 < len(fluxNow) else len(fluxNow)

        yObj = np.arange(objMin, objMax)
        ySky = np.delete(yFlux, yObj)
        fluxSky = np.delete(fluxNow, yObj)

        interceptGuess = fluxSky[0]
        slopeGuess = (fluxSky[-1] - fluxSky[0]) / len(fluxSky)
        try:
            poptSky, pcov = curve_fit(skyFit, ySky, fluxSky, [slopeGuess, interceptGuess])
        except:
            poptSky = [0, 0]
        slopeGuess = poptSky[0]
        interceptGuess = poptSky[1]
        gaussGuessFlux = fluxNow - skyFit(yFlux, *poptSky)

        try:
            poptGauss, pcov = curve_fit(gaussian, yFlux, gaussGuessFlux,
                                        [gaussGuessFlux[peakPixGuess], peakPixGuess, FWHM * gaussian_fwhm_to_sigma],
                                        bounds=((0, peakPixGuess - FWHM, 0),
                                                (gaussGuessFlux[peakPixGuess] * 2, peakPixGuess + FWHM, FWHM)))
        except:
            poptGauss = np.array([0, peakPixGuess, FWHM * gaussian_fwhm_to_sigma])



        amplitudeGuess = poptGauss[0]
        meanGuess = poptGauss[1]
        stddevGuess = poptGauss[2]

        try:
            poptFin, pcov = curve_fit(skyImage, yFlux, fluxNow,
                                      [slopeGuess, interceptGuess, amplitudeGuess, meanGuess, stddevGuess],
                                      bounds=((-np.inf, -np.inf, 0, peakPixGuess - FWHM, 0.001), (
                                      np.inf, np.inf, 2 * gaussGuessFlux[peakPixGuess], peakPixGuess + FWHM, FWHM)))

        except:
            poptSky, pcov = curve_fit(skyFit, yFlux, fluxNow, [slopeGuess, interceptGuess])
            poptFin = np.array([poptSky[0], poptSky[1], 0, peakPixGuess, FWHM * gaussian_fwhm_to_sigma])

        meanObject = poptFin[3]
        sigmaObject = poptFin[4]
        skyLength = len(flux)/10
        objectSigmaRange = 4
        #skyFit Info from small pixel to big pixel
        #estimate sky as outbound of object~1/10 pixel_image
        self.skyFitInfo = np.array( [int (meanObject-objectSigmaRange*sigmaObject - skyLength) , int (meanObject-objectSigmaRange*sigmaObject), int (meanObject+objectSigmaRange*sigmaObject) , int (meanObject+objectSigmaRange*sigmaObject + skyLength),   'polynomial', 1])


    def setFileName(self, fileName):
        self.filename = fileName
        _,self.data = openFitData(self.filename)
        zimshow(self.ax, self.data)
        self.canvas.draw()


    def onPickAtfittingCanvas(self, event):
        if self.isXcutterPicked: return
        if event.artist in (self.xCutter1,self.xCutter2,self.xCutter3,self.xCutter4 ):
            self.isXcutterPicked = True
            self.numXcutterPicked = [self.xCutter1,self.xCutter2,self.xCutter3,self.xCutter4].index(event.artist)
            self.pickedXcutter = event.artist


    def onMoveAtfittingCanvas(self, event):
        if not event.inaxes: return
        if event.inaxes != self.fittingAx: return
        if self.isXcutterPicked:
            self.pickedXcutter.set_xdata(event.xdata)
            [self.yCutter1,self.yCutter2,self.yCutter3,self.yCutter4][self.numXcutterPicked].set_ydata(event.xdata)
            self.skyFitInfo[self.numXcutterPicked] = int(event.xdata)
        self.canvas.draw()
        self.fittingcanvas.draw()

    def onReleaseAtfittingCanvas(self, event):
        if not event.inaxes: return
        if event.inaxes != self.fittingAx: return
        if not self.isXcutterPicked: return
        self.isXcutterPicked = False
        self.skyFitInfo = np.append(np.sort(self.skyFitInfo[:4].astype(int)), self.skyFitInfo[4:])
        print(self.skyFitInfo)


    def onFinish(self):
        self.skyfitSignal.emit(self.skyFitInfo)
        self.close()



class apertureTraceWidget(QMainWindow):
    def __init__(self, fileName, regFactor = None, identificationMethod= None, savePath = None, fileList = None):
        super().__init__()
        self.FWHM = 2 #Todo add method for getting FHWM, gain, readoutnoise from image hdr
        self.fileName = fileName
        self.fileList = fileList
        self.hdr, self.flux = openFitData(self.fileName)
        self.skyFit = linearSky
        self.objFit = gaussian
        self.apertureTracerOpen = False
        self.isXcutterPicked = False
        self.editInfo = editInfo(y0=0, y1=len(self.flux), x0=0, x1=len(self.flux[0]), filename=self.fileName)
        self.norm = znorm(self.flux)
        self.regFactor = regFactor
        self.identificationMethod = identificationMethod
        self.savePath = savePath
        self.skyFitInfo = np.array([0,1,len(self.flux)-2,len(self.flux)-1,'linear',1])
        self.initUI()



    def initUI(self):
        self.layout = QGridLayout()
        self.mainWidget = QWidget()
        self.mainWidget.setLayout(self.layout)
        self.setCentralWidget(self.mainWidget)


        self.imageCanvas = FigureCanvas(Figure()) # contains first image and add points from aperturetrace
        self.mainImageAx = self.imageCanvas.figure.add_subplot(111)


        self.fittingCanvas = FigureCanvas(Figure()) # contains fitting images

        self.apertureExtractWidget = QWidget()
        self.spectrumCanvas = FigureCanvas(Figure())
        self.spectrumSaveFitBtn = QPushButton('Save &Fit')
        self.spectrumSaveFitBtn.clicked.connect(self.onSaveSpectrumInFit)

        self.apertureExtractLayout = QGridLayout()
        self.apertureExtractLayout.addWidget(self.spectrumCanvas, 0, 0, 1, -1)
        self.apertureExtractLayout.addWidget(self.spectrumSaveFitBtn, 1, 0)

        self.apertureExtractWidget.setLayout(self.apertureExtractLayout)




        self.cropBtn = QPushButton('&Crop')
        self.cropWidget = editWidget(self.fileName)
        self.cropBtn.clicked.connect(self.onCrop)
        self.cropWidget.cropCheckWidget.imageCropSignal.connect(self.imageCrop)

        self.skyPickerWidget = apertureSkyPickerWidget(self.fileName)
        self.skyPickerWidget.skyfitSignal.connect(self.onSkyfitSignal)
        self.apertureTracerBtn = QPushButton('&SkySelection')
        self.apertureTracerBtn.clicked.connect(self.onSkySelection)

        self.aperturePointBtn = QPushButton('Aperture&Points')
        self.aperturePointBtn.clicked.connect(self.onFindAperturePoints)

        self.apertureTraceBtn = QPushButton('Aperture&Trace')
        self.apertureTraceBtn.clicked.connect(self.onApertureTrace)

        self.apertureExtractBtn = QPushButton('Aperture&Extract')
        self.apertureExtractBtn.clicked.connect(self.onApertureExtract)


        self.progressBar = QProgressBar(self)



        self.layout.addWidget(self.imageCanvas, 0,0, 1, -1)
        self.layout.addWidget(self.cropBtn, 1,0, 1, -1)
        self.layout.addWidget(self.fittingCanvas, 2,0, 1, -1)
        self.layout.addWidget(self.apertureTracerBtn, 3,0)
        self.layout.addWidget(self.aperturePointBtn, 3,1)
        self.layout.addWidget(self.apertureTraceBtn, 3,2)
        self.layout.addWidget(self.apertureExtractBtn, 3, 3)

        self.layout.addWidget(self.progressBar,4,0,1,-1)

        self.showImage()


        self.fittingCanvas.mpl_connect("pick_event", self.onPickAtfittingCanvas)
        self.fittingCanvas.mpl_connect("motion_notify_event", self.onMoveAtfittingCanvas)
        self.fittingCanvas.mpl_connect("button_release_event", self.onReleaseAtfittingCanvas)

    @pyqtSlot(object)
    def onSkyfitSignal(self, skyfitInfo):
        self.skyFitInfo = skyfitInfo
        if self.apertureTracerOpen:
            self.apertureTracer.skyFitInfo = skyfitInfo


    def onCrop(self):
        self.cropWidget.show()
        self.cropWidget.raise_()

    @pyqtSlot(editInfo)
    def imageCrop(self, crop):
        self.cropWidget.close()
        self.editInfo = crop
        self.showImage()

    def showImage(self):
        data = self.flux[self.editInfo.y0:self.editInfo.y1,
               self.editInfo.x0:self.editInfo.x1]
        zimshow(self.mainImageAx, data, normalize=self.norm)
        self.mainImageAx.figure.canvas.draw()


    def onSkySelection(self):
        self.apertureTracerOpen = True
        self.skyPickerWidget.show()
        self.skyPickerWidget.raise_()
        self.data  = self.flux[self.editInfo.y0:self.editInfo.y1,
                     self.editInfo.x0:self.editInfo.x1]
        self.apertureTracer = apertureTracer(self.data, self.skyFit, self.objFit, fig=self.fittingCanvas.figure,
                                             waveFig = self.spectrumCanvas.figure, skyFitInfo = self.skyFitInfo)
        self.apertureTracer.progressChangeSignal.connect(self.onProgressChanged)


        self.fittingCanvas.draw()


    def onFindAperturePoints(self):
        if not self.apertureTracerOpen : return
        self.apertureTracer.findAperturePoints()
        self.apertureTracer.apertureFittingCut(int(len(self.data[0])/2))
        self.fittingCanvas.draw()

    def onApertureTrace(self):
        if not self.apertureTracerOpen: return
        self.apertureTracer.apertureTrace()
        self.fittingCanvas.draw()

    def onPickAtfittingCanvas(self, event):
        if self.isXcutterPicked: return
        if event.artist == self.apertureTracer.xCutter:
            self.isXcutterPicked = True

    def onMoveAtfittingCanvas(self, event):
        if not event.inaxes: return
        if event.inaxes != self.apertureTracer.apertureAx: return
        if self.isXcutterPicked:
            self.apertureTracer.xCutter.remove()
            self.apertureTracer.xCutter = self.apertureTracer.apertureAx.axvline(int(event.xdata), color='red', picker=True, pickradius=self.apertureTracer.pickDistance)
            self.apertureTracer.apertureFittingCut(int(event.xdata))
        self.fittingCanvas.draw()

    def onReleaseAtfittingCanvas(self, event):
        if not event.inaxes: return
        if event.inaxes != self.apertureTracer.apertureAx: return
        if not self.isXcutterPicked: return
        self.isXcutterPicked = False

    def onApertureExtract(self):
        self.spec = self.apertureTracer.apertureExtract(self.regFactor, self.identificationMethod, sumMethod='gauss')
        self.apertureExtractWidget.show()
        self.apertureExtractWidget.raise_()

    def onSaveSpectrumInFit(self):
        specCCD = CCDData(data=self.spec, header=self.hdr, unit="adu")
        if(self.savePath == None):
            self.savePath = str(QFileDialog.getExistingDirectory(self, "Select Directory to Save"))
        self.filePath  = self.savePath+'/spec_'+self.fileName.split('/')[-1]
        specCCD.write(self.filePath, overwrite=True)
        print(f'Saved in {self.filePath}')
        self.close()


    @pyqtSlot(float)
    def onProgressChanged(self, val):
        self.progressBar.setValue(int(val))

'''
if __name__ == '__main__':
    app = QApplication(sys.argv)
    matchList = pd.read_csv("./matchResult.csv")
    _,flux = openFitData("./Spectroscopy_Example/20181023/combine/comp_15_15.0.fit")
    reId = reIdentifier(matchList, flux, fitMethod='chebyshev')
    reId.doFit()
    regFactor = reId.regFactor
    identificationMethod = 'chebyshev'
    FWHM = 4

    ex = apertureTraceWidget("./Spectroscopy_Example/2021-06-19/HD188350-0001.fit", regFactor, identificationMethod)
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = apertureSkyPickerWidget(currentFileLocation = "./Spectroscopy_Example/20181023/reduced/r_HD18247_60.0_1.fit")
    ex.show()
    sys.exit(app.exec_())
'''