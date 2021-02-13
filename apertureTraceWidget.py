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
from fitInfo import cropInfo
from cropWidget import cropWidget
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

    def __init__(self, flux, skyFit, objFit, fig=None, waveFig = None):
        super().__init__()
        self.flux = flux
        self.skyFit = skyFit
        self.objFit = objFit
        self.FWHM = 4 #Todo add method for getting FHWM, gain, readoutnoise from image hdr
        self.cutExist = False
        self.norm = znorm(flux)
        self.apertureFitMethod = 'chebyshev'
        self.pickDistance = 3
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
        self.progressChangeSignal.emit(100)


        skysubtractedFlux = np.array(skysubtractedFlux)
        skyFlux = np.array(skyFlux)
        self.apertures = np.array(apertures)
        self.apertureCoeffs = np.array(apertureCoeffs)
        self.skysubtractedFlux = skysubtractedFlux.T
        self.skyFlux = skyFlux.T
        zimshow(self.subImgAx, self.skysubtractedFlux, normalize=self.norm)
        zimshow(self.skyImgAx, self.skyFlux, normalize=self.norm)

        self.imgAx.set_xlim(0,len(flux[0]))
        self.imgAx.set_ylim(0,len(flux))
        self.apertureAx.set_xlim(0,len(flux[0]))
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
        self.fitAx.plot(xyFlux, skyImage(xyFlux, *popts), 'b--')
        self.fitAx.plot(xyFlux, skyFit(xyFlux, *popts[:2]), 'y--')
        self.fitAx.plot(xyFlux, objFit(xyFlux, *popts[2:]), 'y--')

        self.residualAx.plot(yFlux, fluxNow-skyImage(yFlux, *popts), 'b--')
        self.residualAx.axhline(0, color='r', linestyle='-')

        # fitAx.plot(yFlux, fluxNow, 'r-')
        # fitAx.plot(xyFlux, skyImage(xyFlux, *poptFin), 'b--')

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
                aduSum.append(1 / np.sqrt(2 * np.pi) * apertureCoeffs[x][2])
            wavelength = np.array(wavelength)
            aduSum = np.array(aduSum)
            ax.plot(wavelength, aduSum)

        elif sumMethod == 'real':
            wavelength = []
            aduSum = []
            for x, y in enumerate(fitVal):
                wavelength.append(fittedWavelength(x, y, regFactor, identificationMethod))
                aduSum.append(np.sum(skysubtractedFlux[:, x]))
            wavelength = np.array(wavelength)
            aduSum = np.array(aduSum)

            ax.plot(wavelength, aduSum)

        spec = np.vstack((wavelength,aduSum))

        return spec

#우선 이미지를 받아와서 다시 크롭한다

#Todo need to implement method that reject some point


class apertureTraceWidget(QMainWindow):
    def __init__(self, fileName):
        super().__init__()
        self.FWHM = 2 #Todo add method for getting FHWM, gain, readoutnoise from image hdr
        self.fileName = fileName
        self.flux = fits.open(fileName)[0].data
        self.skyFit = linearSky
        self.objFit = gaussian
        self.apertureTracerOpen = False
        self.isXcutterPicked = False
        self.cropInfo = cropInfo(y0=0,y1=len(self.flux),x0=0,x1=len(self.flux[0]),filename=self.fileName)
        self.norm = znorm(self.flux)
        self.initUI()



    def initUI(self):
        self.layout = QGridLayout()
        self.mainWidget = QWidget()
        self.mainWidget.setLayout(self.layout)
        self.setCentralWidget(self.mainWidget)


        self.imageCanvas = FigureCanvas(Figure()) # contains first image and add points from aperturetrace
        self.mainImageAx = self.imageCanvas.figure.add_subplot(111)


        self.fittingCanvas = FigureCanvas(Figure()) # contains fitting images

        self.waveCanvas = FigureCanvas(Figure())  # contains fitting images



        self.cropBtn = QPushButton('&Crop')
        self.cropWidget = cropWidget(self.fileName)
        self.cropBtn.clicked.connect(self.onCrop)
        self.cropWidget.cropCheckWidget.imageCropSignal.connect(self.imageCrop)

        self.apertureTracerBtn = QPushButton('&ApertureTracer')
        self.apertureTracerBtn.clicked.connect(self.onApertureTracer)

        self.aperturePointBtn = QPushButton('&AperturePoints')
        self.aperturePointBtn.clicked.connect(self.onFindAperturePoints)

        self.apertureTraceBtn = QPushButton('&ApertureTrace')
        self.apertureTraceBtn.clicked.connect(self.onApertureTrace)

        self.progressBar = QProgressBar(self)



        self.layout.addWidget(self.imageCanvas, 0,0, 1, -1)
        self.layout.addWidget(self.cropBtn, 1,0, 1, -1)
        self.layout.addWidget(self.fittingCanvas, 2,0, 1, -1)
        self.layout.addWidget(self.apertureTracerBtn, 3,0)
        self.layout.addWidget(self.aperturePointBtn, 3,1)
        self.layout.addWidget(self.apertureTraceBtn, 3,2)

        self.layout.addWidget(self.progressBar,4,0,1,-1)

        self.showImage()


        self.fittingCanvas.mpl_connect("pick_event", self.onPickAtfittingCanvas)
        self.fittingCanvas.mpl_connect("motion_notify_event", self.onMoveAtfittingCanvas)
        self.fittingCanvas.mpl_connect("button_release_event", self.onReleaseAtfittingCanvas)


    def onCrop(self):
        self.cropWidget.show()
        self.cropWidget.raise_()

    @pyqtSlot(cropInfo)
    def imageCrop(self, crop):
        self.cropWidget.close()
        self.cropInfo = crop
        self.showImage()

    def showImage(self):
        data = self.flux[self.cropInfo.y0:self.cropInfo.y1,
               self.cropInfo.x0:self.cropInfo.x1]
        zimshow(self.mainImageAx, data, normalize=self.norm)
        self.mainImageAx.figure.canvas.draw()


    def onApertureTracer(self):
        self.apertureTracerOpen = True
        self.data  = self.flux[self.cropInfo.y0:self.cropInfo.y1,
               self.cropInfo.x0:self.cropInfo.x1]
        self.apertureTracer = apertureTracer(self.data, self.skyFit, self.objFit, fig=self.fittingCanvas.figure, waveFig = self.waveCanvas.figure)
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
        self.apertureTracer.apertureExtract()

    @pyqtSlot(float)
    def onProgressChanged(self, val):
        self.progressBar.setValue(int(val))


if __name__ == '__main__':
    app = QApplication(sys.argv)

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
    FWHM = 4

    ex = apertureTraceWidget("./Spectroscopy_Example/20181023/reduced/r_HD18247_60.0_1.fit")
    ex.show()
    sys.exit(app.exec_())








#%%
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
FWHM = 4


C = fits.open("./Spectroscopy_Example/20181023/reduced/r_uranus_30.0_1.fit")

skyFit = linearSky
objFit = gaussian

flux = C[0].data[40:110,300:800]


apt = apertureTracer(flux, skyFit, objFit)
apt.findAperturePoints()
apt.apertureFittingCut(387)
apt.apertureTrace()
spec = apt.apertureExtract(regFactor, identificationMethod, sumMethod = 'gauss')




#작은 x에 대해 값이 작아진다. 왜지? preprocessing 문제?



#M = fits.open("./Spectroscopy_Example/20181023/reduced/r_moon_0.1_1.fit")
#flux = M[0].data[40:110,300:800]