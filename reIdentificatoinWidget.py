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
import math
from scipy import stats






#Todo cross correlation을 사용해 (FFT를 그려서)같은부분을 자동으로 찾아주는 알고리즘을 넣어보자.
#vs
#Todo https://github.com/jveitchmichaelis/rascal/blob/master/rascal/calibrator.py hough transform
# 허프 트렌스폼은 직선을 찾는 알고리즘인데 상대적으로 많은 Wavelength emission값이랑 상대적으로 적은
# pixel peak 값을 매칭해주는 직선(wavelength/pixel 그래프의 직선)을 찾아주는거 같다
# 참고해서 자동매칭되게 하고 그 매칭된 정보를 직접 수정해 사용할 수 있도록 하자

#Todo 다른 피팅방식(체비세프라던가)을 적용해 보자.


def fittedWavelength(x, y, regFactor, fitMethod):
    if fitMethod == 'linear':
        if type(y) is not (list or np.ndarray) : y = [y]

        wavelength = []
        for yVal in y:
            yHigh = math.ceil(yVal)
            yLow = math.floor(yVal)
            yDiff = yVal - yLow
            slopeHigh = regFactor['Slope'][yHigh]
            interceptHigh = regFactor['Intercept'][yHigh]
            slopeLow = regFactor['Slope'][yLow]
            interceptLow = regFactor['Intercept'][yLow]
            wavelengthValHigh = x * slopeHigh + interceptHigh
            wavelengthValLow = x * slopeLow + interceptLow
            wavelength.append((wavelengthValHigh - wavelengthValLow)*yDiff + wavelengthValLow)
        wavelength = np.array(wavelength)

    if len(wavelength==1):wavelength = wavelength[0]
    return wavelength



class reIdentifier(QWidget):
    progressChangeSignal = pyqtSignal(float)

    def __init__(self, matchList: pd.DataFrame, flux: np.array, fitMethod='linear', FWHM=2):
        super().__init__()
        self.matchList = matchList
        self.fitMethod = fitMethod
        self.flux = flux
        self.yStep = 1
        self.regFactor = None
        self.wavelength = None
        self.fittingPoints = pd.DataFrame({'XPixel': [], 'YPixel': [], 'Wavelength': []})
        self.FWHM = FWHM
        '''
        do Reidentification with matching information(matchlist) and array of comp image flux(flux)
        default setting is do linear fitting along the x axis of the image
        need to impliment another way of fiting. 
        '''

    def doFit(self, fitMethod=None, yStep=None):
        if fitMethod is not None:
            self.fitMethod = fitMethod
        if yStep is not None:
            self.yStep = yStep
        if self.fitMethod == 'linear':
            matchList = self.matchList
            matchList = matchList[matchList.Pixel != 0]
            ystep = self.yStep
            data = self.flux
            fwhm = self.FWHM
            regFactor = pd.DataFrame({'Slope': [], 'Intercept': []})
            fittingPoints = pd.DataFrame({'XPixel': [], 'YPixel': [], 'Wavelength': []})
            wavelength = []
            xPixel = np.arange(len(data[0]))
            fitter = LevMarLSQFitter()
            max = int(len(data) / ystep)

            for i in np.arange(len(data) / ystep):
                i = int(i)
                if (ystep == 1):
                    dat = data[i, :]
                elif (i * ystep + ystep - 1 > len(data)):
                    dat = np.average(data[ystep * i:len(data) - 1, :], axis=0)
                else:
                    dat = np.average(data[ystep * i:ystep * i + ystep - 1, :], axis=0)
                ground = np.median(dat[0:100])
                dat_fit = dat - ground
                peakGauss = []
                for peakPix in matchList['Pixel']:
                    peakPix = int(peakPix)
                    g_init = Gaussian1D(amplitude=dat_fit[peakPix],
                                        mean=peakPix,
                                        stddev=fwhm * gaussian_fwhm_to_sigma,
                                        bounds={'amplitude': (dat_fit[peakPix], 2 * dat_fit[peakPix]),
                                                'mean': (peakPix - fwhm, peakPix + fwhm),
                                                'stddev': (0, fwhm)})
                    fitted = fitter(g_init, xPixel, dat_fit)

                    peakGauss.append(fitted.mean.value)
                    dat_fit = dat_fit - fitted(xPixel)
                self.progressChangeSignal.emit(i / max * 100)
                res = stats.linregress(peakGauss, matchList['Wavelength'])
                regAdd = pd.DataFrame({'Slope': [res.slope], 'Intercept': [res.intercept]})
                pointAdd = pd.DataFrame({'XPixel': np.round(peakGauss, 4), 'YPixel': np.full(len(peakGauss), i * ystep),
                                         'Wavelength': matchList['Wavelength']})
                for j in range(ystep):
                    if len(regFactor) >= len(data): continue
                    regFactor = regFactor.append(regAdd, ignore_index=True)
                    fittingPoints = fittingPoints.append(pointAdd, ignore_index=True)
                    wave = xPixel * res.slope + res.intercept
                    wave = np.round(wave, 4)
                    wavelength.append(wave)
        self.progressChangeSignal.emit(100)
        self.fittingPoints = fittingPoints
        self.regFactor = regFactor
        self.wavelength = np.array(wavelength)
        '''
        xPixel, YPixel, wavelength를 X,Y,Z 축으로 가지는 3축 좌표평면에 실제 fitting 이 일어난 각 ypixel과  gauss
        fitting을 통해 직접 fitting에 사용된 xpixel을 x,y 값으로, matchList['Wavelength'] 를 z값으로 사용해 
        처음 데이터를 plot하고 최종값을 plot해 fitting result와 그 residual을 보여준다. 

        '''

    def fitPlot3D(self, fig=None):
        fittingPoints = self.fittingPoints
        wavelength = self.wavelength
        if fig is None:
            fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        xPix = fittingPoints['XPixel']
        yPix = fittingPoints['YPixel']
        zPix = fittingPoints['Wavelength']
        ax.scatter(xPix, yPix, zPix, c='r', marker='o', s=1)
        xSpace = np.arange(len(wavelength[0]))
        ySpace = np.arange(len(wavelength))
        X, Y = np.meshgrid(xSpace, ySpace)
        ax.plot_surface(X, Y, wavelength, cmap='viridis')
        axis = len(wavelength[0]) if len(wavelength[0]) >= len(wavelength) else len(wavelength)

        ax.set_xlim3d(0, axis)
        ax.set_ylim3d(0, axis)

    def fitPlot2D(self, fig=None, Cuts: tuple = None):

        fittingPoints = self.fittingPoints
        wavelength = self.wavelength

        if fig is None:
            fig = plt.figure()
        else:
            fig.clear()
        if Cuts is None:
            imgAx = fig.add_subplot(111)
        else:
            grid = fig.add_gridspec(ncols=2, nrows=3, width_ratios=[10, 2],
                                    height_ratios=[2, 2, 2])
            imgAx = fig.add_subplot(grid[0, 0])
            xCut, yCut = Cuts
        xPix = fittingPoints['XPixel']
        yPix = fittingPoints['YPixel']

        imgAx.plot(xPix, yPix, 'rx', ms=1)

        imgAx.imshow(wavelength, cmap='viridis')

        if Cuts is not None:
            imgAx.set_aspect('auto')
            xCutAx = fig.add_subplot(grid[0, 1], sharey=imgAx)
            xCutAx.plot(wavelength[:, xCut], np.arange(len(wavelength[:, xCut])), 'b')
            self.xLine = imgAx.axvline(xCut, color='black')

            yCutAx = fig.add_subplot(grid[1, 0], sharex=imgAx)
            yResidualAx = fig.add_subplot(grid[2, 0], sharex=imgAx)
            yCutAx.plot(wavelength[yCut], 'b')
            yCutAx.plot(fittingPoints.loc[fittingPoints.YPixel == yCut].XPixel,
                        fittingPoints.loc[fittingPoints.YPixel == yCut].Wavelength, 'rx', ms=5)
            self.yLine = imgAx.axhline(yCut, color='black')
            for x in fittingPoints.loc[fittingPoints.YPixel == yCut].XPixel:
                yResidualAx.plot(x, fittingPoints.loc[np.logical_and(fittingPoints.YPixel == yCut,
                                                                     fittingPoints.XPixel == x)].Wavelength - fittedWavelength(
                    x, yCut, self.regFactor, self.fitMethod), 'rx', ms=5)
            yResidualAx.axhline(0, color='blue')

    def setMatchList(self, matchList: pd.DataFrame):
        self.matchList = matchList
    def setFlux(self, flux: np.ndarray):
        self.flux = flux

class reIdentificationWidget(QWidget):
    def __init__(self, matchList: pd.DataFrame, flux: np.ndarray, fitMethod='linear', FWHM=2):
        super().__init__()
        self.flux = flux
        self.reIdentifier = reIdentifier(matchList=matchList, flux=flux, fitMethod=fitMethod, FWHM=FWHM)
        self.reIdentifier.progressChangeSignal.connect(self.onProgressChanged)
        self.x = len(self.flux[0])
        self.y = len(self.flux)
        self.cuts = int(self.x/2) ,int(self.y/2)
        self.initUI()


    # ui 구성
    # 필요한거 : 피팅 시작 버튼(방법-linearfit, 등등의 fitting method에 따른 버튼을 여러개 만들어서
    # 직접 선택할수 있도록 하자.), 로딩창, 결과 imshow로 보여주는 ax 하나, yCut에서 fitting 그래프랑 residual이랑 보여주는 ax 하나.
    # 일단 가운데에 figcanvas를 놓고, 위쪽에 2dplot 3dplot 선택창. 아래쪽에 피팅 시작버튼이랑 로딩창(세로)

    def initUI(self):
        self.layout = QGridLayout()
        # result Figures
        self.resultFig = plt.Figure(figsize=(10, 6))
        self.resultCanvas = FigureCanvas(self.resultFig)

        # Loding bar and fitting button
        self.loadingBar = QProgressBar(self)
        self.linearButton = QPushButton('Start Fitting with &Line Regression',self)
        self.linearButton.clicked.connect(self.onLinearFittingStart)

        self.yStepSlider = QSlider(Qt.Horizontal, self)
        self.yStepSlider.setValue(1)
        self.yStepSlider.setRange(1, 10)
        self.yStepSliderLabel = QLabel(str(1))
        self.yStepSlider.valueChanged.connect(self.onYStepSliderChanged)

        #2D and 3D plot buttons (not working for now)
        self.plot2dBtn = QPushButton('&2D Plot')
        self.plot3dBtn = QPushButton('&3D Plot')
        self.plot2dBtn.clicked.connect(self.on2DPlot)

        #x, y cut sliders
        self.xCutSlider = QSlider(Qt.Vertical, self)
        self.xCutSlider.setValue(self.cuts[0])
        self.xCutSlider.setRange(0, self.x-1)
        self.xCutSliderLabel = QLabel(f'Xcut : {self.cuts[0]}')
        self.xCutSlider.valueChanged.connect(self.onXCutChanged)

        self.yCutSlider = QSlider(Qt.Vertical, self)
        self.yCutSlider.setValue(self.cuts[1])
        self.yCutSlider.setRange(0, self.y-1)
        self.yCutSliderLabel = QLabel(f'Ycut : {self.cuts[1]}')
        self.yCutSlider.valueChanged.connect(self.onYCutChanged)


        self.layout.addWidget(self.plot2dBtn, 0, 0)
#        self.layout.addWidget(self.plot3dBtn, 0, 1)
        self.layout.addWidget(self.xCutSliderLabel, 0, 1)
        self.layout.addWidget(self.yCutSliderLabel, 0, 2)

        self.layout.addWidget(self.resultCanvas, 1, 0)
        self.layout.addWidget(self.xCutSlider, 1, 1)
        self.layout.addWidget(self.yCutSlider, 1, 2)
        self.layout.addWidget(self.yStepSlider, 2, 0)
        self.layout.addWidget(self.yStepSliderLabel, 2, 1)
        self.layout.addWidget(self.linearButton, 3, 0)
        self.layout.addWidget(self.loadingBar, 4, 0, 1, -1)
        self.setLayout(self.layout)




    def onLinearFittingStart(self):
        self.reIdentifier.doFit(fitMethod='linear')

    def onYStepSliderChanged(self, val):
        self.yStepSliderLabel.setText(str(val))
        self.reIdentifier.yStep = val

    @pyqtSlot(float)
    def onProgressChanged(self, val):
        self.loadingBar.setValue(val)

    def on2DPlot(self):
        self.reIdentifier.fitPlot2D( fig=self.resultFig, Cuts = self.cuts)
        self.resultFig.canvas.draw()
    def onXCutChanged(self, val):
        self.cuts = (val, self.cuts[1])
        self.xCutSliderLabel.setText(f'Xcut : {val}')
        self.reIdentifier.fitPlot2D( fig=self.resultFig, Cuts = self.cuts)
        self.resultFig.canvas.draw()
    def onYCutChanged(self, val):
        val = self.y - val
        self.cuts = (self.cuts[0], val)
        self.yCutSliderLabel.setText(f'Ycut : {val}')
        self.reIdentifier.fitPlot2D( fig=self.resultFig, Cuts = self.cuts)
        self.resultFig.canvas.draw()
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
    ex = reIdentificationWidget(matchList, flux)
    ex.show()
    sys.exit(app.exec_())

