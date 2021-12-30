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
from numpy.polynomial.chebyshev import chebfit, chebval
from matplotlib.figure import Figure









def fittedWavelength(x, y, regFactor, fitMethod):
    if fitMethod =='chebyshev':
        return regFactor(x,y)

    elif fitMethod == 'linear':
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
        self.fitOrder = 1
        '''
        do Reidentification with matching information(matchlist) and array of comp image flux(flux)
        default setting is do linear fitting along the x axis of the image
        need to impliment another way of fiting. 
        '''

    def doFit(self, fitMethod=None, yStep=None, order = None):
        print(order)
        if fitMethod is not None:
            self.fitMethod = fitMethod
        if yStep is not None:
            self.yStep = yStep
        if self.fitMethod == 'chebyshev':
            #Todo Make chebnyshevparameterWidget to change order of chevichev/etc < raise it or make it on the corner
            #Make paramater: N_RIED, STEP_REID, ORDER_ID(It maybe is from identification step), ORDER_WAVELEN_REID,ORDER_SPATIAL_REID, TOL_REID
            matchList = self.matchList
            matchList = matchList[matchList.Pixel != 0]
            flux = self.flux
            print(flux)
            N_SPATIAL, N_WAVELEN = np.shape(flux)
            STEP_REID = 10


            if order is not None and len(order)!=0:
                ORDER_ID, ORDER_WAVELEN_REID, ORDER_SPATIAL_REID = np.array(order.split(',')).astype(int)
                self.fitOrder = order
            else :
                ORDER_ID = 4
                ORDER_WAVELEN_REID, ORDER_SPATIAL_REID = 6,6
                self.fitOrder = '4,6,6'

            N_REID = N_SPATIAL // STEP_REID
            TOL_REID = 5
            fwhm = self.FWHM

            fitter = LevMarLSQFitter()
            #Code from https://github.com/ysBach/SNU_AOclass/blob/master/Notebooks/Spectroscopy_Example.ipynb

            line_REID = np.zeros((len(matchList), N_REID - 1))
            spatialcoord = np.arange(0, (N_REID - 1) * STEP_REID, STEP_REID) + STEP_REID / 2

            print('Reidentify each section by Chebyshev (order {:d})'.format(ORDER_ID))
            print('section      |  found  |  RMS')
            wave = []
            for i in range(0, N_REID - 1):
                lower_cut, upper_cut = i * STEP_REID, (i + 1) * STEP_REID
                reidentify_i = np.sum(flux[lower_cut:upper_cut, :],
                                      axis=0)
                peak_gauss_REID = []

                for peak_pix_init in matchList['Pixel']:
                    search_min = int(np.around(peak_pix_init - TOL_REID))
                    search_max = int(np.around(peak_pix_init + TOL_REID))
                    cropped = reidentify_i[search_min:search_max]
                    x_cropped = np.arange(len(cropped)) + search_min

                    A_init = np.max(cropped)
                    mean_init = peak_pix_init
                    stddev_init = fwhm * gaussian_fwhm_to_sigma
                    g_init = Gaussian1D(amplitude=A_init, mean=mean_init, stddev=stddev_init,
                                        bounds={'amplitude': (0, 2 * np.max(cropped)),
                                                'stddev': (0, TOL_REID)})
                    g_fit = fitter(g_init, x_cropped, cropped)
                    fit_center = g_fit.mean.value
                    if abs(fit_center - peak_pix_init) > TOL_REID:
                        peak_gauss_REID.append(np.nan)
                        continue
                    peak_gauss_REID.append(fit_center)

                peak_gauss_REID = np.array(peak_gauss_REID)
                nonan_REID = np.isfinite(peak_gauss_REID)
                line_REID[:, i] = peak_gauss_REID
                peak_gauss_REID_nonan = peak_gauss_REID[nonan_REID]
                n_tot = len(peak_gauss_REID)
                n_found = np.count_nonzero(nonan_REID)


                coeff_REID1D, fitfull = chebfit(peak_gauss_REID_nonan,
                                                    matchList['Wavelength'][nonan_REID],
                                                    deg=ORDER_WAVELEN_REID,
                                                    full=True)
                fitRMS = np.sqrt(fitfull[0][0] / n_found)
                print('[{:04d}:{:04d}]\t{:d}/{:d}\t{:.3f}'.format(lower_cut, upper_cut,
                                                                  n_found, n_tot, fitRMS))
                self.progressChangeSignal.emit(i / (N_REID - 1) * 100)

            points = np.vstack((line_REID.flatten(),
                                np.tile(spatialcoord, len(line_REID))))
            points = points.T  # list of ()
            nanmask = (np.isnan(points[:, 0]) | np.isnan(points[:, 1]))
            points = points[~nanmask]
            values = np.repeat(matchList['Wavelength'], N_REID - 1)
            values = np.array(values.tolist())
            values = values[~nanmask]
            errors = np.ones_like(values)


            coeff_init = Chebyshev2D(x_degree=ORDER_WAVELEN_REID, y_degree=ORDER_SPATIAL_REID)
            fit2D_REID = fitter(coeff_init, points[:, 0], points[:, 1], values)
            ww, ss = np.mgrid[:N_WAVELEN, :N_SPATIAL]
            xPix = points[:,0]
            yPix = points[:,1]
            fittingPoints = pd.DataFrame({'XPixel': np.round(xPix, 4), 'YPixel': yPix, 'Wavelength': values})
            wavelength = fit2D_REID(ww, ss).T
            regFactor = fit2D_REID

        elif self.fitMethod == 'linear':
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
                                                                     fittingPoints.XPixel == x)].Wavelength.iloc[0] - fittedWavelength(
                    x, yCut, self.regFactor, self.fitMethod), 'rx', ms=5)
            yResidualAx.axhline(0, color='blue')

        fig.suptitle("Reidentify and Wavelength Map\nMethod={}, Order=({})".format(self.fitMethod,self.fitOrder))
        imgAx.set_ylabel('Spatial direction(pix)')
        yResidualAx.set_xlabel('Dispersion direction(pix)')
        yResidualAx.set_ylabel('Wavelength residual(Å)')
        yCutAx.set_ylabel('Wavelength(Å)')

        xCutAx.set_xlabel('Wavelength(Å)')

    def setMatchList(self, matchList: pd.DataFrame):
        self.matchList = matchList
    def setFlux(self, flux: np.ndarray):
        self.flux = flux

    def setAll(self, matchList, flux, fitMethod, FWHM):
        self.matchList = matchList
        self.flux = flux
        self.fitMethod = fitMethod
        self.FWHM = FWHM


class reIdentificationWidget(QWidget):
    identificationDoneSignal = pyqtSignal(object, str)
    def __init__(self, matchList: pd.DataFrame, flux: np.ndarray, fitMethod='linear', FWHM=2):
        super().__init__()
        self.flux = flux
        self.reIdentifier = reIdentifier(matchList=matchList, flux=flux, fitMethod=fitMethod, FWHM=FWHM)
        self.reIdentifier.progressChangeSignal.connect(self.onProgressChanged)
        self.initUI()


    # ui 구성
    # 필요한거 : 피팅 시작 버튼(방법-linearfit, 등등의 fitting method에 따른 버튼을 여러개 만들어서
    # 직접 선택할수 있도록 하자.), 로딩창, 결과 imshow로 보여주는 ax 하나, yCut에서 fitting 그래프랑 residual이랑 보여주는 ax 하나.
    # 일단 가운데에 figcanvas를 놓고, 위쪽에 2dplot 3dplot 선택창. 아래쪽에 피팅 시작버튼이랑 로딩창(세로)

    def initUI(self):
        self.layout = QGridLayout()
        # result Figures

        self.resultCanvas = FigureCanvas(Figure(figsize=(10, 6)))

        # Loding bar and fitting button
        self.loadingBar = QProgressBar(self)
        self.linearButton = QPushButton('Start Fitting with &Line Regression',self)
        self.linearButton.clicked.connect(self.onLinearFittingStart)
        self.chebyButton = QPushButton('Start Fitting with &Chebyshev2D',self)
        self.chebyButton.clicked.connect(self.onChebyshevFittingStart)
        self.chebyArgText = QLineEdit(self)
        self.chebyArgText.setPlaceholderText("1dOrder,xOrder,yOrder")
        #Todo Add Another Method If needed

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
        self.xCutSliderLabel = QLabel(f'Xcut : ')
        self.xCutSlider.valueChanged.connect(self.onXCutChanged)

        self.yCutSlider = QSlider(Qt.Vertical, self)
        self.yCutSliderLabel = QLabel(f'Ycut : ')
        self.yCutSlider.valueChanged.connect(self.onYCutChanged)

        self.finishBtn = QPushButton('&Finish Reidentification')
        self.finishBtn.clicked.connect(self.onFinish)

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
        self.layout.addWidget(self.chebyButton, 4, 0)
        self.layout.addWidget(self.chebyArgText, 4, 1)
        self.layout.addWidget(self.loadingBar, 5, 0, 1, -1)
        self.layout.addWidget(self.finishBtn, 6, 0, 1, -1)
        self.setLayout(self.layout)


    def onChebyshevFittingStart(self):
        args = self.chebyArgText.text()
        self.reIdentifier.doFit(fitMethod='chebyshev', yStep=None, order = args)
        self.on2DPlot()

    def onLinearFittingStart(self):
        self.reIdentifier.doFit(fitMethod='linear')
        self.on2DPlot()

    def onYStepSliderChanged(self, val):
        self.yStepSliderLabel.setText(str(val))
        self.reIdentifier.yStep = val

    @pyqtSlot(float)
    def onProgressChanged(self, val):
        self.loadingBar.setValue(val)

    def on2DPlot(self):
        self.setSlider()
        self.reIdentifier.fitPlot2D( fig=self.resultCanvas.figure, Cuts = self.cuts)
        self.resultCanvas.figure.canvas.draw()

    def onXCutChanged(self, val):
        self.cuts = (val, self.cuts[1])
        self.xCutSliderLabel.setText(f'Xcut : {val}')
        self.reIdentifier.fitPlot2D( fig=self.resultCanvas.figure, Cuts = self.cuts)
        self.resultCanvas.figure.canvas.draw()

    def onYCutChanged(self, val):
        val = self.y - val
        if (val%self.reIdentifier.yStep !=0): return
        self.cuts = (self.cuts[0], val)
        self.yCutSliderLabel.setText(f'Ycut : {val}')
        self.reIdentifier.fitPlot2D( fig=self.resultCanvas.figure, Cuts = self.cuts)
        self.resultCanvas.figure.canvas.draw()

    def onFinish(self):
        reply = QMessageBox.question(self, 'Message', 'Use this as Identification?',
                                    QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

        if reply == QMessageBox.Yes:
            self.identificationDoneSignal.emit(self.reIdentifier.regFactor, self.reIdentifier.fitMethod)
            self.close()

    def setReidentifier(self, matchList, flux, fitMethod, FWHM):
        self.flux = flux
        self.reIdentifier.setAll(matchList=matchList, flux=flux, fitMethod=fitMethod, FWHM=FWHM)

    def setSlider(self):
        self.x = len(self.flux[0])
        self.y = len(self.flux)
        self.cuts = int(self.x/2) ,int(self.y/2)
        self.xCutSlider.setValue(self.cuts[0])
        self.xCutSlider.setRange(0, self.x-1)
        self.yCutSlider.setValue(self.cuts[1])
        self.yCutSlider.setRange(0, self.y-1)




if __name__ == '__main__':
    app = QApplication(sys.argv)

    matchList = pd.DataFrame(dict(Pixel=[403, 374, 362, 265,
                                         245, 238, 223, 213,
                                         178, 169, 144, 131, 53],
                                  Wavelength=[8780.6, 8495.4, 8377.6, 7438.9,
                                              7245.2, 7173.9, 7032.4, 6929.5,
                                              6599.0, 6507.0, 6266.5, 6143.1, 5400.6]))
    _,flux = openFitData("./Spectroscopy_Example/20181023/combine/comp_15_15.0.fit")
    ex = reIdentificationWidget(matchList, flux)
    ex.show()
    sys.exit(app.exec_())

