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
from scipy import stats
#우선 있는 파일 중에 Cali 파일만 열어서 first identification 하는 프로그램
#cali fileinfo를 받아와서 연 다음에

#Todo 지금 있는 rawSpectrum 임시로 그린거니까 다시그리기
# -> neonarc.fit은 지금 있는 cali파일에서 아이덴티피케이션만 해서 그린거라 부정확, 신뢰성있는 스펙트럼으로 바꿔서 넣자.

#Todo 패턴매칭 알고리즘을 사용해 자동으로 비슷한 이미지구간을 찾아주는 기능도 구현해보자


#UI 디자인
#우선 새 창을 열면 위아래로 구성된 plt 캔버스 두개랑 오른쪽에 파일 여는 버튼이 나온다


#위에는 직접 찍은 identification 이미지를 열고 그 오른쪽에 이미지의 y 축에 따라 사용할 데이터의 범위를 정할 수 있게
#이미지를 보여주는 창을 만든다.
#아래는 이미 존재하는 identification 이미지를 불러온다. 이때 이미 저장된 preload 이미지를 사용하는 방법과
#Todo 새로운 스펙트럼을 열어서 비교하는 방법을 사용한다.

class identificationWidget(QMainWindow):
    def __init__(self):
        super().__init__()
        self.peakThreshold = 0.01
        self.peakNumber =  50
        self.peakDistance = 2
        self.initUI()


    def initUI(self):
        self.layout = QGridLayout()
        self.mainWidget = QWidget()

        self.calibrationFileOpenAction = QAction('CalibrationFileOpen', self)
        self.calibrationFileOpenAction.setShortcut('Ctrl+C')
        self.calibrationFileOpenAction.triggered.connect(self.onCalibnationFileOpen)
        self.calibrationFileOpenAction.setStatusTip('Open calibration image')

        #직접 찍은 아이덴티피케이션 이미지의 스펙트럼을 보여주는 fig

        self.selfSpectrumFig = plt.Figure(figsize = (15,5))
        self.selfSpectrumCanvas = FigureCanvas(self.selfSpectrumFig)

        self.selfSpectrumFig.clear()

        self.peakNumberSlider = QSlider(Qt.Horizontal, self)
        self.peakNumberSlider.setValue(50)
        self.peakNumberSlider.setRange(1,100)
        self.peakDistanceSlider = QSlider(Qt.Horizontal, self)
        self.peakDistanceSlider.setValue(2)
        self.peakNumberSlider.setRange(1, 10)
        self.peakThresholdSlider = QSlider(Qt.Horizontal, self)
        self.peakThresholdSlider.setValue(1)
        self.peakNumberSlider.setRange(1, 100)

        self.peakNumberLabel = QLabel(f'Number of Peak = {self.peakNumber}')
        self.peakDistanceLabel = QLabel(f'Distance between Peak = {self.peakDistance}')
        self.peakThresholdLabel = QLabel(f'Threshold of peak = {self.peakThreshold}')

        self.peakNumberSlider.valueChanged.connect(self.onPeakNumberValueChanged)
        self.peakDistanceSlider.valueChanged.connect(self.onPeakDistanceValueChanged)
        self.peakThresholdSlider.valueChanged.connect(self.onPeakThresholdValueChanged)

        self.selfPeakControl = QWidget()
        self.peakControlLayout = QVBoxLayout()
        self.peakControlLayout.addWidget(self.peakNumberLabel)
        self.peakControlLayout.addWidget(self.peakNumberSlider)
        self.peakControlLayout.addWidget(self.peakDistanceLabel)
        self.peakControlLayout.addWidget(self.peakDistanceSlider)
        self.peakControlLayout.addWidget(self.peakThresholdLabel)
        self.peakControlLayout.addWidget(self.peakThresholdSlider)

        self.selfPeakControl.setLayout(self.peakControlLayout)



        # 직접 찍은 아이덴티피케이션 이미지를 보여주는 fig

        self.selfImageFig = plt.Figure(figsize = (5,2))
        self.selfImageCanvas = FigureCanvas(self.selfImageFig)

        self.selfImageFig.clear()

        self.selfImageCanvas.mpl_connect("button_press_event", self.onPressAtImage)
        self.selfImageCanvas.mpl_connect("motion_notify_event", self.onMoveAtImage)
        self.selfImageCanvas.mpl_connect("button_release_event", self.onReleaseAtImage)

        self.NeonArcButton = QPushButton('&Neon')
        self.NeonArcButton.clicked.connect(self.neonSpectrumDraw)
        self.OpenArcButton = QPushButton('&Open')

        self.standardSpectrumButtonLayout = QVBoxLayout()
        self.standardSpectrumButtonLayout.addWidget(self.NeonArcButton)
        self.standardSpectrumButtonLayout.addWidget(self.OpenArcButton)

        self.standardSpectrumButton = QWidget()
        self.standardSpectrumButton.setLayout(self.standardSpectrumButtonLayout)



        #비교할 아이덴티피케이션의 스펙트럼을 보여주는 fig
        self.standardSpectrumFig = plt.Figure(figsize = (15,5))
        self.standardSpectrumCanvas = FigureCanvas(self.standardSpectrumFig)

        self.standardSpectrumFig.clear()
        standardSpectrumAx = self.standardSpectrumFig.add_subplot(111)

        self.setCentralWidget(self.mainWidget)
        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)
        filemenu = menubar.addMenu('&File')#&는 File을 Alt F로 실행하게 해준다
        filemenu.addAction(self.calibrationFileOpenAction)



        self.layout.addWidget(self.selfSpectrumCanvas, 1, 0, 2, 1)
        self.layout.addWidget(self.selfImageCanvas,1,1,1,1)
        self.layout.addWidget(self.selfPeakControl, 2, 1, 1, 1)
        self.layout.addWidget(self.standardSpectrumCanvas,3,0)
        self.layout.addWidget(self.standardSpectrumButton,3,1)


        self.mainWidget.setLayout(self.layout)
        self.setCentralWidget(self.mainWidget)


    #opens file and
    def onCalibnationFileOpen(self):
        filePath = QFileDialog.getOpenFileName(self, 'Open calibration file','./Spectroscopy_Example/20181023/combine/')[0]
        hdr, data = openFitData(filePath)
        self.selfImageFig.clear()
        self.selfImageAx = self.selfImageFig.add_subplot(111)
        zimshow(self.selfImageAx, data)
        self.selfImageCanvas.draw()

        self.imageWidth = int (data.shape[1])
        self.selfData = data

    def onPeakNumberValueChanged(self, val):
        self.peakNumber = val
        self.peakNumberLabel.setText(f'Number of Peak = {self.peakNumber}')
        self.selfSpectrumDraw(ymin = self.selfImageY[0], ymax= self.selfImageY[1], data = self.selfData, args=[self.peakDistance,self.peakThreshold,self.peakNumber])

    def onPeakDistanceValueChanged(self, val):
        self.peakDistance = val
        self.peakDistanceLabel.setText(f'Distance between Peak = {self.peakDistance}')

        self.selfSpectrumDraw(ymin = self.selfImageY[0], ymax= self.selfImageY[1], data = self.selfData, args=[self.peakDistance,self.peakThreshold,self.peakNumber])

    def onPeakThresholdValueChanged(self, val):
        self.peakThreshold = val/100
        self.peakThresholdLabel.setText(f'Threshold of peak = {self.peakThreshold}')
        self.selfSpectrumDraw(ymin = self.selfImageY[0], ymax= self.selfImageY[1], data = self.selfData, args=[self.peakDistance,self.peakThreshold,self.peakNumber])





    def onPressAtImage(self, event):

        self.rect = Rectangle((0, 0), 1, 1, alpha=0.5)
        self.selfImageAx.add_patch(self.rect)
        self.x0 = event.xdata
        self.y0 = event.ydata

    def onMoveAtImage(self, event):
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.imageWidth)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((0, self.y0))
        self.selfImageAx.figure.canvas.draw()

    def onReleaseAtImage(self, event):
        y = int(self.rect.get_y())
        height = int(self.rect.get_height())
        self.rect.remove()
        self.selfImageAx.figure.canvas.draw()
        if (height<0) :
            height = 0-height
        self.selfImageY = np.array([y, y+height])
        self.selfSpectrumDraw(ymin = y, ymax = y+height, data = self.selfData, args = [self.peakDistance, self.peakThreshold, self.peakNumber])




    def selfSpectrumDraw(self, ymin, ymax, data, args):

        MINSEP_PK = args[0]  # minimum separation of peaks
        MINAMP_PK = args[1]  # fraction of minimum amplitude (wrt maximum) to regard as peak
        NMAX_PK = args[2]

        self.selfSpectrumFig.clear()
        self.selfSpectrumAx = self.selfSpectrumFig.add_subplot(111)
        identify = np.average(data[ymin:ymax, :], axis=0)
        ground = np.median(identify[0:200])
        max_intens = np.max(identify)
        peak_pix = peak_local_max(identify, indices=True, num_peaks=NMAX_PK,
                                  min_distance=MINSEP_PK,
                                  threshold_abs=max_intens * MINAMP_PK + ground)
        for i in peak_pix:



            self.selfSpectrumAx.axvline(i, identify[i] / max(identify) + 0.0003, identify[i] / max(identify) + 0.2, color='c')
            if (identify[i] + max(identify) / 2.85> max(identify)):
                self.selfSpectrumAx.text(i, identify[i] + max(identify) / 2000, str(i), clip_on=False)
            else:
                self.selfSpectrumAx.text(i, identify[i] + max(identify) / 2.85, str(i), ha= 'center', va= 'center',
                                         rotation=90, clip_on=True)

        self.selfSpectrumAx.plot(identify, color='r')
        self.selfSpectrumAx.set_ylim(0, )
        plt.tight_layout()
        self.selfSpectrumAx.figure.canvas.draw()

    def neonSpectrumDraw(self):
        filePath = './NeonArcSpectrum.fit'
        hdr, data = openFitData(filePath)
        self.standardSpectrumDraw(data = data, arc = 'Neon')



    def standardSpectrumDraw(self, data, arc, peaks=[]):
        wavelength = data[0]
        flux = data[1]
        self.standardSpectrumFig.clear()
        self.standardSpectrumAx = self.standardSpectrumFig.add_subplot(111)
        if (arc=='Neon'):
            peaks = [5330.8000, 5400.5620, 5764.4180, 5852.4878, 5944.8342, 6029.9971, 6074.3377, 6096.1630,  6143.0623,
                 6163.5939, 6217.2813, 6266.4950,  6304.7892,  6334.4279, 6382.9914,  6402.2460, 6506.5279,  6532.8824,
                 6598.9529,  6717.0428,  6929.4680,  7032.4127,7173.9390 , 7245.1670, 7438.8990 , 7488.8720, 7535.7750,
                 8082.4580, 8377.6070 ]

        for i in peaks:
            self.standardSpectrumAx.axvline(i, flux[np.where(i == wavelength)][0] / max(flux) + 0.003,
                        flux[np.where(i == wavelength)] / max(flux) + 0.2, color='c')
            if (flux[np.where(i == wavelength)][0] + max(flux) / 2.85 > max(flux)):
                self.standardSpectrumAx.text(i, flux[np.where(i == wavelength)][0] + max(flux) / 2000, str(i), clip_on=False)
            else:
                self.standardSpectrumAx.text(i, flux[np.where(i == wavelength)][0] + max(flux) / 2.85, str(i),
                         ha='center', va='center', rotation=90, clip_on=True)
        self.standardSpectrumAx.plot(wavelength, flux, color='r')
        self.standardSpectrumAx.set_ylim(0, )
        plt.tight_layout()
        self.standardSpectrumAx.figure.canvas.draw()

    def onButtonClicked(self, status):
        print(self)
        self.bottonSinal.emit(status)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = identificationWidget()
    ex.show()
    sys.exit(app.exec_())




B = fits.open("./Spectroscopy_Example/20181023/combine/comp_15_15.0.fit")

identify = np.average(B[0].data[50:70,:], axis=0)
max_intens = np.max(identify)
MINSEP_PK = 5   # minimum separation of peaks
MINAMP_PK = 0.01 # fraction of minimum amplitude (wrt maximum) to regard as peak
NMAX_PK = 50
FWHM_ID = 4

peak_pix = peak_local_max(identify, indices=True, num_peaks=NMAX_PK,
                          min_distance=MINSEP_PK,
                          threshold_abs=max_intens * MINAMP_PK)

for i in peak_pix:
    plt.axvline(i,0,100000, color = 'c')
    plt.text(i,40000,str(i), rotation=70)
plt.plot(identify, color = 'r')
plt.xlim(400,800)
plt.show()





ID_init = dict (pixel_init=[403, 374, 362, 265,
                           245, 238, 223, 213,
                           178, 169, 144, 131, 53],
wavelength=[8780.6, 8495.4, 8377.6, 7438.9,
                           7245.2, 7173.9, 7032.4, 6929.5,
                           6599.0, 6507.0, 6266.5, 6143.1, 5400.6])

fitter = LevMarLSQFitter()

x_identify = np.arange(len(identify))
peak_gauss = []
for peak_pix in ID_init['pixel_init']:
    #TODO: put something like "lost_factor" to multiply to FWHM_ID in the bounds.
    g_init = Gaussian1D(amplitude = identify[peak_pix],
                       mean = peak_pix,
                       stddev = FWHM_ID * gaussian_fwhm_to_sigma,
                       bounds={'amplitude': (0, 2*identify_1[peak_pix]),
                               'mean':(peak_pix-FWHM_ID, peak_pix+FWHM_ID),
                               'stddev':(0, FWHM_ID)})
    fitted = fitter(g_init, x_identify, identify)
    peak_gauss.append(fitted.mean.value)

res = stats.linregress(peak_gauss, ID_init['wavelength'])
wavelength  = x_identify*res.slope + res.intersections

for i in peak_gauss:
    plt.axvline(i,0,100000, color = 'c')
    plt.text(i,40000,str(i), rotation=70)
plt.plot(identify, color = 'r')
plt.xlim(400,800)
plt.show()

