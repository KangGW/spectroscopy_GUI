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

from scipy import stats
#우선 있는 파일 중에 Cali 파일만 열어서 first identification 하는 프로그램
#cali fileinfo를 받아와서 연 다음에

#Todo 다른 Arc Spectrum 들 찾아서 넣기 NeonArcSpectrum.fit은 peak wavelength 맞춰놔서 Ientify 용도로는 신뢰성있게 쓸수 있다.




#Todo 용도별로 Widget 묶고 method도 묶어서 보기 편하게 하자.


#UI 디자인
#우선 새 창을 열면 위아래로 구성된 plt 캔버스 두개랑 오른쪽에 파일 여는 버튼이 나온다


#위에는 직접 찍은 identification 이미지를 열고 그 오른쪽에 이미지의 y 축에 따라 사용할 데이터의 범위를 정할 수 있게
#이미지를 보여주는 창을 만든다.
#아래는 이미 존재하는 identification 이미지를 불러온다. 이때 이미 저장된 preload 이미지를 사용하는 방법과

#Todo 새로운 스펙트럼을 열어서 비교하는 방법을 사용한다. 이때 https://physics.nist.gov/의 정보를 받아와서 스펙트럼을 만드는 방법을 생각해 보자(기본
# 적으로 Relative intensity가 있으니 스펙트럼화해서 쓰면 될 것이다!


#호에엥
#1. standardSpectrumCanvas가 선택된 상태에서 axvline을 더블클릭하면 그 라인의 색깔이 바뀌고,
#2. selfSpectrumCanvas에 같은 색깔의 긴 axvline이 생긴다. 이건 마우스 따라 움직이고  selfSpectrumCanvas axvline 근처에 가져다되면 달라붙는다.
#3. 마우스 우클릭시 과정 초기화-> 1. 전으로 이동
#4. 달라붙었을때 더블클릭하면 해당 정보가 저장되고 왼쪽 테이블에 뜬다.



class identificationWidget(QMainWindow):
    def __init__(self):
        super().__init__()
        self.peakThreshold = 0.01
        self.peakNumber =  50
        self.peakDistance = 2
        self.pickDistance = 1
        self.wavelengthPixelList = pd.DataFrame( {'Wavelength':[' '], 'Pixel':[' ']})
        self.selfPeakPixs = []
        self.standardPeakWavelengths = []
        self.matchedPeakPixs = []
        self.isPressed = False
        self.isPicked = False
        self.isMatchFinished = True
        self.standardSpectrum = []
        self.selfSpectrum = []
        self.currentPickedPeakWavelength = 0
        self.selfFWHM = 2
        self.REIDYStep = 2
        self.selfImageY = [0,0]
        self.initUI()

    def initUI(self):
        self.layout = QGridLayout()
        self.mainWidget = QWidget()

        self.calibrationFileOpenAction = QAction('CalibrationFileOpen', self)
        self.calibrationFileOpenAction.setShortcut('Ctrl+C')
        self.calibrationFileOpenAction.triggered.connect(self.onCalibnationFileOpen)
        self.calibrationFileOpenAction.setStatusTip('Open calibration image')

        #직접 찍은 아이덴티피케이션 이미지의 스펙트럼을 보여주는 fig

        self.selfSpectrumFig = plt.Figure(figsize = (13,5))
        self.selfSpectrumCanvas = FigureCanvas(self.selfSpectrumFig)

        self.selfSpectrumFig.clear()

        self.peakNumberSlider = QSlider(Qt.Horizontal, self)
        self.peakNumberSlider.setValue(self.peakNumber)
        self.peakNumberSlider.setRange(1,100)
        self.peakDistanceSlider = QSlider(Qt.Horizontal, self)
        self.peakDistanceSlider.setValue(self.peakDistance)
        self.peakNumberSlider.setRange(1, 10)
        self.peakThresholdSlider = QSlider(Qt.Horizontal, self)
        self.peakThresholdSlider.setValue(int(self.peakThreshold*100))
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

        self.selfImageCanvas.mpl_connect("button_press_event", self.onPressAtImage)
        self.selfImageCanvas.mpl_connect("motion_notify_event", self.onMoveAtImage)
        self.selfImageCanvas.mpl_connect("button_release_event", self.onReleaseAtImage)


        self.selfSpectrumCanvas.mpl_connect('scroll_event', self.onScrollAtSelfSpectrum)
        self.selfSpectrumCanvas.mpl_connect('pick_event', self.onPickPeakAtSelfSpectrum)
        self.selfSpectrumCanvas.mpl_connect("button_press_event", self.onPressAtSelfSpectrum)
        self.selfSpectrumCanvas.mpl_connect("motion_notify_event", self.onMoveAtSelfSpectrum)
        self.selfSpectrumCanvas.mpl_connect("button_release_event", self.onReleaseAtSelfSpectrum)


        self.selfSpectrumGaussFitFig = plt.Figure(figsize=(7, 7))
        self.selfSpectrumGaussFitCanvas = FigureCanvas(self.selfSpectrumGaussFitFig)
        self.gaussFitWidget = QWidget()
        self.gaussFitLayout = QVBoxLayout()
        self.gaussFitButton = QPushButton('&Yes')
        self.gaussFitButton.clicked.connect(self.onGaussFitButtonClicked)
        self.FWHMSlider = QSlider(Qt.Horizontal, self)
        self.FWHMSlider.setValue(self.selfFWHM*10)
        self.FWHMSlider.setRange(1, 100)

        self.FWHMLabel = QLabel(f'FHWM for comp image = {self.selfFWHM}')
        self.FWHMSlider.valueChanged.connect(self.onFWHMChanged)
        self.gaussFitLayout.addWidget(self.selfSpectrumGaussFitCanvas)
        self.gaussFitLayout.addWidget(self.FWHMSlider)
        self.gaussFitLayout.addWidget(self.FWHMLabel)
        self.gaussFitLayout.addWidget(self.gaussFitButton)
        self.gaussFitWidget.setLayout(self.gaussFitLayout)




        self.NeonArcButton = QPushButton('&Neon')
        self.NeonArcButton.clicked.connect(self.neonSpectrumDraw)
        self.OpenArcButton = QPushButton('&Open')

        self.standardSpectrumButtonLayout = QVBoxLayout()
        self.standardSpectrumButtonLayout.addWidget(self.NeonArcButton)
        self.standardSpectrumButtonLayout.addWidget(self.OpenArcButton)

        self.standardSpectrumButton = QWidget()
        self.standardSpectrumButton.setLayout(self.standardSpectrumButtonLayout)



        #비교할 아이덴티피케이션의 스펙트럼을 보여주는 fig
        self.standardSpectrumFig = plt.Figure(figsize = (13,5))
        self.standardSpectrumCanvas = FigureCanvas(self.standardSpectrumFig)

        self.standardSpectrumCanvas.mpl_connect('scroll_event', self.onScrollAtStandardSpectrum)
        self.standardSpectrumCanvas.mpl_connect('pick_event', self.onPickPeakAtStandardSpectrum)
        self.standardSpectrumCanvas.mpl_connect('button_press_event', self.onPressAtStandardSpectrum)

        self.wavelengthPixelTable = QTableView()
        self.wavelengthPixelModel = tableModel(self.wavelengthPixelList)
        self.wavelengthPixelTable.setModel(self.wavelengthPixelModel)
        self.wavelengthPixelTable.setSelectionBehavior(QTableView.SelectRows)
        self.wavelengthPixelTable.doubleClicked.connect(self.onWavelengthPixelTableDoubleClicked)



        self.gaussButton = QPushButton('&GuassFit')
        self.gaussButton.clicked.connect(self.selfSpectrumDrawWithGauss)

        self.matchButton = QPushButton('&Match')
        self.matchButton.clicked.connect(self.onMatch)
        self.abortButton = QPushButton('&Abort')
        self.abortButton.clicked.connect(self.onAbort)
        self.exportButton = QPushButton('&Export')
        self.exportButton.clicked.connect(self.onExport)
        self.importButton = QPushButton('&Import')
        self.importButton.clicked.connect(self.onImport)


        self.tableMatchingButtons = QWidget()
        self.tableMatchingButtonLayout = QVBoxLayout()
        self.tableMatchingButtonLayout.addWidget(self.matchButton)
        self.tableMatchingButtonLayout.addWidget(self.abortButton)
        self.tableMatchingButtonLayout.addWidget(self.exportButton)
        self.tableMatchingButtonLayout.addWidget(self.importButton)
        self.tableMatchingButtons.setLayout(self.tableMatchingButtonLayout)




        self.setCentralWidget(self.mainWidget)
        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)
        filemenu = menubar.addMenu('&File')#&는 File을 Alt F로 실행하게 해준다
        filemenu.addAction(self.calibrationFileOpenAction)


        self.splitter = QSplitter(Qt.Horizontal)
        self.tables = QWidget()
        self.tableLayout = QVBoxLayout()
        self.tableLayout.addWidget(self.gaussButton)
        self.tableLayout.addWidget(self.wavelengthPixelTable)
        self.tableLayout.addWidget(self.tableMatchingButtons)
        self.tables.setLayout(self.tableLayout)


        self.splitter.addWidget(self.tables)
        self.spectrums = QWidget()
        self.spectrumsLayout = QVBoxLayout()
        self.spectrumsLayout.addWidget(self.selfSpectrumCanvas)
        self.spectrumsLayout.addWidget(self.standardSpectrumCanvas)
        self.spectrums.setLayout(self.spectrumsLayout)

        self.splitter.addWidget(self.spectrums)

        self.layout.addWidget(self.splitter, 1, 0, 3, 1)
        self.layout.addWidget(self.selfImageCanvas,1,1,1,1)
        self.layout.addWidget(self.selfPeakControl, 2, 1, 1, 1)
        self.layout.addWidget(self.standardSpectrumButton,3,1)


        self.mainWidget.setLayout(self.layout)
        self.setCentralWidget(self.mainWidget)

    '''
    칼리브레이션 파일(comp 파일)을 열고 Identification에 사용될 Y방향(Wavelength에 수직한 방향) 구간을 결정하는 메소드    
    '''
    def onCalibnationFileOpen(self):
        filePath = QFileDialog.getOpenFileName(self, 'Open calibration file','./Spectroscopy_Example/20181023/combine/')[0]
        hdr, data = openFitData(filePath)
        self.selfImageFig.clear()
        self.selfImageAx = self.selfImageFig.add_subplot(111)
        zimshow(self.selfImageAx, data)
        self.selfImageCanvas.draw()

        self.imageWidth = int (data.shape[1])
        self.selfData = data


    '''
    칼리브레이션파일 스펙트럼에서 Peak을 찾는 과정에 관여하는 3가지 initial value(number of peaks, distance btw peaks,
    threshol of peaks)를 조정해서 적절히 Peak을 찾을 수 있게 하는 메소드 
    Slider의 값을 받아서 canvas에 적용한다.
    '''

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

    '''
    칼리브레이션 스펙트럼 그래프를  확대/축소하는 메소드. 
    스크롤을 내리면 마우스 위치를 중심으로 xlim이 4/5배가 되고,
    스크롤을 올리면 마우스 위치를 중심으로 xlim이 5/4배가 된다.
    '''

    def onScrollAtSelfSpectrum(self, event):
        xmin, xmax = self.selfSpectrumAx.get_xlim()
        xnow = event.xdata
        if (event.button == 'up'):
            xsize = int((xmax-xmin)*0.40)
            xmin = xnow-xsize
            xmax = xnow+xsize
            self.selfSpectrumAx.set_xlim(xmin, xmax)
            self.selfSpectrumAx.figure.canvas.draw()
        elif (event.button == 'down'):
            xsize = int((xmax-xmin)*0.625)
            xmin = xnow-xsize
            xmax = xnow+xsize
            self.selfSpectrumAx.set_xlim(xmin, xmax)
            self.selfSpectrumAx.figure.canvas.draw()

    '''
    칼리브레이션 스펙트럼 그래프에서 pickPeak 메소드를 통해 생성된 peakPicker를 움직이고 그 값을 table에 저장하는 메소드
    두 가지 방식으로 peakPicker를 움직일 수 있다. 
    1. Drag and Drop :
        peakPicker를 클릭한 채로 끌어서 움직일 수 있고 마우스를 놓으면 위치가 고정된다. 
        peak 근처 peakDistance 픽셀에서는 자동으로 peak에 붙고 이때 색깔이 연두색으로 바뀐다. 
        마우스가 이동할때 그 픽셀값이 저장된다. 
    2. Double Click :
        selfSpectrum의 Peak의 text를 더블클릭하면 peakPicker 그 text로 이동하고 그 픽셀값이 저장된다. 
    '''

    def onPickPeakAtSelfSpectrum(self, event):
        if self.isMatchFinished: return
        if self.isPicked: return
        if (event.mouseevent.dblclick and event.artist!=self.peakPicker):
            val = round(float(event.artist.get_text()), 4)
            self.wavelengthPixelList.loc[
                self.wavelengthPixelList.Wavelength == self.currentPickedPeakWavelength, 'Pixel'] = val
            self.peakPicker.remove()
            self.peakPicker = self.selfSpectrumAx.axvline(val, color='green', picker=True, pickradius=self.pickDistance)
            self.selfSpectrumAx.figure.canvas.draw()
            self.onChangedList()
            return
        if not event.mouseevent.button == 1 :return
        self.isPicked = True

    def onMoveAtSelfSpectrum(self, event):
        if not event.inaxes: return
        if event.inaxes != self.selfSpectrumAx: return
        if not self.isPicked: return
        self.peakPicker.remove()
        dist = np.min(np.abs(self.selfPeakPixs - event.xdata))
        val = self.selfPeakPixs[np.argmin(np.abs(self.selfPeakPixs - event.xdata))]

        if (dist<self.peakDistance):
            self.peakPicker = self.selfSpectrumAx.axvline(val, color='green', picker=True , pickradius = self.pickDistance)
            self.wavelengthPixelList.loc[self.wavelengthPixelList.Wavelength == self.currentPickedPeakWavelength, 'Pixel'] = val


        else :
            self.peakPicker = self.selfSpectrumAx.axvline(event.xdata, color='blue', picker=True , pickradius = self.pickDistance)
            self.wavelengthPixelList.loc[self.wavelengthPixelList.Wavelength == self.currentPickedPeakWavelength, 'Pixel'] =  int(
                event.xdata)

        self.selfSpectrumAx.figure.canvas.draw()
        self.onChangedList()

    def onReleaseAtSelfSpectrum(self, event):
        if not event.inaxes: return
        if event.inaxes != self.selfSpectrumAx: return
        if not self.isPicked: return
        self.isPicked = False

    '''
    칼리브레이션 스펙트럼 그래프에서 우클릭을 하면 pickPeak과정을 취소하는 메소드
    우클릭을 하면 self.onPickDisable을 호출한다.     
    '''
    def onPressAtSelfSpectrum(self, event):
        if(event.button == 3):
            self.onPickDisable()

    '''
    스탠다드 스펙트럼 그래프를 확대/축소하는 메소드. 
    스크롤을 내리면 마우스 위치를 중심으로 xlim이 4/5배가 되고,
    스크롤을 올리면 마우스 위치를 중심으로 xlim이 5/4배가 된다.
     '''
    def onScrollAtStandardSpectrum(self, event):
        xmin, xmax = self.standardSpectrumAx.get_xlim()
        xnow = event.xdata
        if (event.button == 'up'):
            xsize = int((xmax-xmin)*0.40)
            xmin = xnow-xsize
            xmax = xnow+xsize
            self.standardSpectrumAx.set_xlim(xmin, xmax)
            self.standardSpectrumAx.figure.canvas.draw()
        elif (event.button == 'down'):
            xsize = int((xmax-xmin)*0.625)
            xmin = xnow-xsize
            xmax = xnow+xsize
            self.standardSpectrumAx.set_xlim(xmin, xmax)
            self.standardSpectrumAx.figure.canvas.draw()

    '''
    pickPeak을 시작하기 위한 조건을 나타낸 메소드들. 
    standard spectrum 그래프에서 peak wavelength text를 더블클릭하거나 왼쪽 테이블에서 wavelength를 더블클릭하면 그
    wavelength에 맞는 pickPeak이 실행된다. 
    '''

    def onPickPeakAtStandardSpectrum(self, event):
        if (event.mouseevent.dblclick):
            self.pickPeak(event.artist.get_position()[0], self.standardSpectrum, event.mouseevent)


    def onWavelengthPixelTableDoubleClicked(self, index):
        row = index.row()
        wavelength = self.standardPeakWavelengths[row]
        self.pickPeak(wavelength, self.standardSpectrum)

    '''
    peak wavelength에 맞는 peak Pixel을 찾기 위한 메소드
    선택된 wavelength와 그 peak의 axvline를 파란색으로 강조해서 보여주고 칼리브레이션 스펙트럼 그래프의 중간이나 그래프상의 
    같은 위치에 움직일수 있는 peakPicker를 생성해 해당 wavelength에 맞는 peak Pixel을 찾을 수 있도록 한다. 
    '''

    def pickPeak(self, waveNow, spectrum, mouse=None):
        if not self.isMatchFinished : return
        wavelength = spectrum[0]
        flux = spectrum[1]
        fluxNow = flux[np.where(wavelength==waveNow)][0]
        self.pickedPeak = self.standardSpectrumAx.axvline(waveNow, color='blue')
        if (fluxNow+max(flux)/2.85 > max(flux)):
            self.pickedText = self.standardSpectrumAx.text(waveNow,
                                                           fluxNow+max(flux)/2000,
                                                           waveNow, c='blue',
                                                           bbox=dict(facecolor='white', ec='none'))
        else:
            self.pickedText = self.standardSpectrumAx.text(waveNow,
                                                           fluxNow+max(flux)/2.85,
                                                           waveNow, ha='center', va='center',
                                                           rotation=90, clip_on=True, c='blue',
                                                           bbox=dict(facecolor='white', ec='none'))
        if (mouse is None) :
            xshift = [(self.selfSpectrumAx.get_xlim()[1] - self.selfSpectrumAx.get_xlim()[0]) / 2 + self.selfSpectrumAx.get_xlim()[0] , 0]
        else:
            xshift = self.selfSpectrumAx.transData.inverted().transform((mouse.x, 0))
        self.peakPicker = self.selfSpectrumAx.axvline(xshift[0], color='blue', picker=True , pickradius = self.pickDistance)
        self.selfSpectrumAx.figure.canvas.draw()
        self.standardSpectrumAx.figure.canvas.draw()
        self.currentPickedPeakWavelength = waveNow
        self.isMatchFinished = False

    '''
    pickPeak을 완료하기 위한 메소드 
    매치가 완료되었으면(onMatch) 해당 내용을 리스트에 저장하고 강제종료시(onAbort) 저장하지 않는다. 
    테이블 아래있는 버튼 (Match, Abort)을 클릭하거나 selfSpectrumFigure 위에서 키(m on Match, a on Abort)를 누르면 완료된다. 

    '''

    def keyPressEvent(self, event):
        if self.isMatchFinished: return
        elif event.key() == Qt.Key_M:
            self.onMatch()
        elif event.key() == Qt.Key_A:
            self.onAatch()


    def onMatch(self):
        print('match')
        self.onPickDisable()
        self.onChangedList()

    def onAbort(self):
        self.wavelengthPixelList.loc[self.wavelengthPixelList.Wavelength == self.currentPickedPeakWavelength, 'Pixel'] = 0
        self.onPickDisable()
        self.onChangedList()

    def onExport(self):
        path = QFileDialog.getSaveFileName(self, 'Choose save file location and name ','./', "CSV files (*.csv)")[0]
        self.wavelengthPixelList.to_csv(path, index=False)

    def onImport(self):
        file = QFileDialog.getOpenFileName(self, 'Choose match file','./', "CSV files (*.csv)")[0]
        self.wavelengthPixelList = pd.read_csv(file)
        self.onChangedList()

    def onPressAtStandardSpectrum(self, event):
        if(event.button == 3):
            self.onPickDisable()

    def onPickDisable(self):
        self.pickedPeak.remove()
        self.pickedText.remove()
        self.peakPicker.remove()

        self.selfSpectrumAx.figure.canvas.draw()
        self.standardSpectrumAx.figure.canvas.draw()
        self.isMatchFinished = True

    '''
    onCalibrationOpen 메소드로 열린 comp 이미지 파일에서 사용할 이미지의 y 축 범위를 찾는 메소드.
    마우스 클릭후 끌어서 범위를 결정하면 selfSpectrumDraw에서 
    '''
    def onPressAtImage(self, event):
        if not event.inaxes: return
        if event.inaxes != self.selfImageAx: return
        self.rect = Rectangle((0, 0), 1, 1, alpha=0.5)
        self.selfImageAx.add_patch(self.rect)
        self.x0 = event.xdata
        self.y0 = event.ydata
        self.isPressed = True

    def onMoveAtImage(self, event):
        if not event.inaxes: return
        if event.inaxes != self.selfImageAx: return
        if not self.isPressed: return
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.imageWidth)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((0, self.y0))
        self.selfImageAx.figure.canvas.draw()

    def onReleaseAtImage(self, event):
        if not event.inaxes: return
        if event.inaxes != self.selfImageAx: return
        if not self.isPressed: return
        y = int(self.rect.get_y())
        height = int(self.rect.get_height())
        self.rect.remove()
        self.selfImageAx.figure.canvas.draw()
        if (height<0) :
            height = 0-height
        self.selfImageY = np.array([y, y+height])
        self.selfSpectrumDraw(ymin = y, ymax = y+height, data = self.selfData, args = [self.peakDistance, self.peakThreshold, self.peakNumber])
        self.isPressed = False


    def selfSpectrumDrawWithGauss(self ):
        self.selfSpectrumGaussShow(self.selfSpectrum, self.selfPeakPixs)


    def selfSpectrumDraw(self, ymin, ymax, data, args):

        MINSEP_PK = args[0]  # minimum separation of peaks
        MINAMP_PK = args[1]  # fraction of minimum amplitude (wrt maximum) to regard as peak
        NMAX_PK = args[2]

        self.selfSpectrumFig.clear()
        self.selfSpectrumAx = self.selfSpectrumFig.add_subplot(111)
        identify = np.average(data[ymin:ymax, :], axis=0)
        ground = np.median(identify[0:200])
        max_intens = np.max(identify)
        peakPixs = peak_local_max(identify, indices=True, num_peaks=NMAX_PK,
                                  min_distance=MINSEP_PK,
                                  threshold_abs=max_intens * MINAMP_PK + ground)
        newPeakPixs = []
        for peakPix in peakPixs:
            newPeakPixs.append(peakPix[0])
        peakPixs = newPeakPixs

        self.selfPeakPixs = np.array(peakPixs)
        self.selfSpectrum = np.array(identify)


        for i in peakPixs:
            self.selfSpectrumAx.axvline(i, identify[i] / max(identify) + 0.0003, identify[i] / max(identify) + 0.2, color='c')
            if (identify[i] + max(identify) / 2.85> max(identify)):
                self.selfSpectrumAx.text(i, identify[i] + max(identify) / 2000, str(i), clip_on=False, picker = self.pickDistance)
            else:
                self.selfSpectrumAx.text(i, identify[i] + max(identify) / 2.85, str(i), ha= 'center', va= 'center',
                                         rotation=90, clip_on=True, picker = self.pickDistance)

        self.selfSpectrumAx.plot(identify, color='r')
        self.selfSpectrumAx.set_xlim(0, len(identify))
        self.selfSpectrumAx.set_ylim(0, )
        self.selfSpectrumAx.figure.canvas.draw()



    # 가우스핏을 통해 peak의 정확한 픽셀값을 찾는다.
    # 값이 제일 큰 3개의 peak 스펙트럼과 그 가우스핏을 예시로 보여주고 특히 FWHM값을 모르거나 추측해야 할 경우
    # FWHM값을 변경하면서 가우스핏이 제대로 되었는지

    def selfSpectrumGaussShow(self, identify, peakPixs):
        self.gaussFitWidget.show()
        self.gaussFitWidget.raise_()
        self.selfSpectrumGaussDraw(identify, peakPixs)

#Todo 이거 잘 빼는 방법(바닥값에 맞게 잘 빼는 방법)을 찾아보자.
    def selfSpectrumGaussDraw(self, identify, peakPixs):
        iterations  = 3

        fitter = LevMarLSQFitter()
        self.selfSpectrumGaussFitFig.clear()
        identify = identify - np.median(identify[0:100]) ##여기!
        for i in np.arange(iterations ):
            a = int(iterations / 5)
            if  iterations % 5 != 0 : a = a+1

            ax = self.selfSpectrumGaussFitFig.add_subplot(a, 5, i+1)
            peakPix = peakPixs[-i]
            xs = np.arange(peakPix - int(self.selfFWHM) * 5, peakPix + int(self.selfFWHM) * 5 + 1)


            g_init = Gaussian1D(amplitude=identify[peakPix],
                                mean=peakPix,
                                stddev=self.selfFWHM * gaussian_fwhm_to_sigma,
                                bounds={'amplitude': (0, 2 * identify[peakPix]),
                                        'mean': (peakPix - self.selfFWHM, peakPix + self.selfFWHM),
                                        'stddev': (0, self.selfFWHM)})

            ax.set_ylim(0,max(identify)*1.1)
            fitted = fitter(g_init, xs, identify[xs])
            ax.set_xlim(peakPix - fitted.stddev/gaussian_fwhm_to_sigma*2 , peakPix + fitted.stddev/gaussian_fwhm_to_sigma*2 )
            ax.plot(xs, identify[xs], 'b')
            xss = np.arange(peakPix - self.selfFWHM * 5, peakPix + self.selfFWHM * 5 + 1, 0.01)
            ax.plot(xss, fitted(xss), 'r--')
            ax.figure.canvas.draw()

    def spectrumGaussFit(self, identify, peakPixs):
        #
        return peak_gauss


    def selfSpectrumGaussFit(self, identify, peakPixs):


        self.selfSpectrumFig.clear()
        self.selfSpectrumAx = self.selfSpectrumFig.add_subplot(111)


        fitter = LevMarLSQFitter()

        sortedPeakPixs = np.sort (peakPixs)
        ground = np.median(identify[0:100])
        identify_fit = identify - ground ## 여기도!

        peak_gauss = []
        i = 0
        x_identify = np.arange(len(identify_fit))
        for peakPix in peakPixs:
            g_init = Gaussian1D(amplitude=identify_fit[peakPix],
                                mean=peakPix,
                                stddev=self.selfFWHM * gaussian_fwhm_to_sigma,
                                bounds={'amplitude': (identify_fit[peakPix], 2 * identify_fit[peakPix]),
                                        'mean': (peakPix - self.selfFWHM, peakPix + self.selfFWHM),
                                        'stddev': (0, self.selfFWHM)}
                                )
            fitted = fitter(g_init, x_identify, identify_fit)
            xss = np.arange(peakPix - int(self.selfFWHM) * 3, peakPix + int(self.selfFWHM) * 3 + 1, 0.01)
            self.selfSpectrumAx.plot(xss, fitted(xss) + ground, 'royalblue')

            peak_gauss.append(fitted.mean.value)
            identify_fit = identify_fit - fitted(np.arange(len(identify)))



        '''
        while i < len(sortedPeakPixs)-1:
            peakPix = sortedPeakPixs[i]
            try:
                peakPix2 = sortedPeakPixs[i + 1]
            except:
                peakPix2 = int(peakPix + self.selfFWHM * 10)

            xs = np.arange(peakPix - int(self.selfFWHM) * 3, peakPix2 + int(self.selfFWHM) * 3 + 1)
            g_init = Gaussian1D(amplitude=identify_fit[peakPix],
                                mean=peakPix,
                                stddev=self.selfFWHM * gaussian_fwhm_to_sigma,
                                bounds={'amplitude': (0, 2 * identify_fit[peakPix]),
                                        'mean': (peakPix - self.selfFWHM, peakPix + self.selfFWHM),
                                        'stddev': (0, self.selfFWHM)}
                                )+\
                     Gaussian1D(amplitude=identify_fit[peakPix2],
                                mean=peakPix2,
                                stddev=self.selfFWHM * gaussian_fwhm_to_sigma,
                                bounds={'amplitude': (0, 2 * identify_fit[peakPix2]),
                                        'mean': (peakPix2 - self.selfFWHM, peakPix2 + self.selfFWHM),
                                        'stddev': (0, self.selfFWHM)}
                                )
            fitted = fitter(g_init, xs, identify_fit[xs])



            # fit 한 두 값이 mean 차이가  시그마의 합의 3배 보다 크면 1D fit
            if (fitted.mean_1.value - fitted.mean_0.value > 3* ( fitted.stddev_1 + fitted.stddev_0) ) :
                xs = np.arange(peakPix - int(self.selfFWHM) * 3, peakPix + int(self.selfFWHM) * 3 + 1)
                g_init = Gaussian1D(amplitude=identify_fit[peakPix],
                                    mean=peakPix,
                                    stddev=self.selfFWHM * gaussian_fwhm_to_sigma,
                                    bounds={'amplitude': (0, 2 * identify_fit[peakPix]),
                                            'mean': (peakPix - self.selfFWHM, peakPix + self.selfFWHM),
                                            'stddev': (0, self.selfFWHM)}
                                    )
                fitted = fitter(g_init, xs, identify_fit[xs])
                xss = np.arange(peakPix - int(self.selfFWHM) * 3, peakPix + int(self.selfFWHM) * 3 + 1, 0.01)
                self.selfSpectrumAx.plot(np.arange(len(identify)), fitted(np.arange(len(identify))) + ground, 'royalblue')
                peak_gauss.append(fitted.mean.value)
                i+=1
            else:
                peak_gauss.append(fitted.mean_0.value)
                peak_gauss.append(fitted.mean_1.value)
                xss = np.arange(peakPix - int(self.selfFWHM) * 3, peakPix2 + int(self.selfFWHM) * 3 + 1, 0.01)
                self.selfSpectrumAx.plot(xss, fitted(xss) + ground, 'yellowgreen')
                i += 2

        print(len(peak_gauss))
        print(len(sortedPeakPixs))
        '''
        peak_gauss = np.round(peak_gauss, 4)


        for i, j in zip(peakPixs, peak_gauss):
            self.selfSpectrumAx.axvline(j, identify[i] / max(identify) + 0.0003, identify[i] / max(identify) + 0.2, color='c')
            if (identify[i] + max(identify) / 2.85> max(identify)):
                self.selfSpectrumAx.text(j, identify[i] + max(identify) / 2000, str(j), clip_on=False, picker = self.pickDistance)
            else:
                self.selfSpectrumAx.text(j, identify[i] + max(identify) / 2.85, str(j), ha= 'center', va= 'center',
                                         rotation=90, clip_on=True, picker = self.pickDistance)

        self.selfSpectrumAx.plot(identify, 'r--')
        self.selfSpectrumAx.set_xlim(0, len(identify))
        self.selfSpectrumAx.set_ylim(0, )
        self.selfPeakPixs = np.array(peak_gauss)
        self.selfSpectrumAx.figure.canvas.draw()
        self.gaussFitWidget.close()

    def onFWHMChanged(self, val):
        self.selfFWHM = val/10
        self.FWHMLabel.setText(f'FHWM for comp image = {self.selfFWHM}')
        self.selfSpectrumGaussDraw(self.selfSpectrum, self.selfPeakPixs)

    def onGaussFitButtonClicked(self):
        self.selfSpectrumGaussFit(self.selfSpectrum, self.selfPeakPixs)

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
                self.standardSpectrumAx.text(i, flux[np.where(i == wavelength)][0] + max(flux) / 2000, str(i), clip_on=False, picker = self.pickDistance)
            else:
                self.standardSpectrumAx.text(i, flux[np.where(i == wavelength)][0] + max(flux) / 2.85, str(i),
                         ha='center', va='center', rotation=90, clip_on=True, picker = self.pickDistance)
        self.standardSpectrumAx.plot(wavelength, flux, 'r--')
        self.standardSpectrumAx.set_ylim(0, )
        self.standardSpectrumAx.set_xlim(min(wavelength), max(wavelength))
        self.standardPeakWavelengths = np.array(peaks)
        self.standardSpectrum = data
        self.matchedPeakPixs = np.zeros((self.standardPeakWavelengths.shape[0]))
        matchInfo = np.column_stack((self.standardPeakWavelengths, self.matchedPeakPixs))
        self.wavelengthPixelList = pd.DataFrame(matchInfo,
                                                columns=['Wavelength', 'Pixel'])
        self.onChangedList()
        self.standardSpectrumAx.figure.canvas.draw()

    def onChangedList(self):

        self.wavelengthPixelModel = tableModel(self.wavelengthPixelList)
        self.wavelengthPixelTable.setModel(self.wavelengthPixelModel)


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




matchList = pd.DataFrame( dict(Pixel=[403, 374, 362, 265,
                           245, 238, 223, 213,
                           178, 169, 144, 131, 53],
Wavelength=[8780.6, 8495.4, 8377.6, 7438.9,
                           7245.2, 7173.9, 7032.4, 6929.5,
                           6599.0, 6507.0, 6266.5, 6143.1, 5400.6]))
flux = B[0].data
data = flux
ystep = 3
reId = reIdentifier(matchList, flux)


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

