import sys
from PyQt5.QtWidgets import *
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
from fitInfo import flags, cropInfo, currentFileInfo
import os







def zimshow(ax, image, **kwargs):
    return ax.imshow(image, norm=znorm(image), origin='lower', **kwargs)

def znorm(image):
    return ImageNormalize(image, interval=ZScaleInterval())



# Todo .ini 나 (세팅파일을 따로 만들거나 main window에서 File>setting 등에서)and/or 읽을 헤더값의 형태를 수정할 수 있게 만들자.
# todo 다른 파일형에 적용할 수 있게 file opener를 수정할 수 있도록 만들자.+하는중
# Todo 최대한 파일 여는거 openFitData로 처리해서 오류처리랑 다른 파일형으로 수정할수 있게 만들자.
# Todo 에러처리 이쁘게(에러나면 추가로 창을 띄워서 처리할 수 있게)
# Todo fileNotFoundError에서 창을 하나 떠 띄워서 원하는 파일을 찾을 수 있게 만들자. 지금은 list로 통째로 받는게 많아서 (과정이)너무 길어져서 비효율적
class fileNotFoundError(Exception):
    def __init__(self):
        super().__init__('file Not Found')


def openFitData(filename, fitInfo='default'):
    try:
        hdu = fits.open(filename)
        #중간에 만든 파일들(combine이나 preprocessing으로) 열때
        if (fitInfo=='default'):
            hdr = hdu[0].header
            data = hdu[0].data
        #for SNU 1m telescope 1D spectrometer
        elif (fitInfo =='SNUO'):
            hdr = hdu[0].header
            data = hdu[0].data
        return hdr, data

    except:
        print(f'could not open {filename}, try other filename')
        raise fileNotFoundError


def fileOpener(filelist, loadinfo=['DATE-OBS', 'EXPTIME', 'IMAGETYP', 'OBJECT']):
    fileInfo = []
    # Run through all the fits files
    for fitsfile in filelist:
        row = []
        row.append(str(fitsfile).split('\\')[-1])
        hdr, data= openFitData(fitsfile)
        # EDIT here for other observatory
        for info in loadinfo:
            row.append(hdr[info])
        fileInfo.append(row)

    return np.array(fileInfo)


class tableModel(QAbstractTableModel):

    def __init__(self, data):
        super(tableModel, self).__init__()
        self._data = data

    def data(self, index, role):
        if role == Qt.DisplayRole:
            value = self._data.iloc[index.row(), index.column()]
            return str(value)

    def rowCount(self, index):
        return self._data.shape[0]

    def columnCount(self, index):
        return self._data.shape[1]

    def headerData(self, section, orientation, role):
        # section is the index of the column/row.
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return str(self._data.columns[section])

            if orientation == Qt.Vertical:
                return str(self._data.index[section])

class pandaTableWidget(QTableWidget):
    def __init__(self, df, parent=None):
        QTableWidget.__init__(self, parent)
        self.df = df
        r = len(self.df.index)
        c = len(self.df.columns)
        self.setRowCount(r)
        self.setColumnCount(c)

        for i in range(self.rowCount()):
            for j in range(self.columnCount()):
                x = self.df.iloc[i, j]
                self.setItem(i, j, QTableWidgetItem(x))
        self.setHorizontalHeaderLabels(self.df.columns)

#아예 fitimage-table 형태로 볼 수 있는 위젯을 하나 만드는게 편할듯 > 만듬
#Todo 다른 형식의 fit 파일도 열수 있도록 수정.

class fitImageTableWidget(QSplitter):

    def __init__(self, currentFileInfo):

        super().__init__()
        self.currentFileInfo = currentFileInfo
        self.currentData = ''
        self.initUI()

    def initUI(self):

        self.fig = plt.Figure()
        self.canvas = FigureCanvas(self.fig)
        # Todo canvas에 vmax vmin을 조절해서 이미지를 이쁘게 보일수 있는 바를 만들자!

        # 테이블들
        # fit파일을 보여주는 테이블
        self.imageLayout = QVBoxLayout()
        self.fitFileTable = QTableView()
        self.fitFileModel = tableModel(self.currentFileInfo.fitFileList)
        self.fitFileTable.setModel(self.fitFileModel)
        # 줄(row)별로 선택할수 있게 하는 기능
        self.fitFileTable.setSelectionBehavior(QTableView.SelectRows)
        # 더블클릭하면 선택한 fit파일을 열어주는 기능
        self.fitFileTable.doubleClicked.connect(self.onFitTableDoubleCliked)
        # FitFile및 그래프를 열기 위한 plt canvas
        self.imageVRangeSlider = imageVRangeSlider()
        self.imageVRangeSlider.vChangeSignal.connect(self.onVChanged)


        self.imageLayout.addWidget(self.canvas)
        self.imageLayout.addWidget(self.imageVRangeSlider)
        self.imageWidget = QWidget()
        self.imageWidget.setLayout(self.imageLayout)
        self.addWidget(self.imageWidget)
        self.addWidget(self.fitFileTable)

    def fileOpen(self):
        filePath = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        path = Path(filePath)
        fileList = list(path.glob("*.fit"))

        if (len(fileList)==0):
            filePath, fileList = self.fitFileErrorOnFileOpen()
        fileInfo = fileOpener(fileList)
        fileInfo = np.hstack((fileInfo, np.zeros((fileInfo.shape[0], 1), str)))
        self.currentFileInfo.fitFileList = pd.DataFrame(np.array(fileInfo),
                                        columns=['FILE-NAME', 'DATE-OBS', 'EXPTIME', 'IMAGETYPE', 'OBJECT', 'REMARKS'])
        self.currentFileInfo.currentFolderLocation = filePath
        self.tableEdit()

    def fitFileErrorOnFileOpen(self):
        filePath = str(QFileDialog.getExistingDirectory(self, "No Fit file on directory, please select another directory"))
        path = Path(filePath)
        fileList = list(path.glob("*.fit"))
        if (len(fileList) == 0):
            self.fitFileErrorOnFileOpen()
            return 0
        return filePath, fileList


    def tableEdit(self):
        self.fig.clear()
        self.canvas.draw()
        self.fitFileModel = tableModel(self.currentFileInfo.fitFileList)
        self.fitFileTable.setModel(self.fitFileModel)

    @pyqtSlot(QModelIndex)
    def onFitTableDoubleCliked(self, index):
        row = index.row()
        file = self.currentFileInfo.fitFileList['FILE-NAME'][row]
        fileloc = self.currentFileInfo.currentFolderLocation + '/' + file
        self.currentFileInfo.currentFileLocation = fileloc
        hdr, data = openFitData(fileloc, fitInfo='SNUO')


        if (self.currentFileInfo.flags.isCropped and not self.currentFileInfo.flags.isReduced):
            data = data[self.currentFileInfo.cropInfo.y0:self.currentFileInfo.cropInfo.y1,
                   self.currentFileInfo.cropInfo.x0:self.currentFileInfo.cropInfo.x1]

        self.fig.clear()
        self.currentData = data
        ax = self.fig.add_subplot(111)
        zimshow(ax, data)

        self.imageVRangeSlider.setRangeLimit(0,data.max())
        self.imageVRangeSlider.setRange(data.min(),data.max())
        self.imageVRangeSlider.setTickInterval(int(data.max()/100))
        self.imageVRangeSlider.update()
        self.currentFileInfo.currentFileLocation = fileloc
        self.canvas.draw()

    @pyqtSlot(tuple)
    def onVChanged(self, vValue):
        VMIN, VMAX = vValue
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        zimshow(ax, self.currentData, vmin=VMIN, vmax=VMAX)
        self.canvas.draw()

class imageVRangeSlider(QWidget):
    vChangeSignal = pyqtSignal(tuple)
    def __init__(self):
        super().__init__()

        self.VMIN = 0
        self.VMAX = 100
        #styleoptionslide를 사용해서 custom slider를 만든다.
        self.opt = QStyleOptionSlider()
        self.opt.minimum = 0
        self.opt.maximum = 100
        #슬라이더가 틱의 어디에 있을 것인가.
        self.setTickPosition(QSlider.TicksAbove)
        self.setTickInterval(1)

        self.setSizePolicy(
            QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed, QSizePolicy.Slider)
        )

    def setRangeLimit(self, minPixelValue, maxPixelValue):
        self.opt.minimum = int(minPixelValue)
        self.opt.maximum = int(maxPixelValue)

    def setRange(self, VMIN, VMAX):
        self.VMIN = int(VMIN)
        self.VMAX = int(VMAX)

    def getRange(self):
        return self.VMIN, self.VMAX

    def setTickPosition(self, position: QSlider.TickPosition):
        self.opt.tickPosition = position

    def setTickInterval(self, ti: int):
        self.opt.tickInterval = ti
    #Todo viridis를 rect에 넣어보자
    #Todo 각 slider에 label이나 값을 붙이자
    #Todo tick에 숫자를 넣어보자
    #여기서 직접 그린다.
    def paintEvent(self, event: QPaintEvent):
        #우선 그릴 페인터를 호출
        painter = QPainter(self)
        #opt 내장함수 불러오기
        self.opt.initFrom(self)
        #중간에 사각형 이거 수정해서 이미지 모양으로 만들자.
        self.opt.rect = self.rect()
        #현재 슬라이더 위치
        self.opt.sliderPosition = 0
        #서브컨트롤, range를 설정하기 위해 2개의 슬라이더를 움직여야되므로 어떤 슬라이더가 움직일지 결정하는 SC slider
        #Groove subcontrol과 슬라이더의 위치를 결정하는 틱마크 서브컨트롤을 사용한다.
        self.opt.subControls = QStyle.SC_SliderGroove | QStyle.SC_SliderTickmarks

        #서브컨트롤을 사용해야하므로 이를 받을수있는 컴플렉스컨트롤 슬라이더를 사용한다.
        self.style().drawComplexControl(QStyle.CC_Slider, self.opt, painter)

        #  Draw INTERVAL

        color = self.palette().color(QPalette.Highlight)
        color.setAlpha(160)
        painter.setBrush(QBrush(color))
        painter.setPen(Qt.NoPen)

        #이건 중간에 선택된 range를 볼수 있게 하는 사각형을 놓는 용도
        self.opt.sliderPosition = self.VMIN
        x_left_handle = (
            self.style()
            .subControlRect(QStyle.CC_Slider, self.opt, QStyle.SC_SliderHandle)
            .right()
        )

        self.opt.sliderPosition = self.VMAX
        x_right_handle = (
            self.style()
            .subControlRect(QStyle.CC_Slider, self.opt, QStyle.SC_SliderHandle)
            .left()
        )

        groove_rect = self.style().subControlRect(
            QStyle.CC_Slider, self.opt, QStyle.SC_SliderGroove
        )

        selection = QRect(
            x_left_handle,
            groove_rect.y(),
            x_right_handle - x_left_handle,
            groove_rect.height(),
        ).adjusted(-1, 1, 1, -1)

        painter.drawRect(selection)

        # Draw first handle

        self.opt.subControls = QStyle.SC_SliderHandle
        self.opt.sliderPosition = self.VMIN
        self.style().drawComplexControl(QStyle.CC_Slider, self.opt, painter)

        # Draw second handle
        self.opt.sliderPosition = self.VMAX
        self.style().drawComplexControl(QStyle.CC_Slider, self.opt, painter)

    def mousePressEvent(self, event: QMouseEvent):

        self.opt.sliderPosition = self.VMIN
        self._first_sc = self.style().hitTestComplexControl(
            QStyle.CC_Slider, self.opt, event.pos(), self
        )

        self.opt.sliderPosition = self.VMAX
        self._second_sc = self.style().hitTestComplexControl(
            QStyle.CC_Slider, self.opt, event.pos(), self
        )

    def mouseMoveEvent(self, event: QMouseEvent):

        distance = self.opt.maximum - self.opt.minimum

        pos = self.style().sliderValueFromPosition(
            0, distance, event.pos().x(), self.rect().width()
        )

        if self._first_sc == QStyle.SC_SliderHandle:
            if pos <= self.VMAX:
                self.VMIN = pos
                self.update()
                return

        if self._second_sc == QStyle.SC_SliderHandle:
            if pos >= self.VMIN:
                self.VMAX = pos
                self.update()

    def mouseReleaseEvent(self, event: QMouseEvent):
        vmin, vmax = self.getRange()
        print(self.opt.tickInterval)
        print(vmin, vmax)
        self.vChangeSignal.emit((vmin,vmax))

    def sizeHint(self):
        """ override """
        SliderLength = 1000
        TickSpace = 5

        w = SliderLength
        h = self.style().pixelMetric(QStyle.PM_SliderThickness, self.opt, self)

        if (
            self.opt.tickPosition & QSlider.TicksAbove
            or self.opt.tickPosition & QSlider.TicksBelow
        ):
            h += TickSpace

        return (
            self.style()
            .sizeFromContents(QStyle.CT_Slider, self.opt, QSize(w, h), self)
            .expandedTo(QApplication.globalStrut())
        )
