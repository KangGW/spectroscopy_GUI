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
from fitImageTableWidget import zimshow, znorm
import os
from fitInfo import flags, cropInfo, currentFileInfo
# plt canvas위의 특정 영역을 선택해서 선택 영역을 emit 하는 widget



class cropWidget(QWidget):
    cropDoneSignal = pyqtSignal(cropInfo)

    def __init__(self):
        super().__init__()
        self.initUI()
        self.cropInfo = cropInfo()

    def initUI(self):

        self.hbox = QHBoxLayout()

        self.fig = plt.Figure()
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(111)

        self.canvas.mpl_connect("button_press_event", self.on_press)
        self.canvas.mpl_connect("motion_notify_event", self.on_move)
        self.canvas.mpl_connect("button_release_event", self.on_release)

        self.hbox.addWidget(self.canvas)
        self.setLayout(self.hbox)

    @pyqtSlot(str)
    def onCropStarted(self, currenFileLocation):
        self.filename = currenFileLocation
        data = fits.open(Path(self.filename))[0].data
        zimshow(self.ax, data)
        self.canvas.draw()

    def on_press(self, event):
        print('press')
        self.rect = Rectangle((0, 0), 1, 1, alpha=0.5)
        self.ax.add_patch(self.rect)
        self.x0 = event.xdata
        self.y0 = event.ydata

    def on_move(self, event):
        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.ax.figure.canvas.draw()

    def on_release(self, event):
        print('release')
        x = int(self.rect.get_x())
        y = int(self.rect.get_y())
        width = int(self.rect.get_width())
        height = int(self.rect.get_height())

        x0 = x
        x1 = x + width
        y0 = y
        y1 = y + height
        if (x0 > x1):
            x0, x1 = x1, x0
        if (y0 > y1):
            y0, y1 = y1, y0

        self.cropInfo.x0 = x0
        self.cropInfo.x1 = x1
        self.cropInfo.y0 = y0
        self.cropInfo.y1 = y1
        self.cropInfo.filename = self.filename
        self.cropDoneSignal.emit(self.cropInfo)
        self.rect.remove()
        self.ax.figure.canvas.draw()




# cropwidget에서 선택된 영역을 crop 할지 물어보고 cropaction을 실행하는 widget
class cropCheckWidget(QWidget):
    imageCropSignal = pyqtSignal(cropInfo)

    def __init__(self):
        super().__init__()
        self.initUI()
        self.cropInfo = cropInfo()

    def initUI(self):
        self.gridLayout = QGridLayout()
        self.setLayout(self.gridLayout)
        self.fig = plt.Figure()
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(111)

        self.gridLayout.addWidget(self.canvas, 0, 0, 1, -1)
        self.yesBtn = QPushButton('&Yes', self)
        self.noBtn = QPushButton('&No', self)
        self.gridLayout.addWidget(self.yesBtn, 1, 0)
        self.gridLayout.addWidget(self.noBtn, 1, 1)
        self.yesBtn.clicked.connect(self.onYes)
        self.noBtn.clicked.connect(self.onNo)

    @pyqtSlot(cropInfo)
    def onCropDone(self, fileinfo):
        self.show()
        self.raise_()
        self.cropInfo = fileinfo

        data = fits.open(Path(self.cropInfo.filename))[0].data[self.cropInfo.y0:self.cropInfo.y1,
               self.cropInfo.x0:self.cropInfo.x1]

        zimshow(self.ax, data)
        self.canvas.draw()

    def onNo(self):
        self.close()

    def onYes(self):
        self.imageCropSignal.emit(self.cropInfo)
        self.close()


