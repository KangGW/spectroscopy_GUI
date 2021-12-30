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
from fitImageTableWidget import zimshow, znorm, openFitData
import os
from fitInfo import flags, editInfo, currentFileInfo



class editWidget(QWidget):
    '''
    Image Edit Widget to crop, rotate(90 degree, clockwise #Todo) or reflect(horizontally #Todo) image to proper(wavelength short to long, Horizontally) Way.
    rotate and reflect function disabled as of now.
    '''

    cropDoneSignal = pyqtSignal(editInfo)
    """
    This signal is emitted when the widget finishes edition.
    It emits editInfo class that change data while handling.
    """

    def __init__(self, currentFileLocation = ''):
        super().__init__()
        self.filename = currentFileLocation
        self.editInfo = editInfo()
        self.isPressed = False
        self.isCropStarted = False
        self.cropCheckWidget = cropCheckWidget(self.editInfo)
        self.initUI()

    def initUI(self):
        self.layout = QGridLayout()

        self.fig = plt.Figure()
        self.canvas = FigureCanvas(self.fig)
        self.ax = self.fig.add_subplot(111)

        self.canvas.mpl_connect("button_press_event", self.on_press)
        self.canvas.mpl_connect("motion_notify_event", self.on_move)
        self.canvas.mpl_connect("button_release_event", self.on_release)

        #self.cropButton = QPushButton(QIcon('crop.png'), "crop", self)
        #self.rotButton = QPushButton(QIcon('rotation.png'), "rotate", self)
        #self.refButton = QPushButton(QIcon('reflection.png'), "reflect", self)


        # adding action to a button
        #self.cropButton.clicked.connect(self.on_imageCrop)
        #self.rotButton.clicked.connect(self.on_imageRotation)
        #self.refButton.clicked.connect(self.on_imageReflection)


        self.layout.addWidget(self.canvas)
        self.setLayout(self.layout)

        if (self.filename !=''):
            _,self.data = openFitData(self.filename)
            zimshow(self.ax, self.data)
        self.canvas.draw()

    def on_imageCrop(self):
        self.isCropStarted = True

    def setFileName(self, fileName):
        self.filename = fileName
        _,self.data = openFitData(self.filename)
        zimshow(self.ax, self.data)
        self.canvas.draw()

    def on_press(self, event):
        #if not self.isCropStarted : return
        if not event.inaxes : return
        if event.inaxes != self.ax: return
        self.rect = Rectangle((0, 0), 1, 1, alpha=0.5)
        self.ax.add_patch(self.rect)
        self.x0 = event.xdata
        self.y0 = event.ydata
        self.isPressed = True

    def on_move(self, event):
        if not event.inaxes : return
        if event.inaxes != self.ax: return
        if not self.isPressed : return

        self.x1 = event.xdata
        self.y1 = event.ydata
        self.rect.set_width(self.x1 - self.x0)
        self.rect.set_height(self.y1 - self.y0)
        self.rect.set_xy((self.x0, self.y0))
        self.ax.figure.canvas.draw()

    def on_release(self, event):
        if not event.inaxes : return
        if event.inaxes != self.ax: return
        if not self.isPressed: return
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

        self.editInfo.x0 = x0
        self.editInfo.x1 = x1
        self.editInfo.y0 = y0
        self.editInfo.y1 = y1
        self.editInfo.filename = self.filename
        self.cropDoneSignal.emit(self.editInfo)
        self.rect.remove()
        self.ax.figure.canvas.draw()
        self.isPressed = False
        self.isCropStarted = False
        self.cropCheckWidget.setCropInfo(self.editInfo)
        self.cropCheckWidget.show()
        self.cropCheckWidget.raise_()

    def on_imageRotation(self):
        self.editInfo.rotAngle= self.editInfo.rotAngle+90
        if self.editInfo.rotAngle>=360:
            self.editInfo.rotAngle = self.editInfo.rotAngle-360
        self.data = np.array(list(zip(*self.data[::-1])))
        zimshow(self.ax, self.data)
        self.ax.figure.canvas.draw()

    def on_imageReflection(self):
        self.editInfo.reflection= not self.editInfo.reflection


class cropCheckWidget(QWidget):
    '''
    sub-Widget of editWidget
    Ask if crop(#Todo Change to work for all edition, not just for crop) is done with result image. click yes to apply edition, click no to edit more.
    '''

    imageCropSignal = pyqtSignal(editInfo)

    def __init__(self, cropInformation):
        super().__init__()
        self.cropInfo = cropInformation
        self.initUI()

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
        if self.cropInfo.filename == '' : return
        else:
            _,self.data = openFitData(self.cropInfo.filename)
            self.data =self.data[self.cropInfo.y0:self.cropInfo.y1,
               self.cropInfo.x0:self.cropInfo.x1]
            zimshow(self.ax, self.data)
            self.canvas.draw()

    def setCropInfo(self, cropInfo):
        self.cropInfo = cropInfo
        _, self.data = openFitData(self.cropInfo.filename)
        self.data = self.data[self.cropInfo.y0:self.cropInfo.y1,
                    self.cropInfo.x0:self.cropInfo.x1]
        zimshow(self.ax, self.data)
        self.canvas.draw()

    def onNo(self):
        self.close()

    def onYes(self):
        self.imageCropSignal.emit(self.cropInfo)
        self.close()



if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = editWidget("./Spectroscopy_Example/20181023/Lan93101-0004sp.fit")
    ex.show()
    sys.exit(app.exec_())

