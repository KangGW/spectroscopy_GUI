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
import os


def zimshow(ax, image, **kwargs):
    return ax.imshow(image, norm=znorm(image), origin='lower', **kwargs)

def znorm(image):
    return ImageNormalize(image, interval=ZScaleInterval())

class flags():
    def __init__(self):
        self.isCropped = False
        self.isReduced = False

# Todo .ini 나 (세팅파일을 따로 만들거나 main window에서 File>setting 등에서)and/or 읽을 헤더값의 형태를 수정할 수 있게 만들자.
# todo 다른 파일형에 적용할 수 있게 file opener를 수정할 수 있도록 만들자.+하는중

def fileOpener(filelist, fitInfo='SNUO', loadinfo=['DATE-OBS', 'EXPTIME', 'IMAGETYP', 'OBJECT']):
    fileInfo = []
    # Run through all the fits files
    for fitsfile in filelist:
        row = []
        if (fitInfo == 'SNUO'):
            row.append(str(fitsfile).split('\\')[-1])
            hdu = fits.open(fitsfile)
            hdr = hdu[0].header
        # EDIT here for other observatory
        for info in loadinfo:
            row.append(hdr[info])
        fileInfo.append(row)
        hdu.close()

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

    def __init__(self, currentFolderLocation, currentFileLocation, fitFileList, cropInfo, flags, fitInfo ='SNUO'):

        super().__init__()
        self.fitFileList = fitFileList
        self.flags = flags
        self.currentFolderLocation = currentFolderLocation
        self.currentFileLocation = currentFileLocation
        self.cropInfo = cropInfo
        self.fitInfo = fitInfo
        self.initUI()

    def initUI(self):

        self.fig = plt.Figure()
        self.canvas = FigureCanvas(self.fig)
        # Todo canvas에 vmax vmin을 조절해서 이미지를 이쁘게 보일수 있는 바를 만들자!

        # 테이블들
        # fit파일을 보여주는 테이블

        self.fitFileTable = QTableView()
        self.fitFileModel = tableModel(self.fitFileList)
        self.fitFileTable.setModel(self.fitFileModel)
        # 줄(row)별로 선택할수 있게 하는 기능
        self.fitFileTable.setSelectionBehavior(QTableView.SelectRows)
        # 더블클릭하면 선택한 fit파일을 열어주는 기능
        self.fitFileTable.doubleClicked.connect(self.onFitTableDoubleCliked)
        # FitFile및 그래프를 열기 위한 plt canvas
        self.addWidget(self.canvas)
        self.addWidget(self.fitFileTable)

    def fileOpen(self):
        filePath = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        path = Path(filePath)
        fileList = list(path.glob("*.fit"))
        fileInfo = fileOpener(fileList)
        fileInfo = np.hstack((fileInfo, np.zeros((fileInfo.shape[0], 1), str)))
        self.fitFileList = pd.DataFrame(np.array(fileInfo),
                                        columns=['FILE-NAME', 'DATE-OBS', 'EXPTIME', 'IMAGETYPE', 'OBJECT', 'REMARKS'])
        self.currentFolderLocation = filePath

        self.tableEdit()

    def tableEdit(self):
        self.fig.clear()
        self.canvas.draw()
        self.fitFileModel = tableModel(self.fitFileList)
        self.fitFileTable.setModel(self.fitFileModel)

    @pyqtSlot(QModelIndex)
    def onFitTableDoubleCliked(self, index):
        row = index.row()
        file = self.fitFileList['FILE-NAME'][row]
        fileloc = self.currentFolderLocation + '/' + file
        self.currentFileLocation = fileloc

        if (self.flags.isCropped and not self.flags.isReduced):
            data = fits.open(Path(fileloc))[0].data[self.cropInfo.y0:self.cropInfo.y1,
                   self.cropInfo.x0:self.cropInfo.x1]
        else:
            data = fits.open(Path(fileloc))[0].data

        self.fig.clear()
        ax = self.fig.add_subplot(111)
        zimshow(ax, data)
        self.currentFileLocation = fileloc
        self.canvas.draw()



