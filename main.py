# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 21:44:52 2021

@author: kangg
@used some of ysbach's work

"""

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

def zimshow(ax, image, **kwargs):
    return ax.imshow(image, norm=znorm(image), origin='lower', **kwargs)

def znorm(image):
    return ImageNormalize(image, interval=ZScaleInterval())

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


#할일: 일단 파일 목록을 받아와서 첫번째 파일에 대한  data = {'FILE-NAME':[' '], 'DATE-OBS':[' '], 'EXPTIME':[' '], 'IMAGETYPE':[' '], 'OBJECT':[' '], 'REMARKS':[' ']}
#목록을 만드는 def를 하나 만들자.

#todo 다른 파일형에 적용할 수 있게 file opener를 수정할 수 있도록 만들자.


def file_opener(filelist):
    fileInfo = []
    # Run through all the fits files
    for fitsfile in filelist:
        row = []
        row.append(str(fitsfile).split('\\')[-1])
        hdu = fits.open(fitsfile)
        hdr = hdu[0].header
        row.append(hdr['DATE-OBS'])
        row.append(hdr['EXPTIME'])
        row.append(hdr['IMAGETYP'])
        row.append(hdr['OBJECT'])
        fileInfo.append(row)
        hdu.close()
    
    return np.array(fileInfo)

class TableModel(QAbstractTableModel):

    def __init__(self, data):
        super(TableModel, self).__init__()
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




class cropInfo:
    def __init__(self):
        self.y0 = 0
        self.y1 = 0
        self.x0 = 0
        self.x1 = 0
        self.filename = ''
    def __repr__(self):
        return str([self.y0,self.y1,self.x0,self.x1])
    

#plt canvas위의 특정 영역을 선택해서 선택 영역을 emit 하는 widget
#ToDo 사각형을 반투명화해서 이쁘게 선택할수 있도록 만들자.

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
    def onCropStarted (self, filename):
        self.filename = filename
        data=fits.open(Path(self.filename))[0].data
        zimshow(self.ax, data)
        self.canvas.draw()
        
    def on_press(self, event):
        print ('press')
        self.rect = Rectangle((0,0), 1, 1)
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
        print ('release')
        x = int(self.rect.get_x())
        y = int(self.rect.get_y())
        width = int(self.rect.get_width())
        height = int(self.rect.get_height())
        
        x0 = x
        x1 = x+width
        y0 = y
        y1 = y+height
        if (x0>x1):
            x0,x1 = x1,x0
        if (y0>y1):
            y0,y1 = y1,y0        
        
        self.cropInfo.x0 = x0
        self.cropInfo.x1 = x1
        self.cropInfo.y0 = y0
        self.cropInfo.y1 = y1
        self.cropInfo.filename = self.filename
        self.cropDoneSignal.emit(self.cropInfo)
        self.rect.remove()
        self.ax.figure.canvas.draw()
        
        
#cropwidget에서 선택된 영역을 crop 할지 물어보고 cropaction을 실행하는 widget
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
        
        self.gridLayout.addWidget(self.canvas,0,0,1,-1)
        self.yesBtn = QPushButton('&Yes', self)
        self.noBtn = QPushButton('&No', self)
        self.gridLayout.addWidget(self.yesBtn,1,0)
        self.gridLayout.addWidget(self.noBtn,1,1)
        self.yesBtn.clicked.connect(self.onYes)
        self.noBtn.clicked.connect(self.onNo)
    
    @pyqtSlot(cropInfo)
    def onCropDone (self, fileinfo):
        self.show()
        self.raise_()
        self.cropInfo = fileinfo
    
        data=fits.open(Path(self.cropInfo.filename))[0].data[self.cropInfo.y0:self.cropInfo.y1,self.cropInfo.x0:self.cropInfo.x1]

        zimshow(self.ax, data)
        self.canvas.draw()
        
    def onNo(self):
        self.close()
    
    def onYes(self):
        self.imageCropSignal.emit(self.cropInfo)
        self.close()

# main windows
class MyApp(QMainWindow):
    imageNameSignal = pyqtSignal(str)
    
    def __init__(self):
        super().__init__() #super는 부모클래스(여기선 QWidget)의 메소드를 불러와주기 위해서 사용
        
        

        #fileTable에 사용할 panda 프레임 만들어놓기
        data = {'FILE-NAME':[' '], 'DATE-OBS':[' '], 'EXPTIME':[' '], 'IMAGETYPE':[' '], 'OBJECT':[' '], 'REMARKS':[' ']}
        self.fitFileList = pd.DataFrame(data)
        data = {'FILE-NAME':[' '], 'REMARKS':[' ']}
        self.graphFileList = pd.DataFrame(data)       
        
        
        
        #Flagges 
        self.isCropped = False
        
        
        self.cropInfo = cropInfo()
        #파일 관리를 용이하게 하기 위한 현재 폴더 위치 지정
        self.currentFolderLocation = ''
        #현재 열린 파일 위치
        self.currentFitFileLocation = ''
        
        
        
        
        
        self.cropWidget = cropWidget()
        self.imageNameSignal.connect(self.cropWidget.onCropStarted)
        self.cropCheckWidget = cropCheckWidget()        
        self.cropWidget.cropDoneSignal.connect(self.cropCheckWidget.onCropDone)
        self.cropCheckWidget.imageCropSignal.connect(self.imageCrop)


        


        
        """self.date = QDate.currentDate()
    
        QToolTip.setFont(QFont('Helvetica', 10))

        """
        #self.move(mv[0], mv[1])
        #self.resize(rs[0], rs[1])
        #self.setGeometry(mv[0], mv[1], rs[0], rs[1])
        self.resize(1000,500)
        self.center()
#        self.setWindowIcon(QIcon('icon.png'))#아이콘 설정
        """
        #button that closes window
        btn = QPushButton('Bye', self)
        btn.move(50,50)
        btn.resize(btn.sizeHint())
        btn.clicked.connect(self.close)
        #btn = 버튼이 cliked = 클릭되었을때 connect = 다음과 
        #같은 역할을 수행한다. self.close = 해당 창을 닫음
        
        #tooltips
        self.setToolTip('<b>Helloooooow</b> there')
        self.setWindowTitle('test')
        btn.setToolTip('Do <b>YOU</b> wanna go away?')
        
        #StatusBar
        self.statusBar().showMessage('idle')
        

        
        #button changes status
        btn2 = QPushButton('Start', self)
        btn2.setFont(QFont('Arial', 30)) 
        btn2.move(300,300)
        btn2.setStyleSheet("color: red;"
                           "border-style: solid;"
                           "border-width: 1px;"
                           "border-color: #FAFAFA;"
                           "border-radius: 3px")
        btn2.clicked.connect(self.Start)
        
        """
        
        #조직도
        #파일
            #오픈(이미지가 들어있는 파일 오픈)
            #에디트(파일 수정(필요없는 파일 제거 등등)
            #exit(끄기)
            
            
            #저장(저장 위치 정하기) _일단 보류
        #프로세싱
            #프리프로세싱(크롭 및 전처리), (아이디피케이션(패턴매칭을 통한 아이디피케이션))
            #리아이디
            #APTRACE
        
        #각각의 과정에서 얻어지는 결과물을 기본 패스 아래 다른 폴더에 넣어서 저장하고, 
        #각각 파일이 없을 경우 에러를 띄워서 순서를 정하자
            
        #
        
        
        #GUI 디자인
        
        
        # 가운데에 이미지를 표시하는 큰 창, 오른쪽에 사용하는 파일 목록을 표시하는 작은 창
        # 오른쪽 위에는 이미지(fits 파일들), 오른쪽 아래에는 프로세싱 중간에 생성된 그래프나 스펙트럼을 확인할 수 있게 하자.
        # 파일 목록을 클릭하면 이미지를 메인창에 표시할 수 있도록 하자.
        # 스플리터로 만들어서 이리저리 움직일 수 있도록
        
        ##스플리터 -> 드래그해서 크기조절 가능!!!!!!
        ##스플리터는 여러개 설정 가능! 스플리터에 다른 스플리터를 addWidget으로 
        ##넣으면 된다
        
        

   
        self.fileSplitter = QSplitter(Qt.Vertical)
        self.mainSplitter = QSplitter(Qt.Horizontal)
        
                
        self.graphFileTable = pandaTableWidget(self.graphFileList)

        self.fig = plt.Figure()
        self.canvas = FigureCanvas(self.fig)
        
        
        self.mainSplitter.addWidget(self.canvas)
        self.mainSplitter.addWidget(self.fileSplitter)
        
        
        #테이블들
        #fit파일을 보여주는 테이블
        
        self.fitFileTable = QTableView()
        self.fitFileModel = TableModel(self.fitFileList)
        self.fitFileTable.setModel(self.fitFileModel)
        #줄(row)별로 선택할수 있게 하는 기능
        self.fitFileTable.setSelectionBehavior(QTableView.SelectRows)
        #더블클릭하면 선택한 fit파일을 열어주는 기능
        self.fitFileTable.doubleClicked.connect(self.onFitTableDoubleCliked)
        #FitFile및 그래프를 열기 위한 plt canvas

        
        
        
        
        
        
        self.fileSplitter.addWidget(self.fitFileTable)
        self.fileSplitter.addWidget(self.graphFileTable)
        
        
        
        
        
        self.setCentralWidget(self.mainSplitter)

        
        # 아래에 현재 상황을 표시하는 스테이터스바 하나
        # 툴바도 가능하면 만들자. 
        
        
        
        
        ##상호작용을 통한 작동은 def로 class 아래 구현해서 사용. connect로 연결한다.
        ##QAction도 있다.
        
        
        
        #File Actions
        
        #Open
        self.fileOpenAction = QAction('Open', self)
        self.fileOpenAction.setShortcut('Ctrl+O')
        self.fileOpenAction.triggered.connect(self.fileOpen)
        self.fileOpenAction.setStatusTip('Open Folder that contains images')
        
        
        #오픈시에 오른쪽에 열린 파일 관련 정보(EXPTIME, OBJECT 등등을 포함한)를 띄우고
        #프리프로세싱 이후에는 라이트프레임만, REID 이후에는 OBJECT 프레임만 남겨서
        #어떤 파일을 지금 사용하고 있는지 확인할 수 있게 하자.
        #관련 위젯을 만들어 띄워줘야할듯
        
        
        
        #Edit
        self.fileEditAction = QAction('Edit', self)
        self.fileEditAction.setShortcut('Ctrl+E')
        self.fileEditAction.triggered.connect(self.close)
        self.fileEditAction.setStatusTip('Edit Images to use')
        
        
        
        #Exit
        self.exitAction = QAction(QIcon('exit.png'), 'Exit', self)
        self.exitAction.setShortcut('Ctrl+Q')
        self.exitAction.setStatusTip('Exit')
        self.exitAction.triggered.connect(self.close)
        
        #Processing Actions
        #Crop
        self.cropAction = QAction('Crop', self)
        self.cropAction.setShortcut('Ctrl+C')
        self.cropAction.triggered.connect(self.onCropAction)
        self.cropAction.setStatusTip('Crop images')       
        #오른쪽에서 선택된 파일을 기반으로 크롭할 수 있게 하자.
        #하나 크롭하면 전부 크롭할 수 있도록.
        #새로운 창을 띄워보자.
        
        
        
        #PreProcessing
        self.preProcessingAction = QAction('Preprocessing', self)
        self.preProcessingAction.setShortcut('Ctrl+P')
        self.preProcessingAction.triggered.connect(self.close)
        self.preProcessingAction.setStatusTip('Preprocess images')
        
        
        
        
        
        
        #Identification
        self.IdentificationAction = QAction('Identificatoin', self)
        self.IdentificationAction.setShortcut('Ctrl+I')
        self.IdentificationAction.triggered.connect(self.close)
        self.IdentificationAction.setStatusTip('match x axis of image and wavelength')

        #ReIdentification
        self.reIdentificationAction = QAction('Reidentificatoin', self)
        self.reIdentificationAction.setShortcut('Ctrl+R')
        self.reIdentificationAction.triggered.connect(self.close)
        self.reIdentificationAction.setStatusTip('match x/y axis of image and wavelength')        
                
        #Aperture Trace
        
        
        
        #001 Menubar를 통한 파일 오픈 및 정보 제공
        #파일 오픈은 Open 버튼(shortcut은 ctrl+O)를 통해 이루어지며 이때 폴더를 Open하게 된다.
        #이때 오픈된 파일에 대한 정보는 새 창으로 띄워주며 이 정보는 File_Information에서 다시 열어볼수 있게 하자.
        #완성!! 
        
        #추가기능으로 Open 버튼을 통해 연 파일 목록중 사용할 것과 사용하지 않을 것을 분류할 수 있도록 하는 기능을 넣어보자
        
        
        
        #전처리부터 하자.
        
    
        

        
        #002 이미지 크롭 및 전처리
        #Processing tap에 preprocessing 으로 실행
        
        
        
        #이후에는 Processing tap에 각각의 과정에 해당하는 버튼(fid, reid, aptrace)을 통해 각 과정을 실행한다.
        #interactive gui를 만들 부분은 처음에 스펙트럼 자르는 부분(전처리과정)-> 이미지 크롭툴
        #first identification에서 패턴매칭하는거(가능하면 처음에 비슷한부분 자동으로 맞춰주기)
        
        #aptrace에서 sky랑 ap 부분 맞춰주는거
        
        
        
        
        #메뉴바 만들기         
        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)
        filemenu = menubar.addMenu('&File')#&는 File을 Alt F로 실행하게 해준다
        filemenu.addAction(self.fileOpenAction)
        filemenu.addAction(self.exitAction)
        
        
        
        
        
#        startAction = QAction('Start', self)
#        startAction.setShortcut('Shift+Space')
#        startAction.setStatusTip('Start')
#        startAction.triggered.connect(self.Start)
        
#        filemenu.addAction(startAction)
        
        
        
        
        filemenu = menubar.addMenu('&Processing')
        filemenu.addAction(self.cropAction)
        filemenu.addAction(self.preProcessingAction)
        filemenu.addAction(self.IdentificationAction)
        
        #ToolBar
        #self.toolbar = self.addToolBar('Exit')
        #self.toolbar.addAction(exitAction)
        
        
        #date and time
        #self.statusBar().showMessage(self.date.toString(Qt.DefaultLocaleLongDate))
        
        
                
    def fileOpen(self):
        file = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        toppath = Path(file)
        files = list(toppath.glob("*.fit"))
        fileInfo = file_opener(files)
        fileInfo = np.hstack((fileInfo, np.zeros((fileInfo.shape[0], 1), str)))
        
        self.fitFileList = pd.DataFrame(np.array(fileInfo),
                   columns=['FILE-NAME', 'DATE-OBS', 'EXPTIME', 'IMAGETYPE', 'OBJECT', 'REMARKS'])
        self.currentFolderLocation = file
        self.onChangedFileList()
        
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
    
    def Start(self):
        self.statusBar().showMessage('processing')   
        
    def onChangedFileList(self):
        self.fitFileModel = TableModel(self.fitFileList)
        self.fitFileTable.setModel(self.fitFileModel)
    
    @pyqtSlot(QModelIndex)
    def onFitTableDoubleCliked(self, index):
        row = index.row()
        file = self.fitFileList['FILE-NAME'][row]
        fileloc = self.currentFolderLocation+'/'+file
        self.currentFitFileLocation = fileloc
        if (self.isCropped):
            data=fits.open(Path(fileloc))[0].data[self.cropInfo.y0:self.cropInfo.y1,self.cropInfo.x0:self.cropInfo.x1]
        else:
            data=fits.open(Path(fileloc))[0].data
        ax = self.fig.add_subplot(111)
        zimshow(ax, data)
        self.canvas.draw()
        
    @pyqtSlot()
    def onCropAction(self):
        self.cropWidget.show()
        self.cropWidget.raise_()
        self.imageNameSignal.emit(self.currentFitFileLocation)
    
    
    
    @pyqtSlot(cropInfo)     
    def imageCrop(self, crop):
        self.cropWidget.close()
        self.cropInfo = crop
        self.isCropped = True

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyApp()
    ex.show()
    sys.exit(app.exec_())