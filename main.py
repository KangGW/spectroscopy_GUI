# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 21:44:52 2021

@author: kangGW

@ used Algorithm and some code from ysbach's work
https://github.com/ysBach/SNU_AOclass/blob/master/Notebooks/Spectroscopy_Example.ipynb

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
from astropy.nddata import CCDData
from cropWidget import cropInfo, cropWidget, cropCheckWidget
from fitImageTableWidget import fitImageTableWidget, flags
from preProccesorWidget import preProccesorWidget
import os






#설계 : 이미지를 볼 수 있는 창과, 진행 상황을 볼 수 있는 창을 만들자.
#그냥 메인 윈도우 창 그대로 만들어도 될듯
#크게 네 단계로 구분하자. 이미지 합치기, 다크-바이아스 플렛-다크 (이미지-다크)/플렛 
#단계를 어떻게 구분할까?
#combine - bias substraction - dark substraction - preprocessing 
#버튼을 일렬로 배치하고 각각의 단계에서 버튼을 온-오프해서 순서를 정하자
#각각의 단계에서 필요한 이미지의 목록을 - combine에선 콤바인이 완료된 이미지의 목록 -bias substraction에선 bias가 빠진 dark - dark substraction에서는 다크가 빠진 플렛 이미지 - preprocessing 후에는 프리프로세싱이 끝난 이미지들


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


##Todo Pycharm 문제인지 astropy 문제인지 모르겠는데 fits open에서 파일이 없으면 에러메세지 없이  튕긴다. 이거 에러창 띄우는 방식으로 해결하는 코드 넣어주면 좋을듯



# main windows
class MyApp(QMainWindow):
    imageNameSignal = pyqtSignal(str)
    
    def __init__(self):
        super().__init__() #super는 부모클래스(여기선 QWidget)의 메소드를 불러와주기 위해서 사용
        #fileTable에 사용할 panda 프레임 만들어놓기
        data = {'FILE-NAME':[' '], 'DATE-OBS':[' '], 'EXPTIME':[' '], 'IMAGETYPE':[' '], 'OBJECT':[' '], 'REMARKS':[' ']}
        self.fitFileList = pd.DataFrame(data)

        #Flagges 
        self.flags = flags()
        self.cropInfo = cropInfo()
        #파일 관리를 용이하게 하기 위한 현재 폴더 위치 지정
        self.currentFolderLocation = ''
        #현재 열린 파일 위치
        self.currentFileLocation = ''

        self.initUI()
        
        
        
    def initUI(self):
        self.cropWidget = cropWidget()
        self.fitImageTableWidget = fitImageTableWidget(currentFolderLocation=self.currentFolderLocation, currentFileLocation = self.currentFileLocation, fitFileList = self.fitFileList, cropInfo=self.cropInfo, flags=self.flags)
        self.imageNameSignal.connect(self.cropWidget.onCropStarted)
        self.cropCheckWidget = cropCheckWidget()        
        self.cropWidget.cropDoneSignal.connect(self.cropCheckWidget.onCropDone)
        self.cropCheckWidget.imageCropSignal.connect(self.imageCrop)
        self.resize(1500,750)
        self.center()
#       self.setWindowIcon(QIcon('icon.png'))#아이콘 설정
        self.setCentralWidget(self.fitImageTableWidget)

        
        # 아래에 현재 상황을 표시하는 스테이터스바 하나
        # 툴바도 가능하면 만들자. 
        
        
        
        
        ##상호작용을 통한 작동은 def로 class 아래 구현해서 사용. connect로 연결한다.
        ##QAction도 있다.

        
        #File Actions
        
        #Open
        self.fileOpenAction = QAction('Open', self)
        self.fileOpenAction.setShortcut('Ctrl+O')
        self.fileOpenAction.triggered.connect(self.fitImageTableWidget.fileOpen)
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
        self.preProcessingAction.triggered.connect(self.onPreprocessing)
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
        
        
                

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
    
    def Start(self):
        self.statusBar().showMessage('processing')   
        

    @pyqtSlot()
    def onCropAction(self):
        self.cropWidget.show()
        self.cropWidget.raise_()
        self.imageNameSignal.emit(self.fitImageTableWidget.currentFileLocation)
    
    
    
    @pyqtSlot(cropInfo)     
    def imageCrop(self, crop):
        self.cropWidget.close()
        self.fitImageTableWidget.cropInfo = crop
        self.fitImageTableWidget.flags.isCropped = True
        self.fitImageTableWidget.tableEdit()
        
    def onPreprocessing(self):
        self.preProcessorWidget \
            = preProccesorWidget(currentFolderLocation=self.fitImageTableWidget.currentFolderLocation, currentFileLocation=self.fitImageTableWidget.currentFileLocation,
                                 fitFileList=self.fitImageTableWidget.fitFileList, cropInfo=self.fitImageTableWidget.cropInfo, flags=self.fitImageTableWidget.flags)
        self.preProcessorWidget.show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyApp()
    ex.show()
    sys.exit(app.exec_())