# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 21:44:52 2021

@author: kangGW

@ used Algorithm and some code from ysbach's work
https://github.com/ysBach/SNU_AOclass/blob/master/Notebooks/Spectroscopy_Example.ipynb

지금 나와있는 파이프라인들
aptrace 및 preprocessing
https://aspired.readthedocs.io/en/latest/tutorials/quickstart.html
identification
https://rascal.readthedocs.io/en/latest/installation/installation.html
에 비해 괜찮은가?

"""
import copy
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
from editWidget import editWidget, cropCheckWidget
from fitImageTableWidget import fitImageTableWidget
from preProcessorWidget import preProcessorWidget
from fitInfo import editInfo, currentFileInfo
from identificationWidget import identificationWidget
from apertureTraceWidget import apertureTraceWidget

import os


#설계 : 이미지를 볼 수 있는 창과, 진행 상황을 볼 수 있는 창을 만들자.
#그냥 메인 윈도우 창 그대로 만들어도 될듯
#크게 네 단계로 구분하자. 이미지 합치기, 다크-바이아스 플렛-다크 (이미지-다크)/플렛 
#단계를 어떻게 구분할까?
#combine - bias substraction - dark substraction - preprocessing 
#버튼을 일렬로 배치하고 각각의 단계에서 버튼을 온-오프해서 순서를 정하자
#각각의 단계에서 필요한 이미지의 목록을 - combine에선 콤바인이 완료된 이미지의 목록 -bias substraction에선 bias가 빠진 dark - dark substraction에서는 다크가 빠진 플렛 이미지 - preprocessing 후에는 프리프로세싱이 끝난 이미지들
"""
Schematic
Toolbar section:
File   
┖ Open > opens folder with images
┖ Edit > edit file list #Todo 
┖ Exit > exit program
Processing  
┖ Crop > crops image to target necessary for 
┖ preProcessing > do preprocessing. 
    It's optional. you can upload preprocessed image instead.
Identification
┖identification > do identification and reIdentification. Need to automized.
ApertureTrace
┖apertureTrace > do aperture trace and extract. make it to fits file
standardization
┖standardization > standardize extracted aprture with standard star fits.
Window section:
[fig showing window] | [header table window]
shows fig with configurable pixel value range |
Information including file name, observation date, exposure time, image type, object and remarks.

Interactive functions:
    For main Window:
    Image visulize range is interactively configurable.
    Can show fits image in process by double clicking table.
    For Processing:
    Can crop image into smaller pieces, especially to limit range into actual target data.For identification:
    Can match peaks of comp image pixel and neon ramp wavelength
    Can check reientification result for each x/y cuts.
    extract useless datapoitns and redo fitting [#Todo]
    Choose fit method by button[#Todo]
    For apertureTrace:
    Can check aptrace result for each x/y cuts.
    extract useless datapoitns and redo fitting [#Todo]
    Choose fit method by button[#Todo]
    
    
    """


    #각각의 과정에서 얻어지는 결과물을 기본 패스 아래 다른 폴더에 넣어서 저장하고,
    #각각 파일이 없을 경우 에러를 띄워서 순서를 정하자

#Todo Widget과 Method를 구분해서 GUI 뿐만 아니라, terminal에서도 사용할 수 있게 + 유지보수가 쉽게 만들자.
#Todo 한꺼번에 여러개 aptrace/apextract 가능하게 수정해 보자
#Todo echelle spectrum 에도 적용 할수 있게 만들어 보자


#GUI 디자인


# 가운데에 이미지를 표시하는 큰 창, 오른쪽에 사용하는 파일 목록을 표시하는 작은 창
# 오른쪽 위에는 이미지(fits 파일들), 오른쪽 아래에는 프로세싱 중간에 생성된 그래프나 스펙트럼을 확인할 수 있게 하자.
# 파일 목록을 클릭하면 이미지를 메인창에 표시할 수 있도록 하자.
# 스플리터로 만들어서 이리저리 움직일 수 있도록



# main windows
class MyApp(QMainWindow):
    def __init__(self):
        super().__init__()
        data  = {'FILE-NAME':[' '], 'DATE-OBS':[' '], 'EXPTIME':[' '], 'IMAGETYPE':[' '], 'OBJECT':[' '], 'REMARKS':[' ']}
        self.nullFitFileList = pd.DataFrame(data)
        self.currentFileInfo = currentFileInfo(self.nullFitFileList)



        self.initUI()
        

        
    def initUI(self):
        #open up widgets
        self.identificationWidget = identificationWidget()
        self.identificationWidget.reidentificationWidget.identificationDoneSignal.connect(self.onIdentificationDone)

        self.editWidget = editWidget()

        self.fitImageTableWidget = fitImageTableWidget(self.currentFileInfo)
        self.editWidget.cropCheckWidget.imageCropSignal.connect(self.imageCrop)

        self.resize(1500,750)
        self.center()
#       self.setWindowIcon(QIcon('icon.png'))#아이콘 설정
        self.setCentralWidget(self.fitImageTableWidget)
        #Todo Make a statusbar to show current progress:ex)doing preprocessing now...etc
        
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
        self.IdentificationAction = QAction('Identification', self)
        self.IdentificationAction.setShortcut('Ctrl+I')
        self.IdentificationAction.triggered.connect(self.onIdentification)
        self.IdentificationAction.setStatusTip('match x axis of image and wavelength')

        #Aperture Trace
        self.apertureTraceAction = QAction('apertureTrace', self)
        self.apertureTraceAction.setShortcut('Ctrl+A')
        self.apertureTraceAction.triggered.connect(self.onApertureTrace)
        self.apertureTraceAction.setStatusTip('aperture Trace and Extract Spectrum from Image')

        #standardization
        #Todo Make it!
        self.standardizationAction = QAction('standardization', self)
        self.standardizationAction.setShortcut('Ctrl+s')
        self.standardizationAction.triggered.connect(self.close)
        self.standardizationAction.setStatusTip('aperture Trace and Extract Spectrum from Image')


        
        
        #001 Menubar를 통한 파일 오픈 및 정보 제공
        #파일 오픈은 Open 버튼(shortcut은 ctrl+O)를 통해 이루어지며 이때 폴더를 Open하게 된다.
        #이때 오픈된 파일에 대한 정보는 새 창으로 띄워주며 이 정보는 File_Information에서 다시 열어볼수 있게 하자.
        #완성!! 
        

        

        #002 이미지 크롭 및 전처리
        #Processing tap에 preprocessing 으로 실행

        #완성!!
        
        
        #이후에는 Processing tap에 각각의 과정에 해당하는 버튼(fid, reid, aptrace)을 통해 각 과정을 실행한다.
        #아에 프리프로세싱이랑 spectrum export 부분은 따로 구분해서 다른 방식으로 프리프로세싱한 이미지를 사용해서 작업할 수 있도록 하자.



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
        
        

        processMenu = menubar.addMenu('&Processing')
        processMenu.addAction(self.cropAction)
        processMenu.addAction(self.preProcessingAction)

        identificationMenu = menubar.addMenu('&Identification')
        identificationMenu.addAction(self.IdentificationAction)

        apertureTraceMenu = menubar.addMenu('&ApertureTrace')
        apertureTraceMenu.addAction(self.apertureTraceAction)



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
        

    def onCropAction(self):
        self.editWidget.setFileName(self.fitImageTableWidget.currentFileInfo.currentFileLocation)
        self.editWidget.show()
        self.editWidget.raise_()

    
    
    @pyqtSlot(editInfo)
    def imageCrop(self, edit):
        self.editWidget.close()
        self.fitImageTableWidget.currentFileInfo.editInfo = edit
        self.fitImageTableWidget.currentFileInfo.flags.isEditted = True
        self.fitImageTableWidget.tableEdit()
        
    def onPreprocessing(self):
        currentFileInfo = copy.copy(self.fitImageTableWidget.currentFileInfo)
        self.preProcessorWidget = preProcessorWidget(currentFileInfo)
        self.preProcessorWidget.preProcessingDoneSignal.connect(self.onPreprocessingDone)
        self.preProcessorWidget.show()

    @pyqtSlot(currentFileInfo)
    def onPreprocessingDone(self, currentFileInfo):
        self.fitImageTableWidget.currentFileInfo = currentFileInfo
        self.preProcessorWidget.close()
        self.fitImageTableWidget.tableEdit()

    def onIdentification(self):
        self.identificationWidget.show()
        self.identificationWidget.raise_()

    @pyqtSlot(pd.DataFrame, str)
    def onIdentificationDone(self, regFactor, identificationMethod):
        self.identificationWidget.close()
        self.regFactor = regFactor
        self.identificationMethod = identificationMethod


    def onApertureTrace(self):
        SpectPath = self.fitImageTableWidget.mainFolderLoc  + '/spectrum'
        Path.mkdir(Path(SpectPath), mode=0o777, exist_ok=True)
        # self.apertureTracerWidget = apertureTraceWidget(self.fitImageTableWidget.currentFileInfo.currentFileLocation, regFactor = self.regFactor, identificationMethod= self.identificationMethod, savePath = SpectPath, imageList = self.fitImageTableWidget.currentFileInfo.fitFileList['FILE-NAME'])
        self.apertureTracerWidget = apertureTraceWidget(self.fitImageTableWidget.currentFileInfo.currentFileLocation,
                                                        regFactor=self.regFactor,
                                                        identificationMethod=self.identificationMethod,
                                                        savePath=SpectPath)
        self.apertureTracerWidget.show()
        self.apertureTracerWidget.raise_()

    def onApertureExtract(self):
        print('')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyApp()
    ex.show()
    sys.exit(app.exec_())