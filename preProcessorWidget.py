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
from fitImageTableWidget import fitImageTableWidget, tableModel, fileOpener, openFitData
from fitInfo import currentFileInfo
import os


# Todo combine method 선택할수 있게 하기
# Todo 이미지 선택기능(잘 안찍힌 사진들 없앨수 있게)넣기


class preProcessorWidget(QWidget):
    preProcessingDoneSignal = pyqtSignal(currentFileInfo)
    def __init__(self, currentFileInfo):
        super().__init__()
        self.fitImageTableWidget = fitImageTableWidget(currentFileInfo)

        # 현재 작업중인 폴더와 이미지가 있는 기본폴더를 구분해줄 필요가 있으므로 따로 만들어 놓자.
        self.rawFitFileList = currentFileInfo.fitFileList
        self.rawFolderLocation = currentFileInfo.currentFolderLocation
        # 프리프로세싱단계에서 바로 열어줄수 있게 완료된 각각의 sub stage의 filelist를 만들어 놓자.
        self.combFileList = ''
        self.darkFileList = ''
        self.flatFileList = ''

        #combined에서 flag가 바뀌어서 중간에 창을 껐을때 문제가 생길 수 있으므로 완료되지 않고 창을 껐을때 수정할수 있는 단계를 만들자.
        self.isFinished = False
        self.initUI()


    def initUI(self):
        self.fitImageTableWidget.tableEdit()
        self.gridLayout = QGridLayout()

        self.combineBtn = QPushButton('&Combine', self)
        self.biasSubstractionBtn = QPushButton('&Bias Substraction', self)
        self.darkSubstractionBtn = QPushButton('&Dark Substraction', self)
        self.preprocessingBtn = QPushButton('&Preprocessing', self)
        self.doneBtn = QPushButton('&Done', self)

        self.combineBtn.setEnabled(True)
        self.biasSubstractionBtn.setEnabled(False)
        self.darkSubstractionBtn.setEnabled(False)
        self.preprocessingBtn.setEnabled(False)
        self.doneBtn.setEnabled(False)

        self.combineBtn.clicked.connect(self.onCombine)
        self.biasSubstractionBtn.clicked.connect(self.onBiasSubstraction)
        self.darkSubstractionBtn.clicked.connect(self.onDarkSubstraction)
        self.preprocessingBtn.clicked.connect(self.onPreProcessing)
        self.doneBtn.clicked.connect(self.onPreProcessingDone)
        self.pb = QProgressBar(self)




        self.gridLayout.addWidget(self.combineBtn, 0, 0)
        self.gridLayout.addWidget(self.biasSubstractionBtn, 0, 1)
        self.gridLayout.addWidget(self.darkSubstractionBtn, 0, 2)
        self.gridLayout.addWidget(self.preprocessingBtn, 0, 3)
        self.gridLayout.addWidget(self.doneBtn, 0, 4)

        self.gridLayout.addWidget(self.fitImageTableWidget, 1, 0, 1, -1)
        self.gridLayout.addWidget(self.pb, 2, 0, 1, -1)


        self.setLayout(self.gridLayout)

        self.resize(1000, 500)

    # 폴더 구조
    # raw Folder > Preprocessing Folder > Combine, Dark, Flat, P

    def onCombine(self):
        # combine bias, dark(per exposure), flat(per expoure)
        self.fitImageTableWidget.currentFileInfo.currentFolderLocation = self.rawFolderLocation + '/combine'
        combPath = Path(self.fitImageTableWidget.currentFileInfo.currentFolderLocation)
        Path.mkdir(combPath, mode=0o777, exist_ok=True)

        grouped = self.rawFitFileList.groupby(["OBJECT", "EXPTIME"])
        ngroup = grouped.ngroups
        i = 0
        for name, group in grouped:
            # Todo 체크창을 띄워서 OBJECT에 있는것중에 합칠것만 선택할수 있게 하자
            # Todo 보고 관측 이미지 합칠 수 있으면 합치자(따로 이미지 align을 안해서 지금 단계에선 합쳐도 될까 싶다.)
            # 나중에 합칠거 생각해서 일단 /combine에 옮기고 합치는 과정만 나중에 추가해 주자
            if name[0] in ["cali", "flat", "comp_10", "comp_15"]:
                savePath = combPath / f"{name[0]}_{float(name[1]):.1f}.fit"
                ccds = []

                for fName in group['FILE-NAME']:
                    filename = self.rawFolderLocation + '/' + fName
                    hdr, data = openFitData(filename, fitInfo='SNUO')
                    if(self.fitImageTableWidget.currentFileInfo.flags.isCropped):
                        ccds.append(
                            data[self.fitImageTableWidget.currentFileInfo.cropInfo.y0:self.fitImageTableWidget.currentFileInfo.cropInfo.y1,
                            self.fitImageTableWidget.currentFileInfo.cropInfo.x0:self.fitImageTableWidget.currentFileInfo.cropInfo.x1])
                    else:
                        ccds.append(data)
                combined = np.median(ccds, axis=0)
                hdr["NCOMBINE"] = (len(ccds), "Number of images combined")
                combined_ccd = CCDData(data=combined, header=hdr, unit="adu")
                combined_ccd.write(savePath, overwrite=True)

            else:
                j = 0
                for fName in group['FILE-NAME']:
                    j=j+1
                    savePath = combPath / f"{name[0]}_{float(name[1]):.1f}_{j}.fit"
                    filename = self.rawFolderLocation + '/' + fName
                    hdr, data = openFitData(filename, fitInfo='SNUO')
                    if(self.fitImageTableWidget.currentFileInfo.flags.isCropped):
                        data = data[self.fitImageTableWidget.currentFileInfo.cropInfo.y0:self.fitImageTableWidget.currentFileInfo.cropInfo.y1,
                            self.fitImageTableWidget.currentFileInfo.cropInfo.x0:self.fitImageTableWidget.currentFileInfo.cropInfo.x1]
                    ccd = CCDData(data=data, header=hdr, unit="adu")
                    ccd.write(savePath, overwrite=True)

            i = i + 1
            self.step = i / ngroup * 100
            self.onProgressChange()

        files = list(combPath.glob("*.fit"))
        fileInfo = fileOpener(files)
        print(fileInfo)
        fileInfo = np.hstack((fileInfo, np.zeros((fileInfo.shape[0], 1), str)))
        combFileList = pd.DataFrame(np.array(fileInfo),
                                                    columns=['FILE-NAME', 'DATE-OBS', 'EXPTIME', 'IMAGETYPE', 'OBJECT', 'REMARKS'])
        self.fitImageTableWidget.currentFileInfo.fitFileList = combFileList
        self.combFileList = combFileList
        self.fitImageTableWidget.tableEdit()
        self.fitImageTableWidget.currentFileInfo.flags.isReduced = True
        self.combineBtn.setEnabled(False)
        self.biasSubstractionBtn.setEnabled(True)


    #Todo 필터등으로 여러개 있을 상황에 대한 수정
    def onBiasSubstraction(self):
        self.fitImageTableWidget.currentFileInfo.currentFolderLocation = self.rawFolderLocation + '/dark'

        darkPath = Path(self.fitImageTableWidget.currentFileInfo.currentFolderLocation)
        Path.mkdir(darkPath, mode=0o777, exist_ok=True)

        darkFitFileList = self.fitImageTableWidget.currentFileInfo.fitFileList[self.fitImageTableWidget.currentFileInfo.fitFileList['IMAGETYPE'] == 'Dark Frame']
        biasFitFileList =  self.fitImageTableWidget.currentFileInfo.fitFileList[self.fitImageTableWidget.currentFileInfo.fitFileList['IMAGETYPE'] == 'Bias Frame']
        #bias data 불러오기
        biasFileName = self.rawFolderLocation + '/combine/' + biasFitFileList['FILE-NAME'].iloc[0]
        biasHdr, biasData = openFitData(biasFileName)

        #dark data 불러오기 및 빼기
        nlist = len(darkFitFileList)
        i = 0
        for darkFile in darkFitFileList['FILE-NAME']:
            savePath = darkPath / darkFile
            darkFileName = self.rawFolderLocation + '/combine/' + darkFile
            darkHdr, darkData = openFitData(darkFileName)
            data = darkData-biasData
            hdr = darkHdr
            hdr['HISTORY'] = 'BiasSub'
            subbedCCD = CCDData(data=data, header=hdr, unit="adu")
            subbedCCD.write(savePath, overwrite=True)
            self.step = i / nlist * 100
            self.onProgressChange()
        files = list(darkPath.glob("*.fit"))
        fileInfo = fileOpener(files)
        print(fileInfo)
        fileInfo = np.hstack((fileInfo, np.zeros((fileInfo.shape[0], 1), str)))
        darkFileList =  pd.DataFrame(np.array(fileInfo),
                                                    columns=['FILE-NAME', 'DATE-OBS', 'EXPTIME', 'IMAGETYPE', 'OBJECT', 'REMARKS'])
        self.fitImageTableWidget.currentFileInfo.fitFileList = darkFileList
        self.darkFileList = darkFileList
        self.fitImageTableWidget.tableEdit()
        self.biasSubstractionBtn.setEnabled(False)
        self.darkSubstractionBtn.setEnabled(True)
#Todo Flat 여러개 있을때(필터별로)수정
    def onDarkSubstraction(self):
        self.fitImageTableWidget.currentFileInfo.currentFolderLocation = self.rawFolderLocation + '/flat'

        flatPath = Path(self.fitImageTableWidget.currentFileInfo.currentFolderLocation)
        Path.mkdir(flatPath, mode=0o777, exist_ok=True)

        darkFitFileList = self.darkFileList
        flatFitFileList =  self.combFileList[self.combFileList['IMAGETYPE'] == 'Flat Field']
        flatFileName = self.rawFolderLocation + '/combine/' + flatFitFileList['FILE-NAME'].iloc[0]
        flatHdr, flatData = openFitData(flatFileName)
        darkFileName = self.rawFolderLocation + '/dark/' + darkFitFileList['FILE-NAME'][darkFitFileList['EXPTIME']==flatFitFileList['EXPTIME'].iloc[0]].iloc[0]
        darkHdr, darkData = openFitData(darkFileName)

        flatFile  = flatFitFileList['FILE-NAME'].iloc[0]
        savePath = flatPath / 'masterFlat.fit'
        data = flatData - darkData
        hdr = flatHdr
        hdr['HISTORY'] = 'BiasSub'
        subbedCCD = CCDData(data=data, header=hdr, unit="adu")
        subbedCCD.write(savePath, overwrite=True)


        files = list(flatPath.glob("*.fit"))
        fileInfo = fileOpener(files)
        print(fileInfo)
        fileInfo = np.hstack((fileInfo, np.zeros((fileInfo.shape[0], 1), str)))
        flatFileList =  pd.DataFrame(np.array(fileInfo),
                                                    columns=['FILE-NAME', 'DATE-OBS', 'EXPTIME', 'IMAGETYPE', 'OBJECT', 'REMARKS'])
        self.fitImageTableWidget.currentFileInfo.fitFileList = flatFileList
        self.flatFileList = flatFileList
        self.fitImageTableWidget.tableEdit()
        self.darkSubstractionBtn.setEnabled(False)
        self.preprocessingBtn.setEnabled(True)

    def onPreProcessing(self):

        self.fitImageTableWidget.currentFileInfo.currentFolderLocation = self.rawFolderLocation + '/reduced'

        reducedPath = Path(self.fitImageTableWidget.currentFileInfo.currentFolderLocation)
        Path.mkdir(reducedPath, mode=0o777, exist_ok=True)

        darkFitFileList = self.darkFileList
        flatFitFileList =  self.flatFileList



        grouped = self.combFileList.groupby(["OBJECT", "EXPTIME"])
        ngroup = grouped.ngroups
        i = 0
        for name, group in grouped:
            # Todo 체크창을 띄워서 OBJECT에 있는것중에 합칠것만 선택할수 있게 하자
            i = i + 1
            if name[0] not in ["cali", "flat", "comp_10", "comp_15"]:
                for fName in group['FILE-NAME']:

                    saveName = 'r_'+fName
                    savePath = reducedPath / saveName
                    print(self.rawFolderLocation + '/combine/' + fName)

                    objFileName = self.rawFolderLocation + '/combine/' + fName
                    objHdr, objData = openFitData(objFileName)
                    objEXPTIME = objHdr['EXPTIME']
                    darkFileName = self.rawFolderLocation + '/dark/' +  f"cali_{objEXPTIME:.1f}.fit"
                    darkHdr, darkData = openFitData(darkFileName)
                    flatFileName = self.rawFolderLocation + '/flat/' +  "masterFlat.fit"
                    flatHdr, flatData = openFitData(flatFileName)
                    reducedData = (objData-darkData)/flatData
                    reducedHdr = objHdr
                    reducedHdr['HISTORY'] = f"Dark Reduced by :{darkFileName}"
                    reducedHdr['HISTORY'] = f"Flat Reduced by :{flatFileName}"
                    reducedCCD = CCDData(data=reducedData, header=reducedHdr, unit="adu")
                    reducedCCD.write(savePath, overwrite=True)
                    self.step = i / ngroup * 100
                    self.onProgressChange()
#Todo comp image는 dark 안빼줘도 되는건가? flat도 안빼줘도 되는건가?
#일단 comp image는 combine folder에 넣자
        files = list(reducedPath.glob("*.fit"))
        fileInfo = fileOpener(files)
        fileInfo = np.hstack((fileInfo, np.zeros((fileInfo.shape[0], 1), str)))
        reducedFileList = pd.DataFrame(np.array(fileInfo),
                                                    columns=['FILE-NAME', 'DATE-OBS', 'EXPTIME', 'IMAGETYPE', 'OBJECT', 'REMARKS'])
        self.fitImageTableWidget.currentFileInfo.fitFileList = reducedFileList
        self.reducedFileList = reducedFileList
        self.fitImageTableWidget.tableEdit()
        self.isFinished = True
        self.preprocessingBtn.setEnabled(False)
        self.doneBtn.setEnabled(True)

    def onPreProcessingDone(self):
        self.preProcessingDoneSignal.emit(self.fitImageTableWidget.currentFileInfo)

    def onProgressChange(self):
        self.pb.setValue(self.step)

        # 작업 끝난 후에 작업이 완료된 이미지 preprocessing widget의 메인 창에서 볼 수 있게 만들자.
        # 어차피 폴더분류도 해놓은 김에 폴더에서 다시 찾아서 받아오는 방향으로 해보자


    def closeEvent(self, event):
        if(self.isFinished == False):
            reply = QMessageBox.question(self, 'Window Close', 'Preprocessing is not finished. Are you still gonna close the window?',
                                 QMessageBox.Yes | QMessageBox.No, QMessageBox.No)

            if reply == QMessageBox.Yes:
                event.accept()
                self.fitImageTableWidget.currentFileInfo.flags.isReduced = False
            else:
                event.ignore()


