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
from fitImageTableWidget import fitImageTableWidget, tableModel, fileOpener
import os


# Todo combine method 선택할수 있게 하기
# Todo 이미지 선택기능(잘 안찍힌 사진들 없앨수 있게)넣기


class preProccesorWidget(QWidget):
    def __init__(self, currentFolderLocation, currentFileLocation, fitFileList, cropInfo, flags):
        super().__init__()
        self.fitImageTableWidget = fitImageTableWidget(currentFolderLocation=currentFolderLocation, currentFileLocation = currentFileLocation, fitFileList = fitFileList, cropInfo=cropInfo, flags=flags)
        self.initUI()
        # 현재 작업중인 폴더와 이미지가 있는 기본폴더를 구분해줄 필요가 있으므로 따로 만들어 놓자.
        self.rawFitFileList = fitFileList
        self.rawFolderLocation = currentFolderLocation

    def initUI(self):
        self.fitImageTableWidget.tableEdit()
        self.gridLayout = QGridLayout()

        self.combineBtn = QPushButton('&Combine', self)
        self.biasSubstractionBtn = QPushButton('&Bias Substraction', self)
        self.darkSubstractionBtn = QPushButton('&Dark Substraction', self)
        self.preprocessingBtn = QPushButton('&Preprocessing', self)
        self.combineBtn.clicked.connect(self.onCombine)
        self.pb = QProgressBar(self)

        self.gridLayout.addWidget(self.combineBtn, 0, 0)
        self.gridLayout.addWidget(self.biasSubstractionBtn, 0, 1)
        self.gridLayout.addWidget(self.darkSubstractionBtn, 0, 2)
        self.gridLayout.addWidget(self.preprocessingBtn, 0, 3)

        self.gridLayout.addWidget(self.fitImageTableWidget, 1, 0, 1, -1)
        self.gridLayout.addWidget(self.pb, 2, 0, 1, -1)

        # self.vbox.addWidget(self.imageWidget)

        self.setLayout(self.gridLayout)

        self.resize(1000, 500)

    # 폴더 구조
    # raw Folder > Preprocessing Folder > Combine, Dark, Flat, P

    def onCombine(self):
        # combine bias, dark(per exposure), flat(per expoure)
        self.fitImageTableWidget.currentFolderLocation = self.rawFolderLocation + '/combine'
        combPath = Path(self.fitImageTableWidget.currentFolderLocation)
        Path.mkdir(combPath, mode=0o777, exist_ok=True)

        grouped = self.rawFitFileList.groupby(["OBJECT", "EXPTIME"])
        ngroup = grouped.ngroups
        i = 0
        for name, group in grouped:
            # Todo 체크창을 띄워서 OBJECT에 있는것중에 합칠것만 선택할수 있게 하자

            if name[0] in ["cali", "flat", "comp_10", "comp_15"]:
                savePath = combPath / f"{name[0]}_{float(name[1]):.1f}.fit"
                ccds = []

                for fpath in group['FILE-NAME']:
                    ccd = fits.open(Path(self.rawFolderLocation + '/' + fpath))
                    if(self.fitImageTableWidget.flags.isCropped):
                        ccds.append(
                            ccd[0].data[self.fitImageTableWidget.cropInfo.y0:self.fitImageTableWidget.cropInfo.y1,
                            self.fitImageTableWidget.cropInfo.x0:self.fitImageTableWidget.cropInfo.x1])
                    else:
                        ccds.append(ccd[0].data)
                combined = np.median(ccds, axis=0)
                hdr = ccd[0].header  # an arbitrary header: the last of the combined ccds
                hdr["NCOMBINE"] = (len(ccds), "Number of images combined")
                # print(combined)
                print(combined)
                combined_ccd = CCDData(data=combined, header=hdr, unit="adu")
                combined_ccd.write(savePath, overwrite=True)
            i = i + 1
            self.step = i / ngroup * 100
            self.onProgressChange()

        files = list(combPath.glob("*.fit"))
        fileInfo = fileOpener(files)
        print(fileInfo)
        fileInfo = np.hstack((fileInfo, np.zeros((fileInfo.shape[0], 1), str)))

        self.fitImageTableWidget.fitFileList = pd.DataFrame(np.array(fileInfo),
                                                    columns=['FILE-NAME', 'DATE-OBS', 'EXPTIME', 'IMAGETYPE', 'OBJECT', 'REMARKS'])
        self.fitImageTableWidget.tableEdit()
        self.fitImageTableWidget.flags.isReduced = True
    def onProgressChange(self):
        self.pb.setValue(self.step)

        # 작업 끝난 후에 작업이 완료된 이미지 preprocessing widget의 메인 창에서 볼 수 있게 만들자.
        # 어차피 폴더분류도 해놓은 김에 폴더에서 다시 찾아서 받아오는 방향으로 해보자


'''        
    def onBiasSubstraction:

    def onDarkSubstraction:

    def onPreprocessing:

'''
