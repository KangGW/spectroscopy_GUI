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

# Work by ysbach(2019) `https://github.com/ysBach/SNU_AOclass/blob/master/Notebooks/Spectroscopy_Example.ipynb
def make_summary(filelist, extension=0, fname_option='relative',
                 output=None, format='ascii.csv',
                 keywords=[], dtypes=[],
                 example_header=None, sort_by='file', verbose=True):
    """ Extracts summary from the headers of FITS files.
    Parameters
    ----------
    filelist: list of str (path-like)
        The list of file paths relative to the current working directory.

    extension: int or str
        The extension to be summarized.

    fname_option: str {'absolute', 'relative', 'name'}
        Whether to save full absolute/relative path or only the filename.

    output: str or path-like
        The directory and file name of the output summary file.

    format: str
        The astropy.table.Table output format.

    keywords: list
        The list of the keywords to extract (keywords should be in ``str``).

    dtypes: list
        The list of dtypes of keywords if you want to specify. If ``[]``,
        ``['U80'] * len(keywords)`` will be used. Otherwise, it should have
        the same length with ``keywords``.

    example_header: str or path-like
        The path including the filename of the output summary text file.

    sort_by: str
        The column name to sort the results. It can be any element of
        ``keywords`` or ``'file'``, which sorts the table by the file name.
    """

    if len(filelist) == 0:
        print("No FITS file found.")
        return

    def _get_fname(path):
        if fname_option == 'relative':
            return str(path)
        elif fname_option == 'absolute':
            return str(path.absolute())
        else:
            return path.name

    options = ['absolute', 'relative', 'name']
    if fname_option not in options:
        raise KeyError(f"fname_option must be one of {options}.")

    skip_keys = ['COMMENT', 'HISTORY']

    if verbose:
        if (keywords != []) and (keywords != '*'):
            print("Extracting keys: ", keywords)
        str_example_hdr = "Extract example header from {:s}\n\tand save as {:s}"
        str_keywords = "All {:d} keywords will be loaded."
        str_keyerror_fill = "Key {:s} not found for {:s}, filling with '--'."
        str_valerror = "Please use 'U80' as the dtype for the key {:s}."
        str_filesave = 'Saving the summary file to "{:s}"'

    # Save example header
    if example_header is not None:
        example_fits = filelist[0]
        if verbose:
            print(str_example_hdr.format(str(example_fits), example_header))
        ex_hdu = fits.open(example_fits)
        ex_hdr = ex_hdu[extension].header
        ex_hdr.totextfile(example_header, overwrite=True)

    # load ALL keywords for special cases
    if (keywords == []) or (keywords == '*'):
        example_fits = filelist[0]
        ex_hdu = fits.open(example_fits)
        ex_hdu.verify('fix')
        ex_hdr = ex_hdu[extension].header
        N_hdr = len(ex_hdr.cards)
        keywords = []

        for i in range(N_hdr):
            key_i = ex_hdr.cards[i][0]
            if (key_i in skip_keys):
                continue
            elif (key_i in keywords):
                str_duplicate = "Key {:s} is duplicated! Only first one will be saved."
                print(str_duplicate.format(key_i))
                continue
            keywords.append(key_i)

        if verbose:
            print(str_keywords.format(len(keywords)))
#            except fits.VerifyError:
#                str_unparsable = '{:d}-th key is skipped since it is unparsable.'
#                print(str_unparsable.format(i))
#                continue

    # Initialize
    if len(dtypes) == 0:
        dtypes = ['U80'] * len(keywords)
        # FITS header MUST be within 80 characters! (FITS standard)

    summarytab = Table(names=keywords, dtype=dtypes)
    fnames = []

    # Run through all the fits files
    for fitsfile in filelist:
        fnames.append(_get_fname(fitsfile))
        hdu = fits.open(fitsfile)
        hdu.verify('fix')
        hdr = hdu[extension].header
        row = []
        for key in keywords:
            try:
                row.append(hdr[key])
            except KeyError:
                if verbose:
                    print(str_keyerror_fill.format(key, str(fitsfile)))
                try:
                    row.append('--')
                except ValueError:
                    raise ValueError(str_valerror.format('U80'))
        summarytab.add_row(row)
        hdu.close()

    # Attache the file name, and then sort.
    fnames = Column(data=fnames, name='file')
    summarytab.add_column(fnames, index=0)
    summarytab.sort(sort_by)

    tmppath = Path('tmp.csv')
    summarytab.write(tmppath, format=format)
    summarytab = Table.read(tmppath, format=format)

    if output is None or output == '':
        tmppath.unlink()

    else:
        output = Path(output)
        if verbose:
            print(str_filesave.format(str(output)))
        tmppath.rename(output)

    return summarytab




# open windows
class MyApp(QMainWindow):
    
    def __init__(self):
        super().__init__() #super는 부모클래스(여기선 QWidget)의 메소드를 불러와주기 위해서 사용
        
        #fileTable에 사용할 panda 프레임 만들어놓기
        data = {'FILE-NAME':[' '], 'DATE-OBS':[' '], 'EXPTIME':[' '], 'IMAGETYPE':[' '], 'OBJECT':[' '], 'REMARKS':[' ']}
        self.fitFileList = pd.DataFrame(data)
        data = {'FILE-NAME':[' '], 'REMARKS':[' ']}
        self.graphFileList = pd.DataFrame(data)       
        
        
        
        
        
        
        
        
        
        
        
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
        
        self.fitFileTable = pandaTableWidget(self.fitFileList)
        
        self.graphFileTable = pandaTableWidget(self.graphFileList)
        
        self.mainSplitter.addWidget(QFrame())
        self.mainSplitter.addWidget(self.fileSplitter)
        
        
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
        
        #PreProcessing
        self.preProcessingAction = QAction('Preprocessing', self)
        self.preProcessingAction.setShortcut('Ctrl+P')
        self.preProcessingAction.triggered.connect(self.close)
        self.preProcessingAction.setStatusTip('Crop and Preprocess images')
        
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

        #ToolBar
        #self.toolbar = self.addToolBar('Exit')
        #self.toolbar.addAction(exitAction)
        
        
        #date and time
        #self.statusBar().showMessage(self.date.toString(Qt.DefaultLocaleLongDate))
        
        
                
    def fileOpen(self):
        file = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        toppath = Path(file)
        files = list(toppath.glob("*.fit"))
        print(files)
        summarytab = make_summary(files, keywords=["DATE-OBS", "EXPTIME", "IMAGETYP", "OBJECT"], output=toppath / "summary_20181023.csv")
        print('='*80)
        
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
    
    def Start(self):
        self.statusBar().showMessage('processing')   
    
        



if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = MyApp()
    ex.show()
    sys.exit(app.exec_())