import pandas as pd

class flags():
    def __init__(self):
        self.isCropped = False
        self.isReduced = False

class cropInfo:
    def __init__(self, y0=0, y1=0, x0=0, x1=0, filename=''):
        self.y0 = y0
        self.y1 = y1
        self.x0 = x0
        self.x1 = x1
        self.filename = filename

    def __repr__(self):
        return str([self.y0, self.y1, self.x0, self.x1])




class currentFileInfo():
    def __init__(self, fitFileList : pd.DataFrame, currentFolderLocation='', currentFileLocation ='' , cropInfo = cropInfo() , flags = flags(), fitInfo='default' ):
        self.fitFileList = fitFileList
        self.flags = flags
        self.currentFolderLocation = currentFolderLocation
        self.currentFileLocation = currentFileLocation
        self.cropInfo = cropInfo
        self.fitInfo = fitInfo

