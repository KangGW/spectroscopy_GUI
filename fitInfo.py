class flags():
    def __init__(self):
        self.isCropped = False
        self.isReduced = False

class cropInfo:
    def __init__(self):
        self.y0 = 0
        self.y1 = 0
        self.x0 = 0
        self.x1 = 0
        self.filename = ''

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

