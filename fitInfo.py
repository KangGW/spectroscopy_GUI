import pandas as pd
import numpy as np

class flags():
    def __init__(self):
        self.isEditted = False
        self.isReduced = False

class editInfo:
    def __init__(self, rotAngle=0, reflection=False, y0=0, y1=0, x0=0, x1=0, filename=''):
        self.y0 = y0
        self.y1 = y1
        self.x0 = x0
        self.x1 = x1
        self.filename = filename
        self.rotAngle = rotAngle
        self.reflection = reflection

    def __repr__(self):
        return str([self.y0, self.y1, self.x0, self.x1])




class currentFileInfo():
    def __init__(self, fitFileList : pd.DataFrame, currentFolderLocation='', currentFileLocation ='', editInfo = editInfo(), flags = flags(), fitInfo='default'):
        self.fitFileList = fitFileList
        self.flags = flags
        self.currentFolderLocation = currentFolderLocation
        self.currentFileLocation = currentFileLocation
        self.editInfo = editInfo
        self.fitInfo = fitInfo


def fitter(fitMethod = 'poly'):
    if fitMethod == 'poly':

        fit = np.polynomial.polynomial.polyfit
        val = np.polynomial.polynomial.polyval

    elif fitMethod == 'legendre':

        fit = np.polynomial.legendre.legfit
        val = np.polynomial.legendre.legval

    elif fitMethod == 'chebyshev':

        fit = np.polynomial.chebyshev.chebfit
        val = np.polynomial.chebyshev.chebval

    return fit, val