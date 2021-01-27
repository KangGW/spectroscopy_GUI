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
from skimage.feature import peak_local_max
from astropy.modeling.fitting import LevMarLSQFitter
from scipy import stats
#우선 있는 파일 중에 Cali 파일만 열어서 first identification 하는 프로그램
#cali fileinfo를 받아와서 연 다음에

#Todo 지금 있는 rawSpectrum 임시로 그린거니까 다시그리기
# -> neonarc.fit은 지금 있는 cali파일에서 아이덴티피케이션만 해서 그린거라 부정확, 신뢰성있는 스펙트럼으로 바꿔서 넣자.
#Todo Neon 말고 다른 칼리브레이션 아크를 쓸 경우 선택할 수 있게 만들자
#Todo 패턴매칭 알고리즘을 사용해 자동으로 비슷한 이미지구간을 찾아주는 기능도 구현해보자








B = fits.open("./Spectroscopy_Example/20181023/combine/comp_15_15.0.fit")

identify = np.average(B[0].data[50:70,:]- C[0].data[50:70,:], axis=0)
max_intens = np.max(identify)
MINSEP_PK = 5   # minimum separation of peaks
MINAMP_PK = 0.01 # fraction of minimum amplitude (wrt maximum) to regard as peak
NMAX_PK = 50
FWHM_ID = 4

peak_pix = peak_local_max(identify, indices=True, num_peaks=NMAX_PK,
                          min_distance=MINSEP_PK,
                          threshold_abs=max_intens * MINAMP_PK)
ID_init = dict (pixel_init=[403, 374, 362, 265,
                           245, 238, 223, 213,
                           178, 169, 144, 131, 53],
wavelength=[8780.6, 8495.4, 8377.6, 7438.9,
                           7245.2, 7173.9, 7032.4, 6929.5,
                           6599.0, 6507.0, 6266.5, 6143.1, 5400.6])

fitter = LevMarLSQFitter()

x_identify = np.arange(len(identify))
peak_gauss = []
for peak_pix in ID_init['pixel_init']:
    #TODO: put something like "lost_factor" to multiply to FWHM_ID in the bounds.
    g_init = Gaussian1D(amplitude = identify[peak_pix],
                       mean = peak_pix,
                       stddev = FWHM_ID * gaussian_fwhm_to_sigma,
                       bounds={'amplitude': (0, 2*identify_1[peak_pix]),
                               'mean':(peak_pix-FWHM_ID, peak_pix+FWHM_ID),
                               'stddev':(0, FWHM_ID)})
    fitted = fitter(g_init, x_identify, identify)
    peak_gauss.append(fitted.mean.value)

res = stats.linregress(peak_gauss, ID_init['wavelength'])
wavelength  = x_identify*res.slope + res.intersections

for i in peak_gauss:
    plt.axvline(i,0,100000, color = 'c')
    plt.text(i,40000,str(i), rotation=70)
plt.plot(identify, color = 'r')
plt.xlim(400,800)
plt.show()

