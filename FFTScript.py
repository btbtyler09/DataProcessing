# -*- coding: utf-8 -*-
"""
Created on Fri Apr 20 11:46:10 2018

@author: BrookeBT
"""

import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy as sci
from mpl_toolkits.mplot3d import Axes3D
params = {'legend.fontsize': '40',
          'figure.figsize': (62/2, 38/2),
         'axes.labelsize': '40',
         'axes.titlesize':'46',
         'xtick.labelsize':'36',
         'ytick.labelsize':'36'}
plt.rcParams.update(params)



def ButterBandpassFilter(Signal,CutoffLow, CutoffHigh, Fs, order=8):
    Nyquist = 0.5*Fs
    LowFraction = CutoffLow/Nyquist
    HighFraction = CutoffHich/Nyquist
    b, a = scipy.signal.butter(order, [LowFraction,HighFraction], btype='band')
    FilteredSignal = filtfilt(b,a,Signal)
    
    
Signal = np.genfromtxt('Signal.txt',autostrip=True)    

window_size = 10;
BlockWidth = 4096;
Overlap = 0.50;

## Filtering data
numberofsamples = length(signal);
b,a = scipy.signal.butter(8,CutoffLow.98);
signal = filtfilt(b,a,signal);