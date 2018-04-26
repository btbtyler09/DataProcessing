# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 05:12:50 2018

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

data = np.genfromtxt('data.txt',autostrip=True)

Length = len(data)

X = np.random.random(Length)
Y = np.random.random(Length)
u = np.random.random(Length)
v = np.random.random(Length)


for ii in range(0,Length):
    X[ii] = data[ii,0]
    Y[ii] = data[ii,1]
    u[ii] = data[ii,2]
    v[ii] = data[ii,3]

alpha_avg = np.average(u)
alpha_std = 2*np.std(u)
beta_avg = np.average(v)
beta_std = 2*np.std(v)
    
plt.figure()
plt.quiver(X,Y,v,u,pivot='mid', width=0.002, headwidth=3, headlength=5, scale=10)
plt.title('Strömungswinkligkeit 140 kph X=-4.5 m Mod 6')
plt.xlabel('Y (mm)')
plt.ylabel('Z (mm)')
plt.annotate('Scale: 0.5°', xy=(-180, 2900), size='28')
plt.annotate('Average α='+str(round(alpha_avg,2))+'°, 2σ='+str(round(alpha_std,2))+'°', xy=(-1600,2900), size='28')
plt.annotate('Average β='+str(round(beta_avg,2))+'°, 2σ='+str(round(beta_std,2))+'°', xy=(600,2900), size='28')
