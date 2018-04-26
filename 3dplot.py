# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 09:05:03 2018

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
          'figure.figsize': (40, 30),
         'axes.labelsize': '40',
         'axes.titlesize':'46',
         'xtick.labelsize':'36',
         'ytick.labelsize':'36'}
plt.rcParams.update(params)

angle = -45


X = np.genfromtxt('X.txt',autostrip=True)
Y = np.genfromtxt('Y.txt',autostrip=True)
#X, Y = np.meshgrid(X, Y)
Z = np.genfromtxt('Z.txt',autostrip=True)

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#surf = ax.plot_surface(Y, Z, X, colors, antialiased=True, shade=0.1)


# fourth dimention - colormap
# create colormap according to x-value (can use any 50x50 array)
color_dimension = Z # change to desired fourth dimension
minn, maxx = color_dimension.min(), color_dimension.max()
norm = matplotlib.colors.Normalize(minn, maxx)
m = plt.cm.ScalarMappable(norm=norm, cmap='jet')
m.set_array([])
fcolors = m.to_rgba(color_dimension)

# plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(Y,Z,X, facecolors=fcolors, vmin=minn, vmax=maxx, shade=False)

#fig.colorbar(surf, shrink=0.5, aspect=20)

ax.set_zlim(-0.1,100)
ax.set_xlabel('Y Coordinate (m)')
ax.set_ylabel('u/U (%)')
ax.set_zlabel('Height (mm)')
ax.yaxis.labelpad=40
ax.xaxis.labelpad=40
ax.zaxis.labelpad=40
v=np.linspace(0,100,10, endpoint=True)




#ax.view_init(30, angle)
#plt.draw()
#plt.pause(.001)


