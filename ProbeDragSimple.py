# -*- coding: utf-8 -*-
"""
ATE Aerotech 4-Axis Probe Manipulator Drag Calculations

Created on Thu Feb 15 08:33:25 2018

@author: BrookeBT
"""

import matplotlib
import numpy as np
import numpy.linalg as la
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy as sci
params = {'legend.fontsize': '40',
          'figure.figsize': (20, 20),
         'axes.labelsize': '40',
         'axes.titlesize':'46',
         'xtick.labelsize':'36',
         'ytick.labelsize':'36'}
plt.rcParams.update(params)


" Variables "
Uinf = 310 #km/hr
T0c = 22 #deg C
T0k = T0c+273 #Kelvin
Pref = 104978 #Pa
Rair = 287.058 
Rv = 461.495
Mair = 0.028964
Mv = 0.018016
R = 8.314
RH = 0.5

" Flow Calculations "
Uinf = Uinf*5/18
Psat = 6.1078*10**((17.27*T0c)/(T0c+237.3))
Pv = RH*Psat
Pair = Pref-Pv
rho = (Pair*Mair+Pv*Mv)/(R*T0k)
u0 = np.array([1, 0, 0])
U0 = Uinf*u0

" Coordinate Transformations "
alpha = np.deg2rad(-90)
beta = np.deg2rad(45)
gamma = np.deg2rad(-30)
Dmax_mag = 1
Mmax_mag = 1
DXmax = 0
DYmax = 0
DZmax = 0
MXmax = 0
MYmax = 0
MZmax = 0

r0 = np.array([0, 0, 0])
r1_0 = np.array([0, 0, -0.160])
U1 = U0
r1 = r0+r1_0
#
R2 = np.array([[1, 0, 0],
              [0, np.cos(alpha), np.sin(alpha)],
              [0, -1*np.sin(alpha), np.cos(alpha)]])
#
r2_1 = np.array([-0.360, 0, -0.650])
U2 = np.dot(R2,U1)
r2 = np.dot(la.inv(R2),r2_1)+r1
#
R3 = np.array([[np.cos(beta), -1*np.sin(beta), 0],
              [np.sin(beta), np.cos(beta), 0],
              [0, 0, 1]])

r3_2 = np.array([0, 0, -0.650])
U3 = np.dot(R3,U2)
r3 = np.dot(np.dot(la.inv(R2),la.inv(R3)),r3_2)+r2
#
R4 = np.array([[np.cos(gamma), 0, -1*np.sin(gamma)],
              [0, 1, 0],
              [np.sin(gamma), 0, np.cos(gamma)]])
#
r4_3 = np.array([-0.025, 0, 0])
U4 = np.dot(R4,U3)
r4 = np.dot(np.dot(np.dot(la.inv(R2),la.inv(R3)),la.inv(R4)),r4_3)+r3
u4 = U4/Uinf
#
" Torpedo Base Drag and Moment "
d1 = .160 #meters
l1 = 1.1 #m
#
Cd1Par = 0.2 #Approximate drag for streamlined cylinder with l/d=7 (Hoerner)
Cd1Per = 1.2 #Approximate drag for cylinder perpendicular to flow
#
D1X = Cd1Par*0.5*rho*U1[0]**2*np.pi*d1**2/4
D1Y = Cd1Per*0.5*rho*U1[1]**2*l1*d1
D1Z = Cd1Per*0.5*rho*U1[2]**2*l1*d1
D1 = np.array([D1X, D1Y, D1Z])
M1 = np.cross(r1,D1)
#
" Cylindrical Arm Drag and Moment "
d2_1 = 0.15 #m
l2_1 = 0.24 #m
d2_2 = (0.09+0.075)/2 #m
l2_2 = 0.925 #m
Cd2Par = 0.2 #Approximate drag for streamlined cylinder (Hoerner)
Cd2Per = 1.2 #Approximate drag for cylinder perpendicular to flow
D2_1X = Cd2Per*0.5*rho*U2[0]**2*l2_1*d2_1
D2_1Y = Cd2Per*0.5*rho*U2[1]**2*l2_1*d2_1
D2_1Z = Cd2Par*0.5*rho*U2[2]**2*np.pi*d2_1**2/4
D2_1 = np.array([D2_1X, D2_1Y, D2_1Z])
D2_2X = Cd2Per*0.5*rho*U2[0]**2*l2_2*d2_2
D2_2Y = Cd2Per*0.5*rho*U2[1]**2*l2_2*d2_2
D2_2Z = Cd2Par*0.5*rho*U2[2]**2*np.pi*d2_2**2/4
D2_2 = np.array([D2_2X, D2_2Y, D2_2Z])
D2 = D2_1+D2_2
M2 = np.cross(r3,np.dot(la.inv(R2),D2))
#
" Probe Holder Drag "
d3 = 0.06 #m
l3 = 0.7 #m
Cd3Par = 0.3 #Approximate drag for streamlined cylinder of L/D = 12
Cd3Per = 1.2 #Approximate drag for cylinder perpendicular to flow
D3X = Cd3Par*0.5*rho*U4[0]**2*np.pi*d3**2/4
D3Y = Cd3Per*0.5*rho*U4[1]**2*l3*d3
D3Z = Cd3Per*0.5*rho*U4[2]**2*l3*d3
D3 = np.array([D3X, D3Y, D3Z])
M3 = np.cross(r4,np.dot(np.dot(np.dot(la.inv(R2),la.inv(R3)),la.inv(R4)),D3))
#
D = D1+np.dot(la.inv(R2),D2)+np.dot(np.dot(np.dot(la.inv(R2),la.inv(R3)),la.inv(R4)),D3)
M = M1+M2+M3
#

    
#print("D = ",D)
#print("M = ",M)
print()
print("D = ", D)
print("M = ", M)
print("alpha = ", np.rad2deg(alpha))
print("beta = ", np.rad2deg(beta))
print("gamma = ", np.rad2deg(gamma))
