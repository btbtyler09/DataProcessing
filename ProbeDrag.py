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
Uinf = 310*0.7 #km/hr
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
FlowAngle = np.deg2rad(10)
Psat = 6.1078*10**((17.27*T0c)/(T0c+237.3))
Pv = RH*Psat
Pair = Pref-Pv
rho = (Pair*Mair+Pv*Mv)/(R*T0k)
u0 = np.array([np.cos(FlowAngle), np.sin(FlowAngle), 0])
U0 = Uinf*u0

" Coordinate Transformations "
#alpha = np.deg2rad(0)
#beta = np.deg2rad(0)
#gamma = np.deg2rad(0)
Dmax_mag = 1
Mmax_mag = 1
DXmax = 0
DYmax = 0
DZmax = 0
MXmax = 0
MYmax = -1000
MZmax = 0
DXmin = 1000
DYmin = 0
DZmin = 0
MXmin = 0
MYmin = 0
MZmin = 0

for ii in range(-160,161):
    print(ii)
    for jj in range(-90,91):
        for kk in range(-45,46):
            alpha = np.deg2rad(ii)
            beta = np.deg2rad(jj)
            gamma = np.deg2rad(kk)            
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
            #
            if la.norm(D) > Dmax_mag:
                Dmax_mag = la.norm(D)
                Dmax = D
                alpha_Dmax = alpha
                beta_Dmax = beta
                gamma_Dmax = gamma
                M_Dmax = M
            
            if la.norm(M) > Mmax_mag:
                Mmax_mag = la.norm(M)
                Mmax = M
                alpha_Mmax = alpha
                beta_Mmax = beta
                gamma_Mmax = gamma
                D_Mmax = D
                
            if D[0] > DXmax:
                DXmax = D[0]
                alpha_DXmax = alpha
                beta_DXmax = beta
                gamma_DXmax = gamma
                M_DXmax = M
                D_DXmax = D
            
            if D[1] > DYmax:
                DYmax = D[1]
                alpha_DYmax = alpha
                beta_DYmax = beta
                gamma_DYmax = gamma
                M_DYmax = M
                D_DYmax = D
                
            if D[2] > DZmax:
                DZmax = D[2]
                alpha_DZmax = alpha
                beta_DZmax = beta
                gamma_DZmax = gamma
                M_DZmax = M
                D_DZmax = D
                
            if D[0] < DXmin:
                DXmin = D[0]
                alpha_DXmin = alpha
                beta_DXmin = beta
                gamma_DXmin = gamma
                M_DXmin = M
                D_DXmin = D
            
            if D[1] < DYmin:
                DYmin = D[1]
                alpha_DYmin = alpha
                beta_DYmin = beta
                gamma_DYmin = gamma
                M_DYmin = M
                D_DYmin = D
                
            if D[2] < DZmin:
                DZmin = D[2]
                alpha_DZmin = alpha
                beta_DZmin = beta
                gamma_DZmin = gamma
                M_DZmin = M
                D_DZmin = D
                
            if M[0] > MXmax:
                MXmax = M[0]
                alpha_MXmax = alpha
                beta_MXmax = beta
                gamma_MXmax = gamma
                D_MXmax = D
                M_MXmax = M
                
            if M[1] > MYmax:
                MYmax = M[1]
                alpha_MYmax = alpha
                beta_MYmax = beta
                gamma_MYmax = gamma
                D_MYmax = D
                M_MYmax = M
                
            if M[2] > MZmax:
                MZmax = M[2]
                alpha_MZmax = alpha
                beta_MZmax = beta
                gamma_MZmax = gamma
                D_MZmax = D
                M_MZmax = M
                
            if M[0] < MXmin:
                MXmin = M[0]
                alpha_MXmin = alpha
                beta_MXmin = beta
                gamma_MXmin = gamma
                D_MXmin = D
                M_MXmin = M
                
            if M[1] < MYmin:
                MYmin = M[1]
                alpha_MYmin = alpha
                beta_MYmin = beta
                gamma_MYmin = gamma
                D_MYmin = D
                M_MYmin = M
                
            if M[2] < MZmin:
                MZmin = M[2]
                alpha_MZmin = alpha
                beta_MZmin = beta
                gamma_MZmin = gamma
                D_MZmin = D
                M_MZmin = M
                
            #print("D = ",D)
            #print("M = ",M)
print("--Maximum Aerodynamic Loads on 4-Axis Probe Holder--")
print("----Forces in Newtons & Moments in Newton Meters----")            

print()
print("--Maximum Force--")
print("Dmax = ", Dmax)
print("M_Dmax = ", M_Dmax)
print("alpha_Dmax = ", np.rad2deg(alpha_Dmax))
print("beta_Dmax = ", np.rad2deg(beta_Dmax))
print("gamma_Dmax = ", np.rad2deg(gamma_Dmax))

print()
print("--Maximum Moment--")
print("Mmax = ", Mmax)
print("D_Mmax = ", D_Mmax)
print("alpha_Mmax = ", np.rad2deg(alpha_Mmax))
print("beta_Mmax = ", np.rad2deg(beta_Mmax))
print("gamma_Mmax = ", np.rad2deg(gamma_Mmax))

print()
print("--Maximum X Force--")
print("DXmax = ", DXmax)
print("D_DXmax = ", D_DXmax)
print("M_DXmax = ", M_DXmax)
print("alpha_DXmax = ", np.rad2deg(alpha_DXmax))
print("beta_DXmax = ", np.rad2deg(beta_DXmax))
print("gamma_DXmax = ", np.rad2deg(gamma_DXmax))

print()
print("--Maximum Y Force--")
print("DYmax = ", DYmax)
print("D_DYmax = ", D_DYmax)
print("M_DYmax = ", M_DYmax)
print("alpha_DYmax = ", np.rad2deg(alpha_DYmax))
print("beta_DYmax = ", np.rad2deg(beta_DYmax))
print("gamma_DYmax = ", np.rad2deg(gamma_DYmax))

print()
print("--Maximum Z Force--")
print("DZmax = ", DZmax)
print("D_DZmax = ", D_DZmax)
print("M_DZmax = ", M_DZmax)
print("alpha_DZmax = ", np.rad2deg(alpha_DZmax))
print("beta_DZmax = ", np.rad2deg(beta_DZmax))
print("gamma_DZmax = ", np.rad2deg(gamma_DZmax))

print()
print("--Minimum X Force--")
print("DXmin = ", DXmin)
print("D_DXmin = ", D_DXmin)
print("M_DXmin = ", M_DXmin)
print("alpha_DXmin = ", np.rad2deg(alpha_DXmin))
print("beta_DXmin = ", np.rad2deg(beta_DXmin))
print("gamma_DXmin = ", np.rad2deg(gamma_DXmin))

print()
print("--Minimum Y Force--")
print("DYmin = ", DYmin)
print("D_DYmin = ", D_DYmin)
print("M_DYmin = ", M_DYmin)
print("alpha_DYmin = ", np.rad2deg(alpha_DYmin))
print("beta_DYmin = ", np.rad2deg(beta_DYmin))
print("gamma_DYmin = ", np.rad2deg(gamma_DYmin))

print()
print("--Minimum Z Force--")
print("DZmin = ", DZmin)
print("D_DZmin = ", D_DZmin)
print("M_DZmin = ", M_DZmin)
print("alpha_DZmin = ", np.rad2deg(alpha_DZmin))
print("beta_DZmin = ", np.rad2deg(beta_DZmin))
print("gamma_DZmin = ", np.rad2deg(gamma_DZmin))

print()
print("--Maximum X Moment--")
print("MXmax = ", MXmax)
print("D_MXmax = ", D_MXmax)
print("M_MXmax = ", M_MXmax)
print("alpha_MXmax = ", np.rad2deg(alpha_MXmax))
print("beta_MXmax = ", np.rad2deg(beta_MXmax))
print("gamma_MXmax = ", np.rad2deg(gamma_MXmax))

print()
print("--Maximum Y Moment--")
print("MYmax = ", MYmax)
print("D_MYmax = ", D_MYmax)
print("M_MYmax = ", M_MYmax)
print("alpha_MYmax = ", np.rad2deg(alpha_MYmax))
print("beta_MYmax = ", np.rad2deg(beta_MYmax))
print("gamma_MYmax = ", np.rad2deg(gamma_MYmax))

print()
print("--Maximum Z Moment--")
print("MZmax = ", MZmax)
print("D_MZmax = ", D_MZmax)
print("M_MZmax = ", M_MZmax)
print("alpha_MZmax = ", np.rad2deg(alpha_MZmax))
print("beta_MZmax = ", np.rad2deg(beta_MZmax))
print("gamma_MZmax = ", np.rad2deg(gamma_MZmax))

print()
print("--Minimum X Moment--")
print("MXmin = ", MXmin)
print("D_MXmin = ", D_MXmin)
print("M_MXmin = ", M_MXmin)
print("alpha_MXmin = ", np.rad2deg(alpha_MXmin))
print("beta_MXmin = ", np.rad2deg(beta_MXmin))
print("gamma_MXmin = ", np.rad2deg(gamma_MXmin))

print()
print("--Minimum Y Moment--")
print("MYmin = ", MYmin)
print("D_MYmin = ", D_MYmin)
print("M_MYmin = ", M_MYmin)
print("alpha_MYmin = ", np.rad2deg(alpha_MYmin))
print("beta_MYmin = ", np.rad2deg(beta_MYmin))
print("gamma_MYmin = ", np.rad2deg(gamma_MYmin))

print()
print("--Minimum Z Moment--")
print("MZmin = ", MZmin)
print("D_MZmin = ", D_MZmin)
print("M_MZmin = ", M_MZmin)
print("alpha_MZmin = ", np.rad2deg(alpha_MZmin))
print("beta_MZmin = ", np.rad2deg(beta_MZmin))
print("gamma_MZmin = ", np.rad2deg(gamma_MZmin))