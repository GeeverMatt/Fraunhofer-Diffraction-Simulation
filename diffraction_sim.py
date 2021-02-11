#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:27:55 2021

@author: matthew
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import special

#setting variables (units in SI)


p = np.linspace(-.01,.01,1000) #size of screen, no. pixels
d = .4 #distance to screen
wl = 700e-9 #wavelength of light
irr = 100 #irradiance (I0)
k = (2 * np.pi)/wl 

#function for generating angles for each pixel
def angle(p):
    return np.arctan(p/d)


theta = angle(p) #generate array of angles

#Frauhofer Diffraction by a single slit

b = 100e-6 #slit width

#beta values in function
def beta(x):
    return (k * b/2)*np.sin(x)

#function for fraunhofer diffraction by single slit using beta function
def I_fhr(x):
    return irr * ((np.sin(beta(x))/beta(x))**2)

diff_fhr = I_fhr(theta) #Intensities array

#plotting the pattern

fig_fhr = plt.figure()
ax_fhr = fig_fhr.add_subplot(111)

pattern_fhr = ax_fhr.plot(p,diff_fhr,'r')

ax_fhr.set_title("Fraunhofer Diffraction by a single slit",fontweight='bold')
ax_fhr.set_xlabel("Distance from centre of pattern (m)")
ax_fhr.set_ylabel("Intensity")

#Fraunhofer Diffraction by a pinhole

a = 100e-6 #Pin hole diameter

#Fraunhofer diffraction by pinhole funtion using first order bessel funtion
def I_fhr_p(x):
    return irr * ( ((2 * scipy.special.j1(k * a * np.sin(x)))/(k * a * np.sin(x)))**2 )

diff_fhr_p = I_fhr_p(theta) #Intensities array

fig_fhr_p = plt.figure()
ax_fhr_p = fig_fhr_p.add_subplot(111)

pattern_fhr_p = ax_fhr_p.plot(p,diff_fhr_p,'r')

ax_fhr_p.set_title("Fraunhofer Diffraction by a pinhole",fontweight='bold')
ax_fhr_p.set_xlabel("Distance from centre of pattern (m)")
ax_fhr_p.set_ylabel("Intensity")

#Fraunhofer Diffraction produced by diffraction grating
b2 = 40e-6 #slit width
a2 = 80e-6 #slit separation
N = 25 #number of slits

#new beta and alpha functions for diffraction grating
def beta_g(x):
    return (k * b2 / 2) * np.sin(x)

def alpha_g(x): 
    return (k * a2 / 2) * np.sin(x)

#Fraunhofer diffraction by diffraction grating function
def I_grating(x):
    return irr * ((np.sin(beta_g(x))/beta_g(x))**2) * ((np.sin(N*alpha_g(x))/np.sin(alpha_g(x)))**2)

diff_grat = I_grating(theta) #Intensities array

fig_grat = plt.figure()
ax_grat = fig_grat.add_subplot(111)

pattern_grat = ax_grat.plot(p,diff_grat,'r')

ax_grat.set_title("Fraunhofer Diffraction by a diffraction grating",fontweight='bold')
ax_grat.set_xlabel("Distance from centre of pattern (m)")
ax_grat.set_ylabel("Intensity")

plt.show()