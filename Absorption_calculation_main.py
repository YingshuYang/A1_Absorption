#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:46:22 2023

@author: yingshuyang
"""



import numpy as np
import matplotlib.pyplot as plt
# import THz_emission_sim_class as cla
# import class_absorption as cla2
import absorption_calculation_class as cla2
from scipy.interpolate import interp1d
from math import pi
from scipy.optimize import minimize





## =============================================================================
## 1.Absorption calculation for a single frequency
## =============================================================================




##########################################################################
frequency       = 374.74e12    #THz
omega_single    = 2*pi*frequency
f_omega_single  = 10   #the freq amp we set manully



eps_Air         = 1   
eps_1           = cla2.epsilon(3.8,	3.89)             
eps_2           = cla2.epsilon(1.4608,	0.0013066)    
eps_3           = cla2.epsilon(1,	0)                  
eps_4           = cla2.epsilon(1,	0)                     
eps_5           = cla2.epsilon(1,	0)               #
epsilon_optical = [eps_Air*8.85e-12,eps_1*8.85e-12,eps_2*8.85e-12,eps_3*8.85e-12,eps_4*8.85e-12,eps_5*8.85e-12]  #permittivity 
mu              = 12.57e-7 



d1min = 0.00001e-9; d1max = 100e-9#d1 = 100e-9
d2 = 0.3e-3
d3 = 0
d4 = 0
d5 = 0
number_of_thickness = 20



##########################################################################




ZERO            = np.linspace(0,0,number_of_thickness)
LAYER1          = np.linspace(d1min,d1max,number_of_thickness)
LAYER2          = np.linspace(d2,d2,number_of_thickness)       
LAYER3          = np.linspace(d3,d3,number_of_thickness)         
LAYER4          = np.linspace(d4,d4,number_of_thickness)
LAYER5          = np.linspace(d5,d5,number_of_thickness)
dnew            = np.transpose([ZERO,LAYER1,LAYER2,LAYER3,LAYER4,LAYER5])

change_Pt_thickness = cla2.Change_thickness(f_omega_single,omega_single,mu,epsilon_optical)
Q_heatloss1, Q_heatloss2, Q_heatloss3, Q_heatloss4, Q_heatloss5, Q_income         = change_Pt_thickness.heatloss(dnew,number_of_thickness)
Absorption1, Absorption2, Absorption3, Absorption4, Absorption5, Absorption_total = change_Pt_thickness.absorption(dnew,number_of_thickness)


reflection      = change_Pt_thickness.reflection(dnew,number_of_thickness)
transmission    = change_Pt_thickness.transmission(dnew,number_of_thickness)

Total           = reflection+transmission+Absorption_total
Absorp          = 1-reflection-transmission
r,t             = change_Pt_thickness.trans_reflect(dnew,number_of_thickness)


plt.figure('test')
plt.plot(LAYER1,transmission,label = 'transmission')
plt.plot(LAYER1,reflection,label = 'reflection')
plt.plot(LAYER1,Absorption_total,label = 'absorption total')
plt.plot(LAYER1,Absorption1,'--',label = 'absorption CoFeB')
plt.ylim(0,1)
plt.legend()

print('Results')
print('Absorption_layer1 = ', Absorption1)
print('Absorption_layer2 = ', Absorption2)
print('Absorption_layer3 = ', Absorption3)
print('Absorption_layer4 = ', Absorption4)
print('Absorption_layer5 = ', Absorption5)
print('Absorption__total = ', Absorption_total)
print('Transmission = ',transmission)
print('Refelction = ',reflection)
print('A+T+R = ',Absorption_total+transmission+reflection)




