#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:46:22 2023

@author: yingshuyang
"""



import numpy as np
import matplotlib.pyplot as plt
import absorption_calculation_class as cla2
from math import pi








## =============================================================================
## 1.Inputs
## =============================================================================

'''Enter the frequency of the laser of which you are interested (units:Hz)''' 
frequency       = 374.74e12    
omega_single    = 2*pi*frequency

#the freq amp we set manully, this do not matter since we are finally going to calculte the percentages
f_omega_single  = 10


''' Enter the dielectric constants for each layer for the frequency you entered above'''
eps_Air         = cla2.epsilon(1,	0)                      # left most Air layer values n = 1, k = 0
eps_1           = cla2.epsilon(3.8,	3.89)                   # Layer1 values n = x, k = y     
eps_2           = cla2.epsilon(1.4608,	0.0013066)          # Layer2 values n = x, k = y    
eps_3           = cla2.epsilon(1,	0)                      # Layer3 values n = x, k = y              
eps_4           = cla2.epsilon(1,	0)                      # Layer4 values n = x, k = y                   
eps_5           = cla2.epsilon(1,	0)                      # Layer5 values n = x, k = y 
epsilon_optical = [eps_Air*8.85e-12,
                   eps_1*8.85e-12,
                   eps_2*8.85e-12,
                   eps_3*8.85e-12,
                   eps_4*8.85e-12,
                   eps_5*8.85e-12]                          #permittivity 
mu              = 12.57e-7                                  #permeability (set as constant)


'''Enter the changing thickness of each layer (units:m):
   Note: if the layer thickness do not change, just set the dmin and dmax as the same thickness
'''
d1min = 0       ; d1max = 100e-9
d2min = 0.3e-3  ; d2max = 0.3e-3
d3min = 1e-9    ; d3max = 1e-9
d4min = 0       ; d4max = 0
d5min = 0       ; d5max = 0

'''Specify the layer that the thickness is changing:
    1- the first layer
    2- the second layer and so on...(in the example case is 1)'''
    
thickness_changing_layer = 1

'''The resolution of the thikness you want'''
number_of_thickness = 20

## =============================================================================
## 2.Runniing of the coe
## =============================================================================




ZERO            = np.linspace(0,0,number_of_thickness)
LAYER1          = np.linspace(d1min,d1max,number_of_thickness)
LAYER2          = np.linspace(d2min,d2max,number_of_thickness)       
LAYER3          = np.linspace(d3min,d3max,number_of_thickness)         
LAYER4          = np.linspace(d4min,d4max,number_of_thickness)
LAYER5          = np.linspace(d5min,d5max,number_of_thickness)

d = [ZERO,LAYER1,LAYER2,LAYER3,LAYER4,LAYER5]
dnew            = np.transpose([ZERO,LAYER1,LAYER2,LAYER3,LAYER4,LAYER5])

change_Pt_thickness = cla2.Change_thickness(f_omega_single,omega_single,mu,epsilon_optical)
Q_heatloss1, Q_heatloss2, Q_heatloss3, Q_heatloss4, Q_heatloss5, Q_income         = change_Pt_thickness.heatloss(dnew,number_of_thickness)
Absorption1, Absorption2, Absorption3, Absorption4, Absorption5, Absorption_total = change_Pt_thickness.absorption(dnew,number_of_thickness)


reflection      = change_Pt_thickness.reflection(dnew,number_of_thickness)
transmission    = change_Pt_thickness.transmission(dnew,number_of_thickness)

Total           = reflection+transmission+Absorption_total
Absorp          = 1-reflection-transmission
r,t             = change_Pt_thickness.trans_reflect(dnew,number_of_thickness)

## =============================================================================
## 3. Result printing
## =============================================================================

plt.figure('Absorption change')
plt.title('ART and absorption change in each layer')
plt.plot(d[1],transmission,linewidth = 3,label = 'transmission')
plt.plot(d[1],reflection,linewidth = 3,label = 'reflection')
plt.plot(d[1],Absorption_total,linewidth = 3,label = 'absorption total')
plt.plot(d[1],Absorption1,'--',label = 'absorption layer1')
plt.plot(d[1],Absorption2,'--',label = 'absorption layer2')
plt.xlabel('Thickness (m)')
plt.ylabel('Ratio')
plt.ylim(-0.1,1.1)
plt.legend()

print('Results:')
print('Absorption_layer1 = ', Absorption1)
print('Absorption_layer2 = ', Absorption2)
print('Absorption_layer3 = ', Absorption3)
print('Absorption_layer4 = ', Absorption4)
print('Absorption_layer5 = ', Absorption5)
print('Absorption__total = ', Absorption_total)
print('Transmission      = ',transmission)
print('Refelction        = ',reflection)
print('A+T+R             = ',Absorption_total+transmission+reflection)




