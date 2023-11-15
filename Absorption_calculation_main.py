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
eps_1           = cla2.epsilon(3.3,	4.8)             #NiFe(Py)
eps_2           = cla2.epsilon(0.57,8.0776)                  # Pt       
eps_3           = cla2.epsilon(1.78,	0.0000049)                  #Sub(Sapphire)
eps_4           = cla2.epsilon(1,	0)                      #quartz
eps_5           = cla2.epsilon(1,	0)               #whatever
epsilon_optical = [eps_Air*8.85e-12,eps_1*8.85e-12,eps_2*8.85e-12,eps_3*8.85e-12,eps_4*8.85e-12,eps_5*8.85e-12]  #permittivity 
mu              = 12.57e-7 



d1 = 3e-9
d2 = 6e-9
d3 = 0.5e-3
d4 = 0
d5 = 0
number_of_thickness = 1
##########################################################################




ZERO            = np.linspace(0,0,number_of_thickness)
LAYER1          = np.linspace(d1,d1,number_of_thickness)
LAYER2          = np.linspace(d2,d1,number_of_thickness)       
LAYER3          = np.linspace(d3,d3,number_of_thickness)         
LAYER4          = np.linspace(d4,d4,number_of_thickness)
LAYER5          = np.linspace(d5,d5,number_of_thickness)
dnew            = np.transpose([ZERO,LAYER1,LAYER2,LAYER3,LAYER4,LAYER5])

change_Pt_thickness = cla2.Change_thickness(f_omega_single,omega_single,mu,epsilon_optical)
Q_heatloss1, Q_heatloss2, Q_heatloss3, Q_heatloss4, Q_heatloss5, Q_income         = change_Pt_thickness.heatloss(dnew,number_of_thickness)
Absorption1, Absorption2, Absorption3, Absorption4, Absorption5, Absorption_total = change_Pt_thickness.absorption(dnew,number_of_thickness)
## sub(quartz)         Pt            Si        Co          Sio2     Total


reflection      = change_Pt_thickness.reflection(dnew,number_of_thickness)
transmission    = change_Pt_thickness.transmission(dnew,number_of_thickness)

Total           = reflection+transmission+Absorption_total
Absorp          = 1-reflection-transmission
r,t             = change_Pt_thickness.trans_reflect(dnew,number_of_thickness)

atran_heat,atrans_in = change_Pt_thickness.transmission_heat(dnew,number_of_thickness)

matrix       = cla2.SpecialMatrix(omega_single,mu,epsilon_optical,dnew[0])
T_0s, T_sinf, T_0inf                     = matrix.Transfer_Matrix()
at        = matrix.Transmission_coeff(T_0inf)#(-0.2126720162058721+0.7503224950023671j)
t_os        = matrix.Transmission_coeff(T_0s)#(-0.258334517143757+0.928450239723685j)
t_sinf        = matrix.Transmission_coeff(T_sinf)
at_noecho =np.multiply(t_os, t_sinf) 

atran = f_omega_single*at/f_omega_single
atran_noecho =f_omega_single*at_noecho/f_omega_single 





print('Py/Pt/Sub:')
print('Absorption_Py = ', Absorption1)
print('Absorption_PtTe2     = ', Absorption2)
print('Absorption_Sub   = ', Absorption3)
print('Absorption_total 0.4462    = ', Absorption_total)
print('Transmission     0.2908   = ',transmission)
print('Refelction       0.2630     = ',reflection)
print('A+T+R             = ',Absorption_total+transmission+reflection)


#####PtTe2/Sub ONly

# print('Py/Pt/Sub:')
# print('Absorption_Py = ', Absorption1)
# print('Absorption_PtTe2     = ', Absorption2)
# print('Absorption_Sub   = ', Absorption3)
# # print('Absorption_Pt     = ', Absorption4)
# # print('Absorption_Sio2 = ', Absorption5)
# print('Absorption_total 0.3079    = ', Absorption_total)
# print('Transmission     0.5127   = ',transmission)
# print('Refelction       0.1794     = ',reflection)
# print('A+T+R             = ',Absorption_total+transmission+reflection)



#####Sub ONly


# print('Py/Pt/Sub:')
# print('Absorption_Py = ', Absorption1)
# print('Absorption_Pt     = ', Absorption2)
# print('Absorption_Sub   = ', Absorption3)
# # print('Absorption_Pt     = ', Absorption4)
# # print('Absorption_Sio2 = ', Absorption5)
# print('Absorption_total  0.038  = ', Absorption_total, 'da = ', Absorption_total-0.0387)
# print('Transmission   0.863456   = ',transmission, 'dt = ', transmission-0.8635)
# print('Refelction    0.0978392     = ',reflection, 'dr = ', reflection-0.0978)
# print('A+T+R             = ',Absorption_total+transmission+reflection)


