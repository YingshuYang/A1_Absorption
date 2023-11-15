#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 15:05:54 2021

@author: yingshuyang
"""




import numpy as np

from scipy import interpolate
from math import pi
from scipy.interpolate import interp1d
from numpy.linalg import multi_dot







#=====================================================

def epsilon(n,k):
    eps = (n+1j*k)**2
    return eps
 
        

class SpecialMatrix(object):
    def __init__(self,omega,mu,epsilon,z):
        self.mu = mu
        self.omega = omega
        self.epsilon = epsilon
        self.z = z



    def a_Matrix(self,material,position):    #special
        a11 = np.exp(1j*self.omega*np.sqrt(self.epsilon[material]*self.mu)*self.z[position])
        a12 = np.exp(-1j*self.omega*np.sqrt(self.epsilon[material]*self.mu)*self.z[position])
        a21 =np.sqrt(self.epsilon[material]/self.mu)*np.exp(1j*self.omega*np.sqrt(self.epsilon[material]*self.mu)*self.z[position])
        a22 =-np.sqrt(self.epsilon[material]/self.mu)*np.exp(-1j*self.omega*np.sqrt(self.epsilon[material]*self.mu)*self.z[position])
        am = np.matrix([[a11, a12],[a21, a22]])
        # am = np.array([[a11, a12],[a21, a22]]).transpose(2,0,1)
        return am 
    
    
    def inv_a_Matrix(self,material,position):
        a = self.a_Matrix(material,position)
        a0 = np.linalg.inv(a)
        return a0  
    


    def M_Matrix(self,n):     #special
        m11 = np.cos(self.omega*np.sqrt(self.epsilon[n]*self.mu)*self.z[n])
        m12 = (1j/np.sqrt(self.epsilon[n]/self.mu))*np.sin(self.omega*np.sqrt(self.epsilon[n]*self.mu)*self.z[n])
        m21 = 1j*np.sqrt(self.epsilon[n]/self.mu)*np.sin(self.omega*np.sqrt(self.epsilon[n]*self.mu)*self.z[n])
        m22 = np.cos(self.omega*np.sqrt(self.epsilon[n]*self.mu)*self.z[n])
        m = np.matrix([[m11, m12],[m21, m22]])
        # m = np.array([[m11, m12],[m21, m22]]).transpose(2,0,1)
        return m



    def Total_M(self):
        m = []
        for n in range(1, len(self.z)):
            m.append(self.M_Matrix(n))
        if (len(m) == 4):
            total_M = multi_dot([m[3],m[2],m[1],m[0]])
        elif (len(m) == 3):
            total_M = multi_dot([m[2],m[1],m[0]])
        elif (len(m) == 2):
            total_M = multi_dot([m[1],m[0]])
        else:
            total_M = m[0]            
        return total_M


    def Noecho_M(self):
        m = []
        for n in range(2, len(self.z)):
            m.append(self.M_Matrix(n))
        if (len(m) == 4):
            total_M = multi_dot([m[3],m[2],m[1],m[0]])
        elif (len(m) == 3):
            total_M = multi_dot([m[2],m[1],m[0]])
        elif (len(m) == 2):
            total_M = multi_dot([m[1],m[0]])
        else:
            total_M = m[0]            
        return total_M

    
    
    def Transfer_Matrix(self):
        a0_0 = self.a_Matrix(0,0)
        as_0 = self.inv_a_Matrix(1,0)
        as_d = self.a_Matrix(1,1)
        ainf_0 = self.inv_a_Matrix(0,0)
        
        a0_0 = self.a_Matrix(0,0)
        ainf_0 = self.inv_a_Matrix(0,0)
        a5_d = self.a_Matrix(5,5)
        a5_0 = self.inv_a_Matrix(5,0)
        a4_d = self.a_Matrix(4,4)
        a4_0 = self.inv_a_Matrix(4,0) 
        a3_d = self.a_Matrix(3,3)
        a3_0 = self.inv_a_Matrix(3,0)
        a2_d = self.a_Matrix(2,2)
        a2_0 = self.inv_a_Matrix(2,0)
        as_d = self.a_Matrix(1,1)
        as_0 = self.inv_a_Matrix(1,0)
        
        M_total = self.Total_M()
        M_noecho = self.Noecho_M()        
        T_0s = multi_dot([as_0,a0_0])
        T_sinf = multi_dot([ainf_0,M_noecho,as_d])
        # T_0inf = multi_dot([ainf_0,M_total,a0_0])
        # T_0inf = multi_dot([ainf_0,a4_d,a4_0,a3_d,a3_0,a2_d,a2_0,as_d,as_0,a0_0])
        T_0inf = multi_dot([ainf_0,a5_d,a5_0,a4_d,a4_0,a3_d,a3_0,a2_d,a2_0,as_d,as_0,a0_0])
        return T_0s, T_sinf, T_0inf       


    def Transmission_coeff(self,t):
        t_coeff = (t[0,0]*t[1,1]-t[0,1]*t[1,0])/t[1,1] 
        return t_coeff
    
    def Reflection_coeff(self,t):
        r_coeff= -t[1,0]/t[1,1]
        return r_coeff
        

    def Transfer_Matrix_special(self):
        a0_0 = self.a_Matrix(0,0)
        ainf_0 = self.inv_a_Matrix(0,0)
        a5_d = self.a_Matrix(5,5)
        a5_0 = self.inv_a_Matrix(5,0)        
        a4_d = self.a_Matrix(4,4)
        a4_0 = self.inv_a_Matrix(4,0) 
        a3_d = self.a_Matrix(3,3)
        a3_0 = self.inv_a_Matrix(3,0)
        a2_d = self.a_Matrix(2,2)
        a2_0 = self.inv_a_Matrix(2,0)
        as_d = self.a_Matrix(1,1)
        as_0 = self.inv_a_Matrix(1,0)
        
        M_1 = self.M_Matrix(1)
        M_2 = self.M_Matrix(2)
        M_3 = self.M_Matrix(3)
        M_4 = self.M_Matrix(4)       
        
        # T_02 = multi_dot([a2_0,as_d,as_0,a0_0])
        # T_2inf = multi_dot([ainf_0,a4_d,a4_0,a3_d,a3_0,a2_d])
        # T_03 = multi_dot([a3_0,a2_d,a2_0,as_d,as_0,a0_0])
        # T_3inf = multi_dot([ainf_0,a4_d,a4_0,a3_d])
        # T_04 = multi_dot([a4_0,a3_d,a3_0,a2_d,a2_0,as_d,as_0,a0_0])
        # T_4inf = multi_dot([ainf_0,a4_d])  
        T_05 = multi_dot([a5_0,a4_d,a4_0,a3_d,a3_0,a2_d,a2_0,as_d,as_0,a0_0])
        T_5inf = multi_dot([ainf_0,a5_d])          
        
        T_02 = multi_dot([a2_0,M_1,a0_0])
        T_2inf = multi_dot([ainf_0,M_4,M_3,a2_d])
        T_03 = multi_dot([a3_0,M_2,M_1,a0_0])
        T_3inf = multi_dot([ainf_0,M_4,a3_d])
        T_04 = multi_dot([a4_0,M_3,M_2,M_1,a0_0])
        T_4inf = multi_dot([ainf_0,a4_d])        
        return T_05,T_5inf,T_04, T_4inf, T_03, T_3inf, T_02, T_2inf


    def E_field(self, material, position, fR, fL):
        ER = fR*np.exp(1j*self.omega*np.sqrt(self.epsilon[material]*self.mu)*self.z[position])
        EL = fL*np.exp(-1j*self.omega*np.sqrt(self.epsilon[material]*self.mu)*self.z[position])
        E_total = ER+EL
        return E_total
    
    def H_field(self,material,position,fR,fL):
        HR = fR*np.sqrt(self.epsilon[material]/self.mu)*np.exp(1j*self.omega*np.sqrt(self.epsilon[material]*self.mu)*self.z[position])
        HL = -fL*np.sqrt(self.epsilon[material]/self.mu)*np.exp(-1j*self.omega*np.sqrt(self.epsilon[material]*self.mu)*self.z[position])
        H_total = HR+HL
        return H_total
    
    
    def f_x(self,t,f_R,f_L):
        f_xR = t[0,0]*f_R+t[0,1]*f_L
        f_xL = t[1,0]*f_R+t[1,1]*f_L
        return f_xR,f_xL
        
        




class Change_thickness(object):
    def __init__(self,f_omega,omega,mu,epsilon):
        self.f_omega = f_omega
        self.omega = omega
        self.mu = mu
        self.epsilon = epsilon


    def heatloss(self,dnew,number_of_thickness):
        Q_heatloss1 = []
        Q_heatloss2 = []
        Q_heatloss3 = []
        Q_heatloss4 = []
        Q_heatloss5 = []
        Q_income = []
        for i in range(number_of_thickness):
            d      = dnew[i]
            matrix       = SpecialMatrix(self.omega,self.mu,self.epsilon,d)
            T_0s, T_sinf, T_0inf                     = matrix.Transfer_Matrix()
            T_05, T_5inf, T_04, T_4inf, T_03, T_3inf, T_02, T_2inf = matrix.Transfer_Matrix_special()            
            r        = matrix.Reflection_coeff(T_0inf)
            t        = matrix.Transmission_coeff(T_0inf)
            
            f_0R   = self.f_omega
            f_0L   = self.f_omega*r 
            
            # Income energy flux                          
            E_field_In = matrix.E_field(0,0,self.f_omega,0)
            H_field_In = matrix.H_field(0,0,self.f_omega,0)
            Q_in  = (pi*E_field_In.real*H_field_In.real)+(pi*E_field_In.imag*H_field_In.imag)
            

    
            #layer1 === substrate
            f_1R, f_1L = matrix.f_x(T_0s,f_0R, f_0L)
            E_field_01 = matrix.E_field(1,0,f_1R,f_1L)
            H_field_01 = matrix.H_field(1,0,f_1R,f_1L)
            E_field_d1 = matrix.E_field(1,1,f_1R,f_1L)
            H_field_d1 = matrix.H_field(1,1,f_1R,f_1L)
            Q1 = (pi*E_field_01.real*H_field_01.real)+(pi*E_field_01.imag*H_field_01.imag)-(pi*E_field_d1.real*H_field_d1.real)-(pi*E_field_d1.imag*H_field_d1.imag)
            # Q1 = heat_loss(self.omega, self.epsilon[1], self.mu, f_1R, f_1L,0,d[1])

      
            
            
    
            #layer2
            f_2R, f_2L = matrix.f_x(T_02,f_0R, f_0L)
            E_field_02 = matrix.E_field(2,0,f_2R,f_2L)
            H_field_02 = matrix.H_field(2,0,f_2R,f_2L)
            E_field_d2 = matrix.E_field(2,2,f_2R,f_2L)
            H_field_d2 = matrix.H_field(2,2,f_2R,f_2L)
            Q2 = (pi*E_field_02.real*H_field_02.real)+(pi*E_field_02.imag*H_field_02.imag)-(pi*E_field_d2.real*H_field_d2.real)-(pi*E_field_d2.imag*H_field_d2.imag)
            # Q2 = heat_loss(self.omega, self.epsilon[2], self.mu, f_2R, f_2L,0,d[2])




                        
            #layer3
            f_3R, f_3L = matrix.f_x(T_03,f_0R, f_0L)
            E_field_03 = matrix.E_field(3,0,f_3R,f_3L)
            H_field_03 = matrix.H_field(3,0,f_3R,f_3L)
            E_field_d3 = matrix.E_field(3,3,f_3R,f_3L)
            H_field_d3 = matrix.H_field(3,3,f_3R,f_3L)
            Q3 = (pi*E_field_03.real*H_field_03.real)+(pi*E_field_03.imag*H_field_03.imag)-(pi*E_field_d3.real*H_field_d3.real)-(pi*E_field_d3.imag*H_field_d3.imag)
            # Q3 = heat_loss(self.omega, self.epsilon[3], self.mu, f_3R, f_3L,0,d[3])
            

            
            
            
            #layer4
            f_4R, f_4L = matrix.f_x(T_04,f_0R, f_0L)
            E_field_04 = matrix.E_field(4,0,f_4R,f_4L)
            H_field_04 = matrix.H_field(4,0,f_4R,f_4L)
            E_field_d4 = matrix.E_field(4,4,f_4R,f_4L)
            H_field_d4 = matrix.H_field(4,4,f_4R,f_4L)
            Q4 = (pi*E_field_04.real*H_field_04.real)+(pi*E_field_04.imag*H_field_04.imag)-(pi*E_field_d4.real*H_field_d4.real)-(pi*E_field_d4.imag*H_field_d4.imag)
            # Q4 = heat_loss(self.omega, self.epsilon[4], self.mu, f_4R, f_4L,0,d[4])



            #layer5
            f_5R, f_5L = matrix.f_x(T_05,f_0R, f_0L)
            E_field_05 = matrix.E_field(5,0,f_5R,f_5L)
            H_field_05 = matrix.H_field(5,0,f_5R,f_5L)
            E_field_d5 = matrix.E_field(5,5,f_5R,f_5L)
            H_field_d5 = matrix.H_field(5,5,f_5R,f_5L)
            Q5 = (pi*E_field_05.real*H_field_05.real)+(pi*E_field_05.imag*H_field_05.imag)-(pi*E_field_d5.real*H_field_d5.real)-(pi*E_field_d5.imag*H_field_d5.imag)
            # Q4 = heat_loss(self.omega, self.epsilon[4], self.mu, f_4R, f_4L,0,d[4])

            
            Q_heatloss1.append(Q1)
            Q_heatloss2.append(Q2)
            Q_heatloss3.append(Q3)
            Q_heatloss4.append(Q4)
            Q_heatloss5.append(Q5)
            Q_income.append(Q_in)
        return Q_heatloss1,Q_heatloss2,Q_heatloss3,Q_heatloss4,Q_heatloss5,Q_income




    def absorption(self,dnew,number_of_thickness):
        Q_heatloss1,Q_heatloss2,Q_heatloss3,Q_heatloss4,Q_heatloss5,Q_income = self.heatloss(dnew,number_of_thickness)
        Q_total = np.array(Q_heatloss1) + np.array(Q_heatloss2) + np.array(Q_heatloss3) + np.array(Q_heatloss4) + np.array(Q_heatloss5)
        Absorption1 = np.array(Q_heatloss1)/np.array(Q_income)   #quartz
        Absorption2 = np.array(Q_heatloss2)/np.array(Q_income)   #NiFe
        Absorption3 = np.array(Q_heatloss3)/np.array(Q_income)   #Ru
        Absorption4 = np.array(Q_heatloss4)/np.array(Q_income)   #NiFe
        Absorption5 = np.array(Q_heatloss5)/np.array(Q_income)   #Al2O3
        Absorption_total = Q_total/Q_income
        return Absorption1,Absorption2,Absorption3,Absorption4,Absorption5,Absorption_total
        
    



    def reflection_heat(self,dnew,number_of_thickness):
        Q_reflect = []
        Q_income = []
        for i in range(number_of_thickness):
            d      = dnew[i]
            matrix       = SpecialMatrix(self.omega,self.mu,self.epsilon,d)
            T_0s, T_sinf, T_0inf                     = matrix.Transfer_Matrix()
            T_05, T_5inf, T_04, T_4inf, T_03, T_3inf, T_02, T_2inf = matrix.Transfer_Matrix_special()
            
            r        = matrix.Reflection_coeff(T_0inf)
            f_0L   = self.f_omega*r*0                            
            E_field_In = matrix.E_field(0,0,self.f_omega,f_0L)
            H_field_In = matrix.H_field(0,0,self.f_omega,f_0L)
            Q_in  = (pi*E_field_In.real*H_field_In.real)+(pi*E_field_In.imag*H_field_In.imag)
            
            f_L = self.f_omega*r
            E_ref = matrix.E_field(0,0,self.f_omega*0,f_L)
            H_ref = matrix.H_field(0,0,self.f_omega*0,f_L)
            Q_ref = -(pi*E_ref.real*H_ref.real)-(pi*E_ref.imag*H_ref.imag)
            
  
            Q_reflect.append(Q_ref)
            Q_income.append(Q_in)  
            
        return Q_reflect,Q_income



    def reflection(self,dnew,number_of_thickness):
        Q_reflect,Q_income = self.reflection_heat(dnew,number_of_thickness)
        reflection = np.array(Q_reflect)/np.array(Q_income)
        return reflection
        



    

    def transmission_heat(self,dnew,number_of_thickness):
        Q_transmission = []
        Q_income = []
        for i in range(number_of_thickness):
            d      = dnew[i]
            matrix       = SpecialMatrix(self.omega,self.mu,self.epsilon,d)
            T_0s, T_sinf, T_0inf                     = matrix.Transfer_Matrix()
            T_05, T_5inf, T_04, T_4inf, T_03, T_3inf, T_02, T_2inf = matrix.Transfer_Matrix_special()
            
            r        = matrix.Reflection_coeff(T_0inf)
            t        = matrix.Transmission_coeff(T_0inf)
            f_0L   = self.f_omega*r*0                            
            E_field_In = matrix.E_field(0,0,self.f_omega,f_0L)
            H_field_In = matrix.H_field(0,0,self.f_omega,f_0L)
            Q_in  = (pi*E_field_In.real*H_field_In.real)+(pi*E_field_In.imag*H_field_In.imag)
            

            f_T = self.f_omega*t
            E_trans = matrix.E_field(0,0,f_T,0)
            H_trans = matrix.H_field(0,0,f_T,0)
            Q_trans = (pi*E_trans.real*H_trans.real)+(pi*E_trans.imag*H_trans.imag)
            
            Q_transmission.append(Q_trans)
            Q_income.append(Q_in)            
    
        return Q_transmission,Q_income






    def transmission(self,dnew,number_of_thickness):
        Q_reflect,Q_income = self.reflection_heat(dnew,number_of_thickness)
        Q_trans,Q_income = self.transmission_heat(dnew,number_of_thickness)
        transmission = np.array(Q_trans)/np.array(Q_income)
        # transmission = (np.array(Q_income)-np.array(Q_trans)-np.array(Q_reflect))/np.array(Q_income)
        return transmission
    

        
    def trans_reflect(self,dnew,number_of_thickness):
        r_coeff = []
        t_coeff = []
        for i in range(number_of_thickness):
            d      = dnew[i]
            matrix       = SpecialMatrix(self.omega,self.mu,self.epsilon,d)
            T_0s, T_sinf, T_0inf                     = matrix.Transfer_Matrix()
            T_05, T_5inf, T_04, T_4inf, T_03, T_3inf, T_02, T_2inf = matrix.Transfer_Matrix_special()
            
            r        = matrix.Reflection_coeff(T_0inf)
            t        = matrix.Transmission_coeff(T_0inf)
            r_coeff.append(r)
            t_coeff.append(t)
            
        return r_coeff,t_coeff     
        




















































    