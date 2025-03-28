#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:46:22 2023
Modified for GUI on March 27, 2025

@author: yingshuyang (original), Grok 3 (GUI adaptation)
"""

import numpy as np
import absorption_calculation_class as cla2
from math import pi
import tkinter as tk
from tkinter import ttk, messagebox

class AbsorptionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Absorption Calculator")
        
        # Frame for inputs
        input_frame = ttk.LabelFrame(root, text="Input Parameters")
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        
        # Frequency
        ttk.Label(input_frame, text="Frequency (Hz):").grid(row=0, column=0, padx=5, pady=5, sticky="e")
        self.freq_entry = ttk.Entry(input_frame)
        self.freq_entry.grid(row=0, column=1, padx=5, pady=5, columnspan=2)
        self.freq_entry.insert(0, "374.74e12")
        
        # Table for layers (excluding Air)
        table_frame = ttk.Frame(input_frame)
        table_frame.grid(row=1, column=0, columnspan=4, padx=5, pady=5)
        
        # Header row
        ttk.Label(table_frame, text="Layer").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        ttk.Label(table_frame, text="n").grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(table_frame, text="k").grid(row=0, column=2, padx=5, pady=5)
        ttk.Label(table_frame, text="d (m)").grid(row=0, column=3, padx=5, pady=5)
        
        # Layer inputs (excluding Air)
        self.layers = ["Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5"]
        self.n_entries = {}
        self.k_entries = {}
        self.d_entries = {}
        
        default_n = ["3.13", "4.9", "1.76", "1", "1"]
        default_k = ["4.76", "5.48", "0", "0", "0"]
        default_d = ["3e-9", "6e-9", "0.5e-3", "0", "0"]
        
        for i, layer in enumerate(self.layers):
            row = i + 1
            ttk.Label(table_frame, text=layer).grid(row=row, column=0, padx=5, pady=2, sticky="w")
            
            self.n_entries[layer] = ttk.Entry(table_frame, width=10)
            self.n_entries[layer].grid(row=row, column=1, padx=5, pady=2)
            self.n_entries[layer].insert(0, default_n[i])
            
            self.k_entries[layer] = ttk.Entry(table_frame, width=10)
            self.k_entries[layer].grid(row=row, column=2, padx=5, pady=2)
            self.k_entries[layer].insert(0, default_k[i])
            
            self.d_entries[layer] = ttk.Entry(table_frame, width=15)
            self.d_entries[layer].grid(row=row, column=3, padx=5, pady=2)
            self.d_entries[layer].insert(0, default_d[i])
        
        # Run button
        ttk.Button(root, text="Calculate", command=self.calculate).grid(row=1, column=0, pady=10)
        
        # Results frame
        self.result_frame = ttk.LabelFrame(root, text="Results")
        self.result_frame.grid(row=2, column=0, padx=10, pady=10, sticky="nsew")
        self.result_text = tk.Text(self.result_frame, height=10, width=60)
        self.result_text.grid(row=0, column=0, padx=5, pady=5)
        
    def calculate(self):
        try:
            # Get frequency
            frequency = float(self.freq_entry.get())  # Already in Hz
            omega_single = 2 * pi * frequency
            f_omega_single = 10  # Fixed as per original code
            
            # Fixed Air layer
            epsilon_optical = [cla2.epsilon(1, 0) * 8.85e-12]  # Air: n=1, k=0
            thicknesses = [0]  # Air: d=0
            
            # Get dielectric constants and thicknesses for editable layers
            for layer in self.layers:
                n = float(self.n_entries[layer].get())
                k = float(self.k_entries[layer].get())
                d = float(eval(self.d_entries[layer].get()))  # eval to handle scientific notation
                eps = cla2.epsilon(n, k) * 8.85e-12
                epsilon_optical.append(eps)
                thicknesses.append(d)
            
            mu = 12.57e-7  # permeability
            
            # Prepare thickness array
            number_of_thickness = 1
            dnew = np.transpose([np.linspace(t, t, number_of_thickness) for t in thicknesses])
            
            # Run calculations
            calc = cla2.Change_thickness(f_omega_single, omega_single, mu, epsilon_optical)
            Q_heatloss1, Q_heatloss2, Q_heatloss3, Q_heatloss4, Q_heatloss5, Q_income = calc.heatloss(dnew, number_of_thickness)
            Absorption1, Absorption2, Absorption3, Absorption4, Absorption5, Absorption_total = calc.absorption(dnew, number_of_thickness)
            reflection = calc.reflection(dnew, number_of_thickness)
            transmission = calc.transmission(dnew, number_of_thickness)
            Total = reflection + transmission + Absorption_total
            
            # Display results
            self.result_text.delete(1.0, tk.END)
            result_str = (
                f"Results:\n"
                f"Absorption_layer1 = {Absorption1[0]:.6f}\n"
                f"Absorption_layer2 = {Absorption2[0]:.6f}\n"
                f"Absorption_layer3 = {Absorption3[0]:.6f}\n"
                f"Absorption_layer4 = {Absorption4[0]:.6f}\n"
                f"Absorption_layer5 = {Absorption5[0]:.6f}\n"
                f"Absorption_total = {Absorption_total[0]:.6f}\n"
                f"Transmission      = {transmission[0]:.6f}\n"
                f"Reflection        = {reflection[0]:.6f}\n"
                f"A+T+R             = {Total[0]:.6f}"
            )
            self.result_text.insert(tk.END, result_str)
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")

def main():
    root = tk.Tk()
    app = AbsorptionGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()