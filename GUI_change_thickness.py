#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 16:46:22 2023
Modified for GUI on March 27, 2025

@author: yingshuyang (original), Grok 3 (GUI adaptation)
"""

import numpy as np
import matplotlib.pyplot as plt
import absorption_calculation_class as cla2
from math import pi
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class AbsorptionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Absorption Calculator with Variable Thickness")
        
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
        table_frame.grid(row=1, column=0, columnspan=6, padx=5, pady=5)
        
        # Header row
        ttk.Label(table_frame, text="Layer").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        ttk.Label(table_frame, text="n").grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(table_frame, text="k").grid(row=0, column=2, padx=5, pady=5)
        ttk.Label(table_frame, text="d_min (m)").grid(row=0, column=3, padx=5, pady=5)
        ttk.Label(table_frame, text="d_max (m)").grid(row=0, column=4, padx=5, pady=5)
        
        # Layer inputs
        self.layers = ["Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5"]
        self.n_entries = {}
        self.k_entries = {}
        self.dmin_entries = {}
        self.dmax_entries = {}
        
        default_n = ["2.83", "3.13", "1", "1", "1"]
        default_k = ["0.1", "0.5", "0", "0", "0"]
        default_dmin = ["6e-9", "3e-9", "0.5e-3", "0", "0"]
        default_dmax = ["6e-9", "13e-9", "0.5e-3", "0", "0"]
        
        for i, layer in enumerate(self.layers):
            row = i + 1
            ttk.Label(table_frame, text=layer).grid(row=row, column=0, padx=5, pady=2, sticky="w")
            
            self.n_entries[layer] = ttk.Entry(table_frame, width=10)
            self.n_entries[layer].grid(row=row, column=1, padx=5, pady=2)
            self.n_entries[layer].insert(0, default_n[i])
            
            self.k_entries[layer] = ttk.Entry(table_frame, width=10)
            self.k_entries[layer].grid(row=row, column=2, padx=5, pady=2)
            self.k_entries[layer].insert(0, default_k[i])
            
            self.dmin_entries[layer] = ttk.Entry(table_frame, width=15)
            self.dmin_entries[layer].grid(row=row, column=3, padx=5, pady=2)
            self.dmin_entries[layer].insert(0, default_dmin[i])
            
            self.dmax_entries[layer] = ttk.Entry(table_frame, width=15)
            self.dmax_entries[layer].grid(row=row, column=4, padx=5, pady=2)
            self.dmax_entries[layer].insert(0, default_dmax[i])
        
        # Thickness changing layer and resolution
        ttk.Label(input_frame, text="Thickness Changing Layer (1-5):").grid(row=2, column=0, padx=5, pady=5, sticky="e")
        self.layer_var = tk.StringVar(value="2")
        ttk.OptionMenu(input_frame, self.layer_var, "2", "1", "2", "3", "4", "5").grid(row=2, column=1, padx=5, pady=5)
        
        ttk.Label(input_frame, text="Number of Thickness Points:").grid(row=3, column=0, padx=5, pady=5, sticky="e")
        self.num_entry = ttk.Entry(input_frame)
        self.num_entry.grid(row=3, column=1, padx=5, pady=5)
        self.num_entry.insert(0, "100")
        
        # Run button
        ttk.Button(root, text="Calculate", command=self.calculate).grid(row=1, column=0, pady=10)
        
        # Results frame
        self.result_frame = ttk.LabelFrame(root, text="Results")
        self.result_frame.grid(row=2, column=0, padx=10, pady=10, sticky="nsew")
        self.result_text = tk.Text(self.result_frame, height=15, width=60)  # Increased height for more data
        self.result_text.grid(row=0, column=0, padx=5, pady=5)
        
        # Plot frame (smaller size)
        self.plot_frame = ttk.LabelFrame(root, text="Absorption Plot")
        self.plot_frame.grid(row=0, column=1, rowspan=3, padx=10, pady=10, sticky="nsew")
        self.fig, self.ax = plt.subplots(figsize=(3, 2))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def calculate(self):
        try:
            # Get frequency
            frequency = float(self.freq_entry.get())
            omega_single = 2 * pi * frequency
            f_omega_single = 10
            
            # Fixed Air layer
            epsilon_optical = [cla2.epsilon(1, 0) * 8.85e-12]  # Air: n=1, k=0
            dmin_values = [0]  # Air: d=0
            dmax_values = [0]  # Air: d=0
            
            # Get dielectric constants and thicknesses
            for layer in self.layers:
                n = float(self.n_entries[layer].get())
                k = float(self.k_entries[layer].get())
                dmin = float(eval(self.dmin_entries[layer].get()))
                dmax = float(eval(self.dmax_entries[layer].get()))
                eps = cla2.epsilon(n, k) * 8.85e-12
                epsilon_optical.append(eps)
                dmin_values.append(dmin)
                dmax_values.append(dmax)
            
            mu = 12.57e-7  # permeability
            thickness_changing_layer = int(self.layer_var.get())
            number_of_thickness = int(self.num_entry.get())
            
            # Prepare thickness array
            d = [np.linspace(dmin_values[i], dmax_values[i], number_of_thickness) for i in range(6)]
            dnew = np.transpose(d)
            
            # Run calculations
            calc = cla2.Change_thickness(f_omega_single, omega_single, mu, epsilon_optical)
            Q_heatloss1, Q_heatloss2, Q_heatloss3, Q_heatloss4, Q_heatloss5, Q_income = calc.heatloss(dnew, number_of_thickness)
            Absorption1, Absorption2, Absorption3, Absorption4, Absorption5, Absorption_total = calc.absorption(dnew, number_of_thickness)
            reflection = calc.reflection(dnew, number_of_thickness)
            transmission = calc.transmission(dnew, number_of_thickness)
            Total = reflection + transmission + Absorption_total
            
            # Display all results
            self.result_text.delete(1.0, tk.END)
            thickness = d[thickness_changing_layer] * 1e9  # Convert to nm
            result_str = f"Results (for {number_of_thickness} points):\n\n"
            result_str += f"Thickness (nm):\n{', '.join([f'{t:.2f}' for t in thickness])}\n\n"
            result_str += f"Absorption_layer1:\n{', '.join([f'{a:.6f}' for a in Absorption1])}\n\n"
            result_str += f"Absorption_layer2:\n{', '.join([f'{a:.6f}' for a in Absorption2])}\n\n"
            result_str += f"Absorption_layer3:\n{', '.join([f'{a:.6f}' for a in Absorption3])}\n\n"
            result_str += f"Absorption_layer4:\n{', '.join([f'{a:.6f}' for a in Absorption4])}\n\n"
            result_str += f"Absorption_layer5:\n{', '.join([f'{a:.6f}' for a in Absorption5])}\n\n"
            result_str += f"Absorption_total:\n{', '.join([f'{a:.6f}' for a in Absorption_total])}\n\n"
            result_str += f"Transmission:\n{', '.join([f'{t:.6f}' for t in transmission])}\n\n"
            result_str += f"Reflection:\n{', '.join([f'{r:.6f}' for r in reflection])}\n\n"
            result_str += f"A+T+R:\n{', '.join([f'{t:.6f}' for t in Total])}\n"
            self.result_text.insert(tk.END, result_str)
            
            # Plot
            self.ax.clear()
            self.ax.plot(thickness, transmission, linewidth=2, label='Transmission')
            self.ax.plot(thickness, reflection, linewidth=2, label='Reflection')
            self.ax.plot(thickness, Absorption_total, linewidth=2, label='Absorption Total')
            self.ax.plot(thickness, Absorption1, '--', label='Absorption Layer 1')
            self.ax.plot(thickness, Absorption2, '--', label='Absorption Layer 2')
            self.ax.set_title('ART and Absorption Change', fontsize=10)
            self.ax.set_xlabel('Thickness (nm)', fontsize=8)
            self.ax.set_ylabel('Ratio', fontsize=8)
            self.ax.set_ylim(-0.1, 1.1)
            self.ax.legend(fontsize=6)
            self.ax.tick_params(axis='both', labelsize=6)
            self.fig.tight_layout()
            self.canvas.draw()
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")

def main():
    root = tk.Tk()
    app = AbsorptionGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()