#!/usr/bin/env python
# coding: utf-8

# ## Transition Energies of the eleven C15 Defect Cluster in Tungsten (Ramp–Arrhenius Fitting)

# In[ ]:


import json
import numpy as np

# -----------------------------
# READ GROUPED JSON
# -----------------------------
with open("combined_rampflip_data.json", "r") as f:
    grouped_data = json.load(f)

# -----------------------------
# ARRHENIUS PARAMETERS
# -----------------------------
k_B = 1.38e-23          # J/K
eV_to_J = 1.602e-19
r = (2100 - 300) / 1e-07
w0_list = [1e12, 1e13]

# -----------------------------
# PRINT HEADER
# -----------------------------
header = (
    f"{'JSON File':<20} "
    f"{'E_t (eV) @1e12':>15} {'ΔE_t (eV) @1e12':>18} "
    f"{'E_t (eV) @1e13':>15} {'ΔE_t (eV) @1e13':>18} "
    f"{'Avg Temp (K)':>12}"
)
print(header)
print("-" * len(header))

# -----------------------------
# CALCULATE Eₜ FOR EACH JSON
# -----------------------------
for json_file, entries in grouped_data.items():
    # Extract valid flip temperatures
    temps = []
    for e in entries:
        try:
            T = float(e["flipTemp"])
            if 300 < T < 2100:
                temps.append(T)
        except:
            pass

    if not temps:
        print(f"{json_file:<20} {'nan':>15} {'nan':>18} {'nan':>15} {'nan':>18} {'nan':>12} ")
        continue

    temps = np.array(temps)
    avg_temp = np.mean(temps)
    delta_T = np.std(temps)

    results = []
    for w_0 in w0_list:
        E_t = -k_B * avg_temp * (np.log(r / w_0) - np.log(avg_temp))
        dE_t_dT = -k_B * (np.log(r / w_0) - np.log(avg_temp) - 1)
        delta_E_t = np.abs(dE_t_dT) * delta_T

        E_t_eV = E_t / eV_to_J
        delta_E_t_eV = delta_E_t / eV_to_J
        results.extend([E_t_eV, delta_E_t_eV])

    print(f"{json_file:<20} {results[0]:15.3f} {results[1]:18.3f} {results[2]:15.3f} {results[3]:18.3f} {avg_temp:12.2f} ")


