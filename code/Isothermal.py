#!/usr/bin/env python
# coding: utf-8

# ## Isothermal Simulation–Based Transition Energy Analysis for 2C, TRIPOD-3B, and 5C Clusters Using Arrhenius Fitting
print("\n==================== Arrhenius Fitting ====================\n")
# 

# In[4]:


import json
import numpy as np
from scipy import stats

# ---------------- CONFIG ----------------
json_file = "combined_transflip_data.json"
kB = 8.617333262e-5  # eV/K
MIN_POINTS = 3
TAU_MAX = 1e6      
# --------------------------------------

with open(json_file, "r") as f:
    combined = json.load(f)

def fit_cluster(cluster_name, entries):

    invT = []
    ln_tau = []
    skipped = 0

    for e in entries:
        try:
            T = float(e.get("flipTemp", e.get("temp")))
            tau = float(e["flipTime"])

            # ---------- VALID DATA CONDITION ----------
            if (
                T > 0 and
                tau > 0 and
                tau < TAU_MAX and
                np.isfinite(T) and
                np.isfinite(tau)
            ):
                invT.append(1.0 / T)
                ln_tau.append(np.log(tau))
            else:
                skipped += 1

        except Exception:
            skipped += 1

    invT = np.array(invT)
    ln_tau = np.array(ln_tau)

    if len(invT) < MIN_POINTS:
        print(f"⚠️ {cluster_name}: not enough valid points ({len(invT)})")
        return

    # -------- Linear fit --------
    slope, intercept, r, p, stderr = stats.linregress(invT, ln_tau)

    # -------- Arrhenius params --------
    Ea = slope * kB
    Ea_err = stderr * kB
    tau0 = np.exp(intercept)
    w0 = 1.0 / tau0

    # -------- Output --------
    print(f"\n================ {cluster_name} =================")
    print(f"Valid points   : {len(invT)}")
    print(f"Skipped points : {skipped}")
    print(f"R²             : {r**2:.4f}")
    print(f"Eₐ             : {Ea:.6f} ± {Ea_err:.6f} eV")
    print(f"τ₀             : {tau0:.3e} s")
    print(f"ω₀             : {w0:.3e} s⁻¹")

# ---------------- RUN SEPARATELY ----------------
fit_cluster("2C", combined["2C"])
fit_cluster("Tripod-3B", combined["3T"])
fit_cluster("5C", combined["5C"])


# #### Eyring Analysis of Enthalpy and Entropy for 2C, TRIPOD-3B, and 5C Clusters
print("\n==================== Eyring Fit ====================\n")

# In[7]:



import json
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# -------- Constants --------
kB = 8.617333262e-5   # eV/K
h  = 4.135667696e-15  # eV·s
TAU_MAX = 9.9e6
MIN_POINTS = 3

# -------- Display names & colors --------
display_names = {
    "2C": "2C",
    "3T": "Tripod-3B",
    "5C": "5C"
}

colors = {
    "2C": "blue",
    "3T": "green",
    "5C": "tab:blue"
}

# -------- Load JSON --------
with open("combined_transflip_data.json", "r") as f:
    combined = json.load(f)

# -------- Loop over clusters --------
for cluster_key, data in combined.items():
    cluster_name = display_names.get(cluster_key, cluster_key)

    x = []
    y = []
    skipped = 0

    for entry in data:
        try:
            T = float(entry.get("temp", entry.get("flipTemp")))
            tau = float(entry["flipTime"])

            if 0 < tau < TAU_MAX and np.isfinite(tau):
                omega = 1.0 / tau
                x.append(1.0 / T)
                y.append(np.log(omega / T))
            else:
                skipped += 1

        except Exception:
            skipped += 1

    x = np.array(x)
    y = np.array(y)

    if len(x) < MIN_POINTS:
        print(f"⚠️ {cluster_name}: Not enough valid points")
        continue

    # -------- Linear fit --------
    slope, intercept, r, p, stderr = stats.linregress(x, y)

    # -------- Eyring parameters --------
    Delta_H = -kB * slope
    Delta_H_err = kB * stderr
    Delta_S = kB * (intercept - np.log(kB / h))

    # -------- Print results --------
    print(f"\n=========== {cluster_name} ===========")
    print(f"Valid points : {len(x)}")
    print(f"Skipped      : {skipped}")
    print(f"R²           : {r**2:.4f}")
    print(f"ΔH‡          : {Delta_H:.6f} ± {Delta_H_err:.6f} eV")
    print(f"ΔS‡          : {Delta_S:.3e} eV/K "
          f"({Delta_S * 96485:.2f} J/mol·K)")

    # -------- Plot --------
    plt.figure(figsize=(5.5, 4.5), dpi=600)

    plt.scatter(
        x,
        y,
        s=15,
        color=colors[cluster_key],
        label=f"Data-{cluster_name}",
        zorder=2
    )

    # 🔧 ONLY FIX: extend fit slightly so last point is not cut
    xfit = np.linspace(x.min(), x.max() * 1.01, 300)

    plt.plot(
        xfit,
        slope * xfit + intercept,
        color="red",
        linewidth=1.6,
        label="Eyring fit",
        zorder=3
    )
    dx = x.max() - x.min()
    plt.xlim(x.min() - 0.02 * dx, x.max() + 0.02 * dx)
    

    plt.xlabel("1 / T (1/K)", fontsize=9)
    plt.ylabel("ln(ω / T)", fontsize=9)
    plt.legend(fontsize=7, frameon=False)
    plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    plt.tight_layout()

    plt.savefig(f"eyring_fit_{cluster_name}.jpg", dpi=600, bbox_inches="tight")
    plt.savefig(f"eyring_fit_{cluster_name}.pdf", bbox_inches="tight")
    plt.show()
# ## Eyring-Based Transition Times of 2C, TRIPOD-3B, and 5C Clusters at Various Temperatures
# 

# In[3]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================================================
#                    CONSTANTS
# =========================================================
kB = 8.617e-5         # eV/K
h  = 4.135667696e-15  # eV·s

# =========================================================
#                    TEMPERATURES
# =========================================================
T_list = [77, 194.65, 300, 500, 900]

# Custom x spacing ONLY for visualization
x_positions = [0, 2.5, 4, 5.5, 7]
x_labels = ["77 (LN₂)", "194.65 (Dry Ice)", "300", "500", "900"]

# =========================================================
#               EYRING PARAMETERS (FROM FITS)
# =========================================================
clusters_info = {
    "2C":         {"Em": 0.53, "dS": -4.00e-4},
    "Tripod-3B":  {"Em": 0.11, "dS": -7.90e-4},
    "5C":         {"Em": 0.58, "dS": -7.844e-4}
}

# =========================================================
#          COMPUTE RATES + TRANSITION TIMES (ONCE)
# =========================================================
records = []

for cluster, params in clusters_info.items():
    Em, dS = params["Em"], params["dS"]
    for T in T_list:
        w = (kB * T / h) * np.exp(dS / kB) * np.exp(-Em / (kB * T))
        t_s = 1.0 / w
        t_days = t_s / 86400.0

        records.append([cluster, T, w, t_s, t_days])

df = pd.DataFrame(
    records,
    columns=["Cluster", "T (K)", "w (1/s)", "t (s)", "t (days)"]
)

# =========================================================
#            MANUSCRIPT-STYLE NUMERICAL OUTPUT
# =========================================================
print("\n==================== Transition Times from Eyring Fits ====================\n")

for cluster in df["Cluster"].unique():
    print(f"--- {cluster} ---")
    sub = df[df["Cluster"] == cluster]

    for _, row in sub.iterrows():
        print(
            f"T = {row['T (K)']:>7.2f} K   "
            f"w = {row['w (1/s)']:.3e}   "
            f"t = {row['t (s)']:.3e} s   "
            f"t = {row['t (days)']:.3e} days"
        )
    print()

print("===========================================================================\n")

# =========================================================
#                       PLOT
# =========================================================
plt.figure(figsize=(7, 4))

colors  = {"2C": "#0b1d51", "Tripod-3B": "#9f2b68", "5C": "#237a57"}
markers = {"2C": "o",       "Tripod-3B": "s",       "5C": "^"}

for cluster in clusters_info.keys():
    subset = df[df["Cluster"] == cluster]

    plt.plot(
        x_positions,
        subset["t (days)"],
        label=cluster,
        color=colors[cluster],
        marker=markers[cluster],
        linewidth=2.2,
        markersize=6
    )

plt.yscale("log")
plt.xticks(x_positions, x_labels)

plt.xlabel("Temperature (K)", fontsize=12)
plt.ylabel("Transition Time (days)", fontsize=12)
plt.title("Cluster Transition Time vs Temperature", fontsize=14)

plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend(frameon=False)
plt.tight_layout()

plt.savefig("cluster_transition_times_eyring.jpg", dpi=600)
plt.show()


# In[ ]:





# In[ ]:




