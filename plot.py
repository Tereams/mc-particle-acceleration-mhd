import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read Data
df = pd.read_csv("results.csv")

# =========================
# 1. Energy spectrum
# =========================
W = df["W_final_eV"].values
W = W[W > 0]

plt.figure()
bins = np.logspace(np.log10(W.min()), np.log10(W.max()), 50)
plt.hist(W, bins=bins, histtype="step", density=True)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("W [eV]")
plt.ylabel("PDF")
plt.title("Energy spectrum")
plt.tight_layout()
plt.savefig("energy_spectrum.png", dpi=200)

# =========================
# 2. Escape Time
# =========================
df_esc = df[df["escaped"] == 1]
t = df_esc["t_final"].values

plt.figure()
bins = np.logspace(np.log10(t.min()), np.log10(t.max()), 50)
plt.hist(t, bins=bins, histtype="step", density=True)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Escape time [s]")
plt.ylabel("PDF")
plt.title("Escape time distribution")
plt.tight_layout()
plt.savefig("escape_time.png", dpi=200)

# =========================
# 3. Energy vs kicks
# =========================
plt.figure()
plt.scatter(df["kkicks"], df["W_final_eV"], s=3, alpha=0.3)
plt.yscale("log")
plt.xlabel("Number of kicks")
plt.ylabel("Final energy [eV]")
plt.title("Energy vs kicks")
plt.tight_layout()
plt.savefig("energy_vs_kicks.png", dpi=200)

plt.show()
