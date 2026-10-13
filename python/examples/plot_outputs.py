"""An EXAMPLE of plotting PACE's Python output. Not part of the package.

PACE for Python returns numbers; drawing them is left to you. These panels are
here to show what the bound functions give you, not to be a plotting API, and
they deliberately do not try to reproduce the R package's figures. The R plots
carry locked conventions -- which blocks appear, how pairs are normalised --
and a second implementation that drifted from them would be worse than none.

Run from the python/ directory:  python examples/plot_outputs.py
"""
import gzip
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import pace

FIX = Path(__file__).parent.parent / "tests" / "fixtures"
OUT = Path(__file__).parent


def col(name):
    return np.ascontiguousarray(np.loadtxt(FIX / (name + ".gz")).ravel())


# --- driver scores, the same Macrophage <- Tumour pair the R figure shows ---
genes = ["MRC1", "APOC1", "CD163", "AIF1", "FCER1G"]
ds = pace.driver_scores(
    estimate_shrunk=col("ds_estimate.txt"), spec=col("ds_spec.txt"),
    focal_mean=col("ds_focal_mean.txt"), mu_bar=col("ds_mu_bar.txt"),
    alpha=col("ds_alpha.txt"), var_n=float(np.loadtxt(FIX / "ds_var_n.txt.gz")))
order = np.argsort(ds["driver_score"])[::-1][:5]

# --- mash, the Dendritic_Cell slice ---
n_cond = int(np.loadtxt(FIX / "bc_Dendritic_Cell_dim.txt.gz"))
bhat = np.asfortranarray(np.loadtxt(FIX / "bc_Dendritic_Cell_b.txt.gz").reshape(-1, n_cond))
shat = np.asfortranarray(np.loadtxt(FIX / "bc_Dendritic_Cell_s.txt.gz").reshape(-1, n_cond))
m = pace.mash(bhat, shat)

fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))

ax[0].barh([genes[i] for i in order][::-1], ds["driver_score"][order][::-1], color="#F5A623")
ax[0].set_xlabel("Driver score (spatial cell state)")
ax[0].set_title("Macrophage <- Tumour", loc="left")

shrunk = m["posterior_mean"].ravel()
raw = bhat.ravel()
ax[1].scatter(raw, shrunk, s=6, alpha=0.3, color="#1f4e5f")
lim = np.percentile(np.abs(raw), 99.5)
ax[1].plot([-lim, lim], [-lim, lim], lw=0.8, color="grey", ls="--")
ax[1].set_xlim(-lim, lim); ax[1].set_ylim(-lim, lim)
ax[1].set_xlabel("raw slope"); ax[1].set_ylabel("shrunken")
ax[1].set_title("what mash does", loc="left")

called = m["lfsr"].ravel() < 0.05
# a linear count axis buries the calls: 12 of them against 2188 nulls piled at
# lfsr = 1. The log axis is what makes both visible at once.
ax[2].hist(np.log10(np.clip(m["lfsr"].ravel(), 1e-20, 1)), bins=50, color="#8c2f39")
ax[2].axvline(np.log10(0.05), color="black", lw=1, ls="--")
ax[2].set_yscale("log")
ax[2].set_xlabel("log10 lfsr"); ax[2].set_ylabel("gene-condition pairs (log)")
ax[2].set_title(f"{called.sum()} calls at lfsr < 0.05", loc="left")

for a in ax: a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "pace_python_output.png", dpi=110)
print("drew", OUT / "pace_python_output.png")
print("top driver:", genes[order[0]], f"{ds['driver_score'][order[0]]:.3f}")
print("mash calls:", int(called.sum()), "of", called.size)
