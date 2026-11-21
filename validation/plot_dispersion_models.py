#!/usr/bin/env python3
"""Standing-wave period error and damping vs kh: NHFOAM and NHWAVE with the
same basin, resolution (64 cells per wavelength, T/dt = 200) and number of
uniform layers.

usage: plot_dispersion_models.py <nhfoam_dispersion.csv> <nhwave_dir> <out.png>
    nhfoam csv: from run_dispersion_study.sh (kh,layers,T_lin,T_nh,ratio,amp_ratio)
    nhwave_dir: runs of other_models.nhwave_standing named kh<kh>_L<n>
(also writes out.pdf and prints the NHWAVE table)
"""
import sys, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from latex_style import use_latex, save
from analyse_standing import fit

use_latex(10)
nf = np.genfromtxt(sys.argv[1], delimiter=",", names=True)
nwdir, out = sys.argv[2], sys.argv[3]
g = 9.81
khs = [0.5, 1, 2, 3, 5, 7]
layers = [1, 2, 3, 5, 10, 20]

nw = {}
print("NHWAVE: kh, layers, period ratio, amplitude ratio")
for kh in khs:
    for nl in layers:
        fn = os.path.join(nwdir, f"kh{kh:g}_L{nl}", "results", "probe_0001")
        if not os.path.exists(fn):
            continue
        d = np.loadtxt(fn)
        T = 2 * np.pi / np.sqrt(g * kh * np.tanh(kh))
        Tzc, peaks = fit(d[:, 0], d[:, 1], T)
        amp = peaks[-2, 1] / peaks[0, 1] if len(peaks) > 2 else np.nan
        nw[(kh, nl)] = (Tzc / T, amp)
        print(f"{kh:4g} {nl:3d} {Tzc / T:.5f} {amp:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(10, 3.9))
cols = {3: "#c0392b", 5: "#e67e22", 10: "#1f5fa8", 20: "#2e8b57"}
for nl, col in cols.items():
    m = nf["layers"] == nl
    o = np.argsort(nf["kh"][m])
    axes[0].plot(nf["kh"][m][o], 100 * (nf["ratio"][m][o] - 1), "o-", color=col, ms=4,
                 label=f"NHFOAM, {nl} layers")
    axes[1].plot(nf["kh"][m][o], 100 * (1 - nf["amp_ratio"][m][o]), "o-", color=col, ms=4)
    k = [kh for kh in khs if (kh, nl) in nw]
    axes[0].plot(k, [100 * (nw[(kh, nl)][0] - 1) for kh in k], "s--", color=col, ms=4,
                 mfc="none", label=f"NHWAVE, {nl} layers")
    axes[1].plot(k, [100 * (1 - nw[(kh, nl)][1]) for kh in k], "s--", color=col, ms=4, mfc="none")
axes[0].axhline(0, color="k", lw=0.6)
axes[0].set_ylabel(r"period error [\%]")
axes[1].set_ylabel(r"amplitude loss after 4.5 periods [\%]")
for ax in axes:
    ax.set_xlabel(r"$kh$")
    ax.grid(alpha=0.3)
axes[0].set_ylim(-3, 8)
axes[0].legend(fontsize=7.5, ncol=2)
fig.suptitle(r"Standing waves, 64 cells per wavelength, $T/\Delta t = 200$, uniform layers")
fig.tight_layout()
save(fig, out, dpi=150)
