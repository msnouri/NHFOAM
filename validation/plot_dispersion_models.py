#!/usr/bin/env python3
"""Standing-wave period error and damping vs kh: NHFOAM (pressure on the
layer interfaces), NHWAVE and NHFOAM with cell-centred pressure, with the
same basin, resolution (64 cells per wavelength, T/dt = 200) and number of
uniform layers.

usage: plot_dispersion_models.py <interface.csv> <nhwave_dir> <out.png> [<cell.csv>]
    csv files: from run_dispersion_study.sh (PRESSURE=interface or cell)
               (kh,layers,T_lin,T_nh,ratio,amp_ratio)
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
nc = (np.genfromtxt(sys.argv[4], delimiter=",", names=True)
      if len(sys.argv) > 4 else None)
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

fig, axes = plt.subplots(1, 2, figsize=(10, 4.1))
cols = {2: "#8e44ad", 3: "#c0392b", 10: "#1f5fa8"}


def series(tab, nl, key):
    m = tab["layers"] == nl
    o = np.argsort(tab["kh"][m])
    return tab["kh"][m][o], tab[key][m][o]


for nl, col in cols.items():
    kh, r = series(nf, nl, "ratio")
    axes[0].semilogy(kh, 100*np.abs(r - 1), "o-", color=col, ms=4,
                     label=f"NHFOAM, {nl} layers")
    kh, a = series(nf, nl, "amp_ratio")
    axes[1].plot(kh, 100*(1 - a), "o-", color=col, ms=4)

    k = [x for x in khs if (x, nl) in nw]
    axes[0].semilogy(k, [100*abs(nw[(x, nl)][0] - 1) for x in k], "s--",
                     color=col, ms=4, mfc="none", label=f"NHWAVE, {nl} layers")
    axes[1].plot(k, [100*(1 - nw[(x, nl)][1]) for x in k], "s--", color=col,
                 ms=4, mfc="none")

    if nc is not None and nl != 2:
        kh, r = series(nc, nl, "ratio")
        axes[0].semilogy(kh, 100*np.abs(r - 1), "^:", color=col, ms=4,
                         alpha=0.7,
                         label=f"NHFOAM cell-centred, {nl} layers")
        kh, a = series(nc, nl, "amp_ratio")
        axes[1].plot(kh, 100*(1 - a), "^:", color=col, ms=4, alpha=0.7)

axes[0].set_ylabel(r"$|$period error$|$ [\%]")
axes[1].set_ylabel(r"amplitude loss after 4.5 periods [\%]")
axes[0].set_ylim(0.01, 60)
for ax in axes:
    ax.set_xlabel(r"$kh$")
    ax.grid(alpha=0.3, which="both")
axes[0].legend(fontsize=7, ncol=1, loc="upper left")
fig.suptitle(r"Standing waves, 64 cells per wavelength, $T/\Delta t = 200$,"
             r" uniform layers")
fig.tight_layout()
save(fig, out, dpi=150)
