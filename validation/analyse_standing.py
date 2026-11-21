#!/usr/bin/env python3
"""Period and amplitude decay of a standing wave from NHFOAM gauge output."""
import sys, numpy as np

def load(case):
    d = np.loadtxt(f"{case}/postProcessing/gauges/0/eta.dat", comments="#")
    return d[:, 0], d[:, 1]

def fit(t, e, T_guess):
    """Least-squares fit e = A exp(-gamma t) cos(omega t) near the expected frequency."""
    from numpy.linalg import lstsq
    # zero up-crossings (and down) -> period
    s = np.sign(e)
    idx = np.where(s[:-1] * s[1:] < 0)[0]
    tz = t[idx] - e[idx] * (t[idx + 1] - t[idx]) / (e[idx + 1] - e[idx])
    if len(tz) >= 3:
        T_zc = 2 * np.mean(np.diff(tz))
    else:
        T_zc = np.nan
    # peaks (abs max per half period) -> damping
    peaks = []
    bounds = np.concatenate(([t[0]], tz, [t[-1]]))
    for a, b in zip(bounds[:-1], bounds[1:]):
        m = (t >= a) & (t <= b)
        if m.sum() > 3:
            peaks.append((t[m][np.argmax(np.abs(e[m]))], np.max(np.abs(e[m]))))
    peaks = np.array(peaks)
    return T_zc, peaks

if __name__ == "__main__":
    case = sys.argv[1]; T_lin = float(sys.argv[2])
    t, e = load(case)
    e0 = e[0] if abs(e[0]) > 0 else 1
    Tzc, peaks = fit(t, e, T_lin)
    print(f"period NHFOAM = {Tzc:.5f} s, linear theory = {T_lin:.5f} s, ratio = {Tzc/T_lin:.5f}")
    if len(peaks) > 2:
        # amplitude at last full half period relative to first
        nT = (peaks[-2, 0] - peaks[0, 0]) / T_lin
        print(f"amplitude ratio over {nT:.2f} periods: {peaks[-2,1]/peaks[0,1]:.5f}")
