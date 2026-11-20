#!/usr/bin/env python3
"""Wave height (zero-crossing mean over the last periods) at each gauge."""
import sys, numpy as np
case = sys.argv[1]; T = float(sys.argv[2]); nPer = int(sys.argv[3]) if len(sys.argv) > 3 else 5
fn = f"{case}/postProcessing/gauges/0/eta.dat"
d = np.loadtxt(fn, comments="#"); t = d[:, 0]
names = [l.split()[1] for l in open(fn) if l.startswith("# G")]
m = t > t[-1] - nPer * T
out = []
for i, n in enumerate(names):
    e = d[m, i + 1]
    out.append(f"{n}:{(e.max() - e.min()):.4f}")
print(" ".join(out))
