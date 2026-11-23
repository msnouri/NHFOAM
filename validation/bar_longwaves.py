#!/usr/bin/env python3
"""Submerged bar: the slow drift of a run (long waves of the flume).

usage: bar_longwaves.py <T> [kind:]<case>[:label] ...
For windows of 5 periods ending at t = 45, 50, ..., 65 s, prints the time
lag of the first harmonic at x = 35.7 m behind x = 22 m and the spread of the
mean water level over the gauges (max - min of the window means).
kind = nhfoam (default) | nhwave | reef3d (see other_models.py).
"""
import sys
import numpy as np
import other_models as om

T = float(sys.argv[1])
runs = [om.parse_case(a) for a in sys.argv[2:]]
xs = [22, 24, 30.5, 32.5, 33.5, 34.5, 35.7, 37.3, 39.0, 41.0]
w = 2 * np.pi / T
ends = [45.0, 50.0, 55.0, 60.0, 65.0]
print(f"{'run':38s}" + "".join(f"  t={e:g}: lag  spread" for e in ends))
for kind, path, label in runs:
    t, eta = om.read_any(kind, path, xs)
    row = f"{label[:38]:38s}"
    lags = []
    for e in ends:
        if t[-1] < e - 1e-6:
            row += "        --       --"
            continue
        m = (t > e - 5 * T) & (t <= e)
        A = {x: 2 * np.mean(eta[x][m] * np.exp(1j * w * t[m])) for x in (22, 35.7)}
        lag = (np.angle(A[35.7] / A[22]) % (2 * np.pi)) / w
        if lags:  # unwrap across the period
            lag += T * np.round((lags[-1] - lag) / T)
        lags.append(lag)
        means = [eta[x][m].mean() for x in xs]
        row += f"  {lag:8.3f} s {1000*(max(means)-min(means)):5.2f} mm"
    if len(lags) > 1:
        row += f"   drift {lags[-1] - lags[0]:+.3f} s"
    print(row)
