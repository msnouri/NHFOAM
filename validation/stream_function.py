#!/usr/bin/env python3
"""Steady periodic waves by the Fourier approximation method of Fenton (1988),
for the stream-function wave theory of NHFOAM's regularWave (forcing zones and
nhSetWaves).

In a frame moving with the wave (X = x - c t, z above the bed)

    psi(X, z) = -U z + sum_j B_j sinh(j k z)/cosh(j k d) cos(j k X),

with psi = -Q and (u^2 + w^2)/2 + g eta = R on the surface z = eta(X).  The
unknowns are k, c, U, Q, R, eta at N+1 points over half a wave and B_1..B_N;
the closing conditions are the mean depth d, the height H, the period T and
the mean current: zero mass flux (current='mass', a closed flume) or zero
Eulerian current (current='euler').  The height is raised in steps from a
small value (Newton with a numerical Jacobian at each step).

In the fixed frame
    eta(x, t) = sum_j E_j cos(j theta),          theta = k x - omega t
    u = uE + sum_j j k B_j cosh(j k z)/cosh(j k d) cos(j theta)
    w =      sum_j j k B_j sinh(j k z)/cosh(j k d) sin(j theta)
with uE = c - U the Eulerian mean current and E_j the cosine coefficients of
the surface relative to the mean level.

usage: stream_function.py H T d [N]       prints the NHFOAM dictionary entries

Fenton, J.D. (1988) The numerical solution of steady water wave problems.
Computers & Geosciences 14(3), 357-368.
"""
import sys
import numpy as np

G = 9.81


def _residual(x, N, H, T, current):
    k = x[0]
    eta = x[1:N + 2]
    B = x[N + 2:2 * N + 2]
    U, Q, R, c = x[2 * N + 2:2 * N + 6]
    j = np.arange(1, N + 1)
    m = np.arange(N + 1)
    th = np.outer(m, j) * np.pi / N                     # (N+1, N)
    kz = k * np.outer(eta, j)                            # (N+1, N)
    sh = np.sinh(kz) / np.cosh(j * k)
    ch = np.cosh(kz) / np.cosh(j * k)
    psi = -U * eta + (sh * np.cos(th)) @ B
    u = -U + (ch * np.cos(th)) @ (j * k * B)
    w = (sh * np.sin(th)) @ (j * k * B)
    f = np.empty(2 * N + 6)
    f[:N + 1] = psi + Q
    f[N + 1:2 * N + 2] = 0.5 * (u ** 2 + w ** 2) + eta - R
    f[2 * N + 2] = (0.5 * eta[0] + eta[1:N].sum() + 0.5 * eta[N]) / N - 1.0
    f[2 * N + 3] = eta[0] - eta[N] - H
    f[2 * N + 4] = k * c * T - 2 * np.pi
    f[2 * N + 5] = (c - Q) if current == 'mass' else (c - U)
    return f


def _newton(x, N, H, T, current, tol=1e-12, maxit=50):
    for it in range(maxit):
        f = _residual(x, N, H, T, current)
        if np.max(np.abs(f)) < tol:
            return x
        J = np.empty((len(x), len(x)))
        for i in range(len(x)):
            dx = 1e-7 * max(1.0, abs(x[i]))
            xp = x.copy(); xp[i] += dx
            J[:, i] = (_residual(xp, N, H, T, current) - f) / dx
        x = x - np.linalg.solve(J, f)
    raise RuntimeError('stream function: Newton did not converge '
                       '(max residual %.2e)' % np.max(np.abs(f)))


def fenton(H, T, d, N=20, current='mass', g=G, steps=10):
    """Dimensional solution: dict with k, c, omega, uE, E (0..N), B (1..N),
    eta(x, t), velocity(x, z, t) (z above the still-water level)."""
    # nondimensional with d and g
    Hn, Tn = H / d, T * np.sqrt(g / d)
    om = 2 * np.pi / Tn
    k = om ** 2
    for _ in range(100):                                 # linear dispersion
        k -= (k * np.tanh(k) - om ** 2) / (np.tanh(k) + k / np.cosh(k) ** 2)
    c = om / k
    m = np.arange(N + 1)
    x = np.zeros(2 * N + 6)
    h0 = Hn / steps
    x[0] = k
    x[1:N + 2] = 1 + 0.5 * h0 * np.cos(m * np.pi / N)
    x[N + 2] = 0.5 * h0 * c / np.tanh(k)
    x[2 * N + 2:2 * N + 6] = [c, c, 0.5 * c ** 2 + 1, c]
    sols = []
    for s in range(1, steps + 1):
        Hs = Hn * s / steps
        if len(sols) >= 2:                               # extrapolate in H
            x = 2 * sols[-1] - sols[-2]
        elif sols:
            x = sols[-1]
        sols.append(_newton(x, N, Hs, Tn, current))
    x = sols[-1]
    k, etan, Bn = x[0], x[1:N + 2], x[N + 2:2 * N + 2]
    U, Q, R, c = x[2 * N + 2:2 * N + 6]
    # cosine coefficients of eta - d (DCT-I of the half-wave values)
    w = np.ones(N + 1); w[0] = w[-1] = 0.5
    E = np.array([2.0 / N * np.sum(w * (etan - 1.0) * np.cos(jj * m * np.pi / N))
                  for jj in range(N + 1)])
    E[0] *= 0.5; E[N] *= 0.5
    sg = np.sqrt(g * d)
    out = dict(H=H, T=T, d=d, N=N, current=current,
               k=k / d, c=c * sg, omega=2 * np.pi / T, uE=(c - U) * sg,
               Q=Q * d * sg, R=R * g * d, E=E * d, B=Bn * d * sg,
               crest=(etan[0] - 1) * d, trough=(etan[N] - 1) * d)
    kk, dd, Bd, Ed = out['k'], d, out['B'], out['E']
    jj = np.arange(1, N + 1)

    def eta(xx, t):
        th = kk * np.asarray(xx, float) - out['omega'] * t
        return sum(Ed[i] * np.cos(i * th) for i in range(N + 1))

    def velocity(xx, z, t):
        th = kk * np.asarray(xx, float) - out['omega'] * t
        zb = np.maximum(np.asarray(z, float) + dd, 0)
        u = out['uE'] + sum(i * kk * Bd[i - 1] * np.cosh(i * kk * zb) / np.cosh(i * kk * dd)
                            * np.cos(i * th) for i in jj)
        ww = sum(i * kk * Bd[i - 1] * np.sinh(i * kk * zb) / np.cosh(i * kk * dd)
                 * np.sin(i * th) for i in jj)
        return u, ww

    out['eta'], out['velocity'] = eta, velocity
    return out


def foam_entries(sol, indent='            '):
    """Entries of a forcing zone (or nhSetWaves) for this wave."""
    E = ' '.join('%.10g' % v for v in sol['E'])
    B = ' '.join('%.10g' % v for v in sol['B'])
    return (f"{indent}theory      streamFunction;\n"
            f"{indent}waveHeight  {sol['H']:.10g};\n"
            f"{indent}wavePeriod  {sol['T']:.10g};\n"
            f"{indent}depth       {sol['d']:.10g};\n"
            f"{indent}waveNumber  {sol['k']:.10g};\n"
            f"{indent}uMean       {sol['uE']:.10g};\n"
            f"{indent}etaCoeffs   ({E});\n"
            f"{indent}psiCoeffs   ({B});\n")


if __name__ == '__main__':
    H, T, d = (float(a) for a in sys.argv[1:4])
    N = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    s = fenton(H, T, d, N)
    print('// stream-function wave (Fenton 1988), N = %d: L = %.4f m, c = %.4f m/s,'
          ' crest %.4f m, trough %.4f m, Eulerian mean current %.4f m/s'
          % (N, 2 * np.pi / s['k'], s['c'], s['crest'], s['trough'], s['uE']))
    print(foam_entries(s))
