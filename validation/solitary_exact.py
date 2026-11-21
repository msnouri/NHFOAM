#!/usr/bin/env python3
"""Fully nonlinear solitary wave on water of constant depth (exact solution of
the Euler equations, the same wave as Tanaka's 1986 method gives).

Method: conformal mapping of the fluid layer onto a strip of depth d, in
which the surface elevation Y(alpha) of a steady wave satisfies Babenko's
equation (finite depth, g = d = 1)

    F^2 C(Y) - Y - Y C(Y) - C(Y^2)/2 = 0,     C = k coth(k) in Fourier space,

with F = c/sqrt(g d).  It is solved with Petviashvili's iteration on a long
periodic domain (Clamond & Dutykh 2013, Computers & Fluids 84, 35-38).  The
map z(zeta) = zeta + Z(zeta) gives the physical coordinates, and the complex
velocity in the frame of the wave is u - i w = -c / z_zeta.

usage (python):
    sw = SolitaryWave(a=0.5, d=1.0)     # a: crest height above still water
    sw.c                                # celerity [m/s]
    sw.eta(x)                           # surface elevation, crest at x = 0
    u, w = sw.velocity(x, z)            # fixed frame, wave moving to +x
"""
import math
import numpy as np


class SolitaryWave:
    def __init__(self, a, d=1.0, g=9.81, L=120.0, N=8192, tol=1e-13):
        self.d, self.g = d, g
        self.L, self.N = L, N
        self.alpha = (np.arange(N) - N // 2) * (L / N)
        k = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
        self.k = k
        kc = np.ones_like(k)
        nz = k != 0
        kc[nz] = k[nz] / np.tanh(k[nz])
        self.Csym = kc                          # symbol of C, C(0) = 1
        eps = a / d
        # secant iteration on F for the crest height
        F0 = math.sqrt(1 + eps) * 0.995
        F1 = math.sqrt(1 + eps) * 0.985
        Y0 = self._solve(F0, None, tol)
        Y1 = self._solve(F1, None, tol)
        f0, f1 = Y0.max() - eps, Y1.max() - eps
        for _ in range(60):
            F2 = F1 - f1 * (F1 - F0) / (f1 - f0)
            Y2 = self._solve(F2, Y1, tol)
            F0, f0, F1, f1, Y1 = F1, f1, F2, Y2.max() - eps, Y2
            if abs(f1) < 1e-12:
                break
        self.F = F1
        self.Y = Y1
        self.c = F1 * math.sqrt(g * d)
        self.residual = self._residual(F1, Y1)
        self._surface()

    # -- operators -------------------------------------------------------
    def _C(self, f):
        return np.real(np.fft.ifft(self.Csym * np.fft.fft(f)))

    def _residual(self, F, Y):
        CY = self._C(Y)
        r = F * F * CY - Y - Y * CY - 0.5 * self._C(Y * Y)
        return np.abs(r).max()

    def _solve(self, F, Y, tol):
        """Petviashvili iteration for L Y = N(Y), L = F^2 C - 1."""
        if Y is None:
            A = F * F - 1
            kap = math.sqrt(0.75 * A)
            Y = A / np.cosh(kap * self.alpha) ** 2
        Lsym = F * F * self.Csym - 1
        for _ in range(2000):
            CY = self._C(Y)
            Nl = Y * CY + 0.5 * self._C(Y * Y)
            LY = np.real(np.fft.ifft(Lsym * np.fft.fft(Y)))
            S = np.dot(Y, LY) / np.dot(Y, Nl)
            Yn = S * S * np.real(np.fft.ifft(np.fft.fft(Nl) / Lsym))
            if np.abs(Yn - Y).max() < tol:
                return Yn
            Y = Yn
        return Y

    # -- geometry and kinematics -------------------------------------------
    def _surface(self):
        """Physical surface (X(alpha), Y(alpha)), crest at x = 0."""
        Yh = np.fft.fft(self.Y) / self.N
        self.Yh = Yh
        k = self.k
        nz = k != 0
        # X_alpha = 1 + C(Y): mean 1 + Ybar, periodic part integrated spectrally
        CYh = self.Csym * Yh
        Ph = np.zeros_like(Yh)
        Ph[nz] = CYh[nz] / (1j * k[nz])
        P = np.real(np.fft.ifft(Ph) * self.N)
        Xs = (1 + Yh[0].real) * self.alpha + P
        i0 = np.argmax(self.Y)
        self.x0 = Xs[i0]
        self.Xs = Xs - self.x0

    def _map(self, beta):
        """x, y, z_zeta along the line Im(zeta) = beta (-1 <= beta <= 0),
        nondimensional, for all alpha."""
        k = self.k
        q = np.abs(k)
        nz = k != 0
        Yh = self.Yh
        # coefficients of Z (map) and Z' (map derivative), k != 0
        g = np.zeros_like(k)
        gp = np.zeros_like(k)
        pos = k > 0
        neg = k < 0
        den = 1 - np.exp(-2 * q[nz])
        e = np.zeros_like(k)
        e[pos] = 2 * np.exp(-k[pos] * (beta + 2)) / (1 - np.exp(-2 * k[pos]))
        e[neg] = -2 * np.exp(q[neg] * beta) / (1 - np.exp(-2 * q[neg]))
        # Z_k(zeta) = -i e^{ik(zeta+i)}/sinh(k) = -i e^{ik alpha} * e
        Zh = np.where(nz, -1j * Yh * e, 0)
        Zph = np.where(nz, Yh * k * e, 0)     # dZ/dzeta = k e^{ik(zeta+i)}/sinh(k)
        Z = np.fft.ifft(Zh) * self.N
        Zp = np.fft.ifft(Zph) * self.N
        zeta = self.alpha + 1j * beta
        z = zeta + Yh[0].real * (zeta + 1j) + Z - self.x0
        zz = 1 + Yh[0].real + Zp
        return z.real, z.imag, zz

    def eta(self, x):
        """Surface elevation [m] at x [m] (crest at x = 0)."""
        xs = np.asarray(x) / self.d
        return self.d * np.interp(xs, self.Xs, self.Y)

    def field(self, nBeta=81):
        """Nondimensional x, y, u, w on lines of constant beta
        (arrays nBeta x N); u, w in the fixed frame."""
        betas = np.linspace(-1, 0, nBeta)
        X, Yy, U, W = [], [], [], []
        for b in betas:
            x, y, zz = self._map(b)
            cw = -self.F / zz          # u - i w in the frame of the wave
            X.append(x); Yy.append(y)
            U.append(cw.real + self.F); W.append(-cw.imag)
        return [np.array(v) for v in (X, Yy, U, W)]

    def velocity(self, x, z, nBeta=161):
        """Velocity (u, w) [m/s] in the fixed frame at points (x, z) [m],
        z measured from still water (z = -d at the bed).  Along each line of
        constant beta the map alpha -> x is monotonic: interpolate there to
        the target x, then across the lines to the target z."""
        if not hasattr(self, "_lines") or self._lines[0].shape[0] != nBeta:
            self._lines = self.field(nBeta)
        X, Yy, U, W = self._lines
        xs = np.asarray(x, float).ravel() / self.d
        zs = np.asarray(z, float).ravel() / self.d
        u = np.zeros_like(xs)
        w = np.zeros_like(xs)
        inside = np.abs(xs) < 0.45 * self.L          # far field: at rest
        xi = xs[inside]
        yb = np.array([np.interp(xi, X[b], Yy[b]) for b in range(nBeta)])
        ub = np.array([np.interp(xi, X[b], U[b]) for b in range(nBeta)])
        wb = np.array([np.interp(xi, X[b], W[b]) for b in range(nBeta)])
        zi = zs[inside]
        ui = np.empty_like(xi)
        wi = np.empty_like(xi)
        for j in range(len(xi)):
            ui[j] = np.interp(zi[j], yb[:, j], ub[:, j])
            wi[j] = np.interp(zi[j], yb[:, j], wb[:, j])
        u[inside] = ui
        w[inside] = wi
        s = math.sqrt(self.g * self.d)
        shape = np.shape(x)
        return (u * s).reshape(shape), (w * s).reshape(shape)


if __name__ == "__main__":
    for a in [0.1, 0.3, 0.5, 0.7]:
        sw = SolitaryWave(a)
        print(f"a/d = {a:.2f}: F = {sw.F:.6f}, residual {sw.residual:.1e}, "
              f"F^2 = {sw.F**2:.6f}, KdV sqrt(1+a) = {math.sqrt(1+a):.6f}")
